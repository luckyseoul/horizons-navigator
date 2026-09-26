#!/usr/bin/env python3
"""Tests for horizons_ui.py.

Everything here is offline by default: vector parsing, CENTER
classification, COMMAND encoding, input validation and - when a display is
available - the Tk application itself.

    python3 -m unittest test_horizons_ui -v

The GUI tests need a display; on a headless machine run them under Xvfb:

    xvfb-run -a python3 -m unittest test_horizons_ui -v

The single test that talks to JPL is opt-in:

    HORIZONS_LIVE=1 python3 -m unittest test_horizons_ui -v
"""

import csv
import os
import tempfile
import threading
import time
import unittest
from unittest import mock

import tkinter as tk

import numpy as np

import horizons_ui
from horizons_ui import HorizonsAPI, HorizonsUI


# A trimmed copy of a real VECTORS reply: every timestamp carries a position
# row followed by a velocity row.
LABELLED = """\
*******************************************************************************
Ephemeris / API_USER ...
*******************************************************************************
Target body name: Mars (499)                      {source: mar099}
Center body name: Sun (10)                        {source: DE441}
Start time      : A.D. 2026-Sep-26 00:00:00.0000 TDB
Stop  time      : A.D. 2026-Oct-01 00:00:00.0000 TDB
Step-size       : 7200 minutes
*******************************************************************************
$$SOE
2461309.500000000 = A.D. 2026-Sep-26 00:00:00.0000 TDB
 X = 1.941968158932753E-01 Y = 1.539603497469368E+00 Z = 2.750271703788334E-02
 VX=-1.335339041450438E-02 VY= 2.939882132815847E-03 VZ= 3.890369976799800E-04
2461314.500000000 = A.D. 2026-Oct-01 00:00:00.0000 TDB
 X = 1.272603245715422E-01 Y = 1.552781205446242E+00 Z = 2.942016831701569E-02
 VX=-1.341669120657371E-02 VY= 2.331553317528713E-03 VZ= 3.778401905247213E-04
$$EOE
*******************************************************************************
"""

# The same two timestamps with CSV_FORMAT=YES.
CSV = """\
$$SOE
2461309.500000000, A.D. 2026-Sep-26 00:00:00.0000, 1.941968158932753E-01, 1.539603497469368E+00, 2.750271703788334E-02, -1.335339041450438E-02, 2.939882132815847E-03, 3.890369976799800E-04,
2461314.500000000, A.D. 2026-Oct-01 00:00:00.0000, 1.272603245715422E-01, 1.552781205446242E+00, 2.942016831701569E-02, -1.341669120657371E-02, 2.331553317528713E-03, 3.778401905247213E-04,
$$EOE
"""

# A comma-formatted data block, as CSV_FORMAT=YES produces.
CSV_BLOCK = """\
 Date__(UT)__HR:MN, , ,            delta,     deldot,  1-way_down_LT,
$$SOE
 2026-Sep-26 00:00, , , 1.69933032252892,-11.4546799,    14.13289934,
$$EOE
"""


class TestVectorParsing(unittest.TestCase):
    def test_one_point_per_timestamp(self):
        """Regression: a reply with 2 timestamps must plot exactly 2 points."""
        x, y, z = HorizonsUI._parse_vectors(LABELLED)
        self.assertEqual(len(x), 2)
        self.assertEqual(len(y), 2)
        self.assertEqual(len(z), 2)

    def test_velocity_components_are_not_positions(self):
        """Regression: VX/VY/VZ must never be read as X/Y/Z."""
        x, _, _ = HorizonsUI._parse_vectors(LABELLED)
        self.assertNotIn(-1.335339041450438E-02, x)  # VX of the first timestamp
        self.assertNotIn(-1.341669120657371E-02, x)  # VX of the second

    def test_position_values(self):
        x, y, z = HorizonsUI._parse_vectors(LABELLED)
        np.testing.assert_allclose(
            x, [1.941968158932753E-01, 1.272603245715422E-01])
        np.testing.assert_allclose(
            z, [2.750271703788334E-02, 2.942016831701569E-02])

    def test_csv_layout(self):
        x, y, z = HorizonsUI._parse_vectors(CSV)
        self.assertEqual(len(x), 2)
        np.testing.assert_allclose(
            x, [1.941968158932753E-01, 1.272603245715422E-01])
        np.testing.assert_allclose(
            y, [1.539603497469368E+00, 1.552781205446242E+00])

    def test_compact_negative_notation(self):
        """Horizons writes 'X =-7.2E-03' with no space after the sign."""
        text = ("$$SOE\n"
                " X =-7.227630149399285E-03 Y = 1.570046147944983E+00"
                " Z = 3.307963746710238E-02\n"
                "$$EOE")
        x, y, z = HorizonsUI._parse_vectors(text)
        self.assertEqual(len(x), 1)
        np.testing.assert_allclose(x, [-7.227630149399285E-03])
        np.testing.assert_allclose(z, [3.307963746710238E-02])

    def test_missing_block_returns_none(self):
        self.assertIsNone(HorizonsUI._parse_vectors("no data here"))
        self.assertIsNone(HorizonsUI._parse_vectors("$$SOE\n$$EOE"))


class TestCenterOrigin(unittest.TestCase):
    def test_sun_is_heliocentric(self):
        self.assertEqual(HorizonsAPI.center_origin("500@10"), "Heliocentric")

    def test_barycenter_is_barycentric(self):
        self.assertEqual(HorizonsAPI.center_origin("500@0"), "Barycentric")

    def test_geocentric_is_not_solar(self):
        """Regression: '500@399' contains a '0', which used to be enough."""
        self.assertIsNone(HorizonsAPI.center_origin("500@399"))

    def test_planet_centred_is_not_solar(self):
        for code in ("500@499", "500@599", "500@699"):
            self.assertIsNone(HorizonsAPI.center_origin(code))

    def test_observatory_is_not_solar(self):
        self.assertIsNone(HorizonsAPI.center_origin("568"))

    def test_blank_input(self):
        self.assertIsNone(HorizonsAPI.center_origin(""))
        self.assertIsNone(HorizonsAPI.center_origin(None))


class TestCommandEncoding(unittest.TestCase):
    def test_small_body_semicolon(self):
        self.assertEqual(HorizonsAPI.encode_command("Apophis;"), "Apophis%3B")

    def test_designation_lookup(self):
        self.assertEqual(
            HorizonsAPI.encode_command("DES=2024 YR4;"), "DES%3D2024%20YR4%3B")

    def test_plain_code_is_unchanged(self):
        self.assertEqual(HorizonsAPI.encode_command("499"), "499")


class TestQueryDoesNotMutateParams(unittest.TestCase):
    def test_command_survives_a_failed_query(self):
        """Regression: query() used to pop COMMAND out of the caller's dict."""
        class Unreachable(HorizonsAPI):
            BASE_URL = "http://127.0.0.1:9/api/horizons.api"

        params = {"COMMAND": "499", "CENTER": "500@10"}
        with self.assertRaises(Exception):
            Unreachable().query(params)
        self.assertEqual(params, {"COMMAND": "499", "CENTER": "500@10"})


class TestValidation(unittest.TestCase):
    def test_valid_ephemeris_request(self):
        self.assertIsNone(HorizonsUI._validate_inputs({
            "COMMAND": "499", "CENTER": "500@10", "MAKE_EPHEM": "YES",
            "START_TIME": "2026-01-01", "STOP_TIME": "2026-02-01",
            "STEP_SIZE": "1 d"}))

    def test_missing_target(self):
        msg = HorizonsUI._validate_inputs({"COMMAND": "", "CENTER": "500@10"})
        self.assertIn("target body", msg)

    def test_missing_time_fields_are_named(self):
        msg = HorizonsUI._validate_inputs({
            "COMMAND": "499", "CENTER": "500@10", "MAKE_EPHEM": "YES",
            "START_TIME": "", "STOP_TIME": "2026-02-01", "STEP_SIZE": ""})
        self.assertIn("START_TIME", msg)
        self.assertIn("STEP_SIZE", msg)

    def test_object_lookup_needs_no_time_span(self):
        self.assertIsNone(HorizonsUI._validate_inputs(
            {"COMMAND": "Apophis;", "CENTER": "500@10", "MAKE_EPHEM": "NO"}))


class TestReplySummary(unittest.TestCase):
    def test_header_fields_are_extracted(self):
        joined = "\n".join(HorizonsUI._summarize_reply(LABELLED))
        self.assertIn("Target  : Mars (499)", joined)
        self.assertIn("Center  : Sun (10)", joined)
        self.assertIn("Start   : A.D. 2026-Sep-26", joined)
        self.assertIn("Stop    : A.D. 2026-Oct-01", joined)
        self.assertIn("Step    : 7200 minutes", joined)

    def test_sample_count_counts_timestamps_not_printed_lines(self):
        """Regression: 2 timestamps must not be counted as 6 printed rows.

        A labelled VECTORS sample spans a date row, a position row and a
        velocity row, so counting lines over-reports by 3x.
        """
        self.assertIn("Samples : 2", HorizonsUI._summarize_reply(LABELLED))
        self.assertEqual(len(HorizonsUI._parse_vectors(LABELLED)[0]), 2)

    def test_whitespace_in_values_is_collapsed(self):
        joined = "\n".join(HorizonsUI._summarize_reply(LABELLED))
        self.assertIn("Mars (499) {source: mar099}", joined)

    def test_plain_text_has_no_summary(self):
        self.assertIsNone(HorizonsUI._summarize_reply("nothing to see here"))


class TestCsvRows(unittest.TestCase):
    def test_labelled_vectors_become_position_rows(self):
        rows = HorizonsUI._csv_rows(LABELLED)
        self.assertEqual(rows[0], ["x_au", "y_au", "z_au", "distance_au"])
        self.assertEqual(len(rows), 3)  # header + 2 timestamps

    def test_distance_column_is_consistent(self):
        for x, y, z, distance in HorizonsUI._csv_rows(LABELLED)[1:]:
            expected = (float(x) ** 2 + float(y) ** 2 + float(z) ** 2) ** 0.5
            self.assertAlmostEqual(float(distance), expected, places=9)

    def test_velocity_rows_are_not_exported_as_positions(self):
        rows = HorizonsUI._csv_rows(LABELLED)
        self.assertEqual(len(rows), 3)
        self.assertNotAlmostEqual(float(rows[1][0]), -1.335339041450438E-02)

    def test_comma_reply_passes_through(self):
        rows = HorizonsUI._csv_rows(CSV_BLOCK)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0].strip(), "2026-Sep-26 00:00")
        self.assertEqual(rows[0][3].strip(), "1.69933032252892")

    def test_no_block(self):
        self.assertIsNone(HorizonsUI._csv_rows("prose only"))
        self.assertIsNone(HorizonsUI._csv_rows(""))


class TestApplication(unittest.TestCase):
    """Drives the real Tk application.  Skipped when there is no display."""

    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError as exc:  # no DISPLAY / no X server
            raise unittest.SkipTest(f"no display available: {exc}")
        cls.root.withdraw()
        cls.app = HorizonsUI(cls.root)

    @classmethod
    def tearDownClass(cls):
        try:
            cls.root.destroy()
        except Exception:
            pass

    def test_widgets_exist(self):
        self.assertIsNotNone(self.app.plotter)
        self.assertIsNotNone(self.app.output_text)
        self.assertIsNotNone(self.app.animate_btn)

    def test_target_label_follows_the_command_code(self):
        original = self.app.command_var.get()
        try:
            self.app.preset_var.set("Mars")
            self.app.command_var.set("599")
            self.assertEqual(self.app._target_label(), "Jupiter")
            self.app.command_var.set("Apophis;")
            self.assertEqual(self.app._target_label(), "Apophis")
        finally:
            self.app.command_var.set(original)

    def test_plot_consumes_each_timestamp_once(self):
        self.app.command_var.set("499")
        self.app._plot_result({"result": LABELLED},
                              {"EPHEM_TYPE": "VECTORS", "CENTER": "500@10"})
        x, y, z = self.app.last_vectors
        self.assertEqual(len(x), 2)
        self.assertEqual(self.app.status_var.get(), "Plotted 2 points")

    def test_solar_reference_rings_only_for_solar_centers(self):
        self.app.command_var.set("499")
        self.app._plot_result({"result": LABELLED},
                              {"EPHEM_TYPE": "VECTORS", "CENTER": "500@399"})
        self.assertEqual(self.app.plotter.orbits["Mars"][0].shape, (2,))
        self.assertNotIn("Sun", self.app.plotter.orbits)
        self.assertIn("center 500@399", self.app.plotter.ax.get_title())

    def test_stop_animation_does_not_stack_artists(self):
        theta = np.linspace(0, 2 * np.pi, 32)
        self.app.last_vectors = (np.cos(theta), np.sin(theta),
                                 np.zeros_like(theta))
        before = len(self.app.plotter.ax.get_lines())
        for _ in range(3):
            self.app._toggle_animation()
            self.app._stop_animation()
        self.assertEqual(len(self.app.plotter.ax.get_lines()), before)

    def test_worker_results_are_delivered_on_the_main_thread(self):
        seen = {}
        worker = threading.Thread(
            target=lambda: self.app._post(
                lambda: seen.setdefault("thread", threading.current_thread())))
        worker.start()
        worker.join()
        deadline = time.time() + 5
        while "thread" not in seen and time.time() < deadline:
            self.root.update()
            time.sleep(0.02)
        self.assertEqual(seen.get("thread"), threading.main_thread())

    def test_display_result_puts_summary_above_raw_text(self):
        self.app._display_result({"result": LABELLED})
        text = self.app.output_text.get(1.0, tk.END)
        self.assertIn("Target  : Mars (499)", text)
        self.assertIn("Samples : 2", text)
        self.assertIn("$$SOE", text)
        self.assertLess(text.index("Target  :"), text.index("$$SOE"))
        self.assertEqual(self.app.raw_reply, LABELLED)

    def test_quantity_picker_appends_and_switches_to_observer(self):
        self.app.quantities_entry_var.set("1")
        self.app.ephem_type_var.set("VECTORS")
        self.app.quantity_picker.set("9  Visual mag. & Surf Brt")
        self.app._on_quantity_pick()
        self.assertEqual(self.app.quantities_entry_var.get(), "1,9")
        self.assertEqual(self.app.ephem_type_var.get(), "OBSERVER")

    def test_quantity_picker_does_not_duplicate(self):
        self.app.quantities_entry_var.set("1")
        self.app.quantity_picker.set("1  Astrometric RA & DEC")
        self.app._on_quantity_pick()
        self.app._on_quantity_pick()
        self.assertEqual(self.app.quantities_entry_var.get(), "1")

    def test_save_as_csv_writes_real_csv(self):
        self.app._display_result({"result": LABELLED})
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "orbit.csv")
            with mock.patch.object(horizons_ui.filedialog, "asksaveasfilename",
                                   return_value=path):
                self.app._save_results()
            self.assertTrue(os.path.exists(path))
            with open(path, newline="") as handle:
                rows = list(csv.reader(handle))
        self.assertEqual(rows[0], ["x_au", "y_au", "z_au", "distance_au"])
        self.assertEqual(len(rows), 3)

    def test_save_as_csv_falls_back_to_text_when_not_tabular(self):
        self.app._display_result({"result": "prose with no data block"})
        warned = []
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "note.csv")
            with mock.patch.object(horizons_ui.filedialog, "asksaveasfilename",
                                   return_value=path), \
                 mock.patch.object(horizons_ui.messagebox, "showwarning",
                                   side_effect=lambda *a, **k: warned.append(a)):
                self.app._save_results()
            self.assertTrue(os.path.exists(path))
        self.assertEqual(len(warned), 1)
        self.assertIn("Saved", self.app.status_var.get())


@unittest.skipUnless(os.environ.get("HORIZONS_LIVE") == "1",
                     "set HORIZONS_LIVE=1 to query the JPL service")
class TestLiveAPI(unittest.TestCase):
    def test_mars_state_vectors_round_trip(self):
        api = HorizonsAPI()
        result = api.query({
            "COMMAND": "499", "CENTER": "500@10",
            "START_TIME": "2026-09-26", "STOP_TIME": "2026-12-26",
            "STEP_SIZE": "5 d", "EPHEM_TYPE": "VECTORS", "OBJ_DATA": "YES",
            "MAKE_EPHEM": "YES", "CSV_FORMAT": "NO", "VEC_TABLE": "2",
            "REF_PLANE": "ECLIPTIC", "OUT_UNITS": "AU-D"})
        x, y, z = HorizonsUI._parse_vectors(result["result"])
        self.assertGreater(len(x), 10)
        distance = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        self.assertTrue(1.3 < float(np.mean(distance)) < 1.8,
                        f"Mars mean heliocentric distance {np.mean(distance):.3f} AU")


if __name__ == "__main__":
    unittest.main(verbosity=2)
