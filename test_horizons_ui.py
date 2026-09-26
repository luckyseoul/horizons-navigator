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

import os
import threading
import time
import unittest

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

# A trimmed OBSERVER reply in the CSV layout used for link timing.
LINK_CSV = """\
***********************************************************************************************
 Date__(UT)__HR:MN, , ,            delta,     deldot,  1-way_down_LT,       VmagSn,     VmagOb,
***********************************************************************************************
$$SOE
 2026-Sep-26 00:00, , , 1.69933032252892,-11.4546799,    14.13289934,   29.6901475, 32.5118381,
 2026-Sep-28 00:00, , , 1.68600676490452,-11.6126242,   14.02209069,   29.7071702, 32.1666973,
$$EOE
***********************************************************************************************
"""

# The same row with the columns shuffled, to prove the positions come from
# the header rather than from a hard-coded layout.
LINK_CSV_REORDERED = """\
 Date__(UT)__HR:MN,     deldot,       VmagSn,  1-way_down_LT,            delta,     VmagOb,
$$SOE
 2026-Sep-26 00:00,-11.4546799,   29.6901475,    14.13289934, 1.69933032252892, 32.5118381,
$$EOE
"""

# Constants Horizons itself reports in the reply header.
AU_KM = 149597870.700
C_KMS = 299792.458


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


class TestLinkTiming(unittest.TestCase):
    def test_rows_are_parsed(self):
        rows = HorizonsUI._parse_link_times(LINK_CSV)
        self.assertEqual(len(rows), 2)
        date, range_au, light_min = rows[0]
        self.assertEqual(date, "2026-Sep-26 00:00")
        self.assertAlmostEqual(range_au, 1.69933032252892)
        self.assertAlmostEqual(light_min, 14.13289934)

    def test_light_time_equals_range_over_c(self):
        """The reported light time must be range / c - a real physics check."""
        for _, range_au, light_min in HorizonsUI._parse_link_times(LINK_CSV):
            expected = range_au * AU_KM / C_KMS / 60.0
            self.assertAlmostEqual(light_min, expected, delta=0.001)

    def test_columns_are_found_by_header_not_position(self):
        rows = HorizonsUI._parse_link_times(LINK_CSV_REORDERED)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0][1], 1.69933032252892)   # delta
        self.assertAlmostEqual(rows[0][2], 14.13289934)        # 1-way_down_LT

    def test_reply_without_link_columns(self):
        self.assertIsNone(HorizonsUI._parse_link_times("no header here"))
        self.assertIsNone(HorizonsUI._parse_link_times(LABELLED))


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

    def test_link_tab_and_chart(self):
        tabs = [self.app.view_tabs.tab(i, "text")
                for i in range(self.app.view_tabs.index("end"))]
        self.assertEqual(tabs, ["3D Orbit", "DTN Link Timing"])

    def test_link_summary_renders_and_charts(self):
        self.app._show_link_summary("Earth", "Mars", {"result": LINK_CSV})
        text = self.app.output_text.get(1.0, tk.END)
        self.assertIn("DTN link timing: Earth <-> Mars", text)
        self.assertIn("One-way light time", text)
        self.assertIn("2026-Sep-26 00:00", text)
        self.assertEqual(self.app.view_tabs.select(), str(self.app.link_tab))
        self.assertEqual(len(self.app.link_plotter.ax.get_lines()), 1)
        self.assertIn("Earth <-> Mars", self.app.link_plotter.ax.get_title())

    def test_link_summary_handles_a_reply_without_timing(self):
        self.app._show_link_summary("Earth", "Mars", {"result": "nope"})
        self.assertEqual(self.app.status_var.get(),
                         "No range/light-time columns in the reply")


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

    def test_earth_mars_light_time_round_trip(self):
        api = HorizonsAPI()
        result = api.query({
            "COMMAND": "399", "CENTER": "500@499",
            "START_TIME": "2026-09-26", "STOP_TIME": "2026-12-26",
            "STEP_SIZE": "10 d", "EPHEM_TYPE": "OBSERVER",
            "QUANTITIES": "20,21,22", "CSV_FORMAT": "YES",
            "OBJ_DATA": "NO", "MAKE_EPHEM": "YES"})
        rows = HorizonsUI._parse_link_times(result["result"])
        self.assertGreater(len(rows), 5)
        for _, range_au, light_min in rows:
            self.assertAlmostEqual(light_min, range_au * AU_KM / C_KMS / 60.0,
                                   delta=0.01)
        mean = float(np.mean([row[2] for row in rows]))
        self.assertTrue(3.0 < mean < 25.0,
                        f"Earth-Mars mean light time {mean:.2f} min")


if __name__ == "__main__":
    unittest.main(verbosity=2)
