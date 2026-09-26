# Progress — horizons-navigator

## 2026-09-26

Verified the tool end-to-end against the live JPL service, then fixed the
defects that verification turned up. The DTN framing turned out to be
buildable rather than aspirational, so it now carries a real tool.

### Found and fixed

- **Orbit geometry was wrong.** Horizons prints a position row (`X= Y= Z=`)
  and a velocity row (`VX= VY= VZ=`) under every timestamp, and the parser
  was reading both. A 74-step Mars query plotted 148 points, so the drawn
  polyline jumped up to 1.65 AU between vertices where the real per-step
  motion is at most 0.07 AU. Position rows are now anchored to the start of
  the line, so velocity components can never be mistaken for positions.
- **The Sun and the planetary reference rings were drawn for every center.**
  The heliocentric test was a substring check (`"10" in center or "0" in
  center`), which is true for every `500@…` code — including geocentric and
  Mars-centred queries. CENTER codes are now classified exactly, and the
  title says which origin was used.
- **Stopping the animation stacked duplicate orbit lines.** The static orbit
  was re-plotted on every stop; three start/stop cycles left eight lines on
  the axes instead of five.
- **`query()` deleted `COMMAND` from the caller's dict** (it popped in place).
- Plot labels followed the preset combo rather than the command actually
  entered, so a custom designation was still labelled with the previous
  preset.
- Tk is no longer touched from worker threads; background results are handed
  back through a queue drained on the main thread.
- `datetime.utcnow()` (deprecated) replaced with an aware UTC timestamp,
  COMMAND values are fully percent-encoded, saved output is written as UTF-8,
  and missing time fields now produce a named message instead of a round trip
  to JPL.

### Added

- **DTN link timing.** A second tab computing expected bundle travel time
  between two nodes. It asks Horizons for the range and the one-way down-leg
  light time (`QUANTITIES=20,21,22`, CSV), charts light time across the query
  window with the best and worst samples marked, and reports round-trip time
  and range. Columns are located by CSV header name, so the parse does not
  depend on the order the quantities were requested in.
- `test_horizons_ui.py` — 35 tests. Offline by default (parsing, CENTER
  classification, encoding, validation, link timing, and the Tk application
  under Xvfb); two opt-in live queries against JPL via `HORIZONS_LIVE=1`. All
  of the regression tests fail against the previous revision, so they are
  known to bite.

### Verified

- Full 1-year Mars query against the live service: 74 points, heliocentric
  distance 1.494–1.666 AU (Mars is 1.381–1.666 AU), title
  "Mars Orbit (Heliocentric)", animation starts and stops without disturbing
  the plot.
- Earth↔Mars link over the same window: one-way light time 5.64–17.09 min
  (round trip 11.28–34.18 min), range 0.678–2.055 AU. The shortest travel
  time lands on 2027-Feb-18, which is the February 2027 Mars opposition —
  the curve finds it unaided. The reported light time matches range / c to
  better than 1e-6 minutes at every sample, which checks the parse and the
  physics at once.

## 2026-07-14

- 3D orbit visualization for NASA/JPL Horizons ephemerides remains the public tool baseline.
- Progress checkpoint: useful for deep-space / Solar System Internet routing narratives and mission geometry intuition alongside DTN work.
- Next: optional ephemeris query polish and screenshot assets for presentations when needed.

_Updated 2026-09-26 09:26 UTC_
