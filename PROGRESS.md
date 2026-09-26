# Progress — horizons-navigator

## 2026-09-26

Verified the tool end-to-end against the live JPL service, then fixed the
defects that verification turned up.

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

- `test_horizons_ui.py` — 41 tests. Offline by default (parsing, CENTER
  classification, encoding, validation, reply summarising, CSV conversion and
  the Tk application under Xvfb); two opt-in live queries against JPL via
  `HORIZONS_LIVE=1`. All of the regression tests fail against the revision
  before this work, so they are known to bite.

### Improved

- **Observer quantities are pickable.** The 33-entry `HorizonsAPI.QUANTITIES`
  map had never been wired to anything, so the field was free text and you had
  to know the codes. A dropdown now lists them by name; picking one appends
  its code and selects OBSERVER, which is the only type the field applies to.
- **Replies open with a summary.** Target, center, time span, step size and
  sample count, above the raw text. The sample count counts timestamps rather
  than printed lines, so a labelled VECTORS reply with 74 timestamps reports
  74 and not 222.
- **Save honours a `.csv` filename.** It previously wrote the raw text whatever
  the extension claimed. It now writes a real CSV data block — position and
  distance per timestamp for a VECTORS reply — and says so when a reply has
  nothing tabular to write.

### Removed

- DTN link timing. It was added, worked, and was taken back out so the tool
  stays an ephemeris viewer rather than a link-budgeting one. The revert is a
  single commit and the orbit fixes underneath were left untouched.

### Verified

- Full 1-year Mars query against the live service: 74 points, heliocentric
  distance 1.494–1.666 AU (Mars is 1.381–1.666 AU), title
  "Mars Orbit (Heliocentric)", animation starts and stops without disturbing
  the plot.
- The whole app driven headlessly against the live service after these
  changes: summary header rendered, 74-point orbit at 1.494–1.666 AU, a
  75-row CSV export (`x_au, y_au, z_au, distance_au`), the quantity picker
  appending its code and switching type to OBSERVER, and three animate/stop
  cycles leaving the line count at 5.

## 2026-07-14

- 3D orbit visualization for NASA/JPL Horizons ephemerides remains the public tool baseline.
- Progress checkpoint: useful for mission geometry intuition and ephemeris exploration.
- Next: optional ephemeris query polish and screenshot assets for presentations when needed.

_Updated 2026-09-26 09:34 UTC_
