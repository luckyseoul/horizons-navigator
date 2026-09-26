# JPL Horizons Navigator

<p align="center">
  <img src="logo.png" alt="Horizons orbit navigator" width="168" />
</p>


Ephemeris viewer and 3D orbit visualiser for NASA/JPL's Horizons system.

Query any body or designation, plot its trajectory in three dimensions, and
read the reply without wading through Horizons' raw fixed-width text.

## Install (Ubuntu 24.04)
```bash
sudo apt install python3-tk
pip install matplotlib numpy --break-system-packages
python3 horizons_ui.py
```

## Features

- **3D Orbit Visualization** - Interactive plot with rotation/zoom
- **Animation** - Watch objects move along their orbits in real-time
- **Preset Bodies** - Planets, moons, barycenters
- **Custom Queries** - Asteroids, comets by designation
- **Reference Orbits** - Mercury/Venus/Earth/Mars shown for scale
- **Observer Quantities** - Pick Horizons' quantity codes by name
- **Reply Summary** - Target, center, time span, step and sample count up front
- **CSV Export** - Save the data block as real CSV, not a text dump

## Usage

1. Select target (preset or custom like `Apophis;`)
2. Set Center to "Sun (heliocentric)"
3. Set Type to **VECTORS**
4. Click "Get Ephemeris & Plot"
5. Click **▶ Animate** to watch it move

The 3D view draws the Sun and the planetary reference rings only when the
plot origin really is the Sun or the solar-system barycenter. Geocentric,
planet-centred and observatory-centred queries are plotted on their own.

Every reply opens with a summary — target, center, time span, step size and
sample count — above the raw Horizons text.

**Quantities** apply only to OBSERVER ephemerides. The dropdown beside the
field lists Horizons' quantity codes by name; picking one appends its code
and selects the OBSERVER type for you.

**Save** writes the displayed text. Give the file a `.csv` extension and it
writes the data block as real CSV instead — for a VECTORS reply that is one
`x_au, y_au, z_au, distance_au` row per timestamp. Replies with no data
block fall back to text with a warning.

## Tests

The suite runs offline — parsing, CENTER classification, COMMAND encoding,
validation and the Tk application are all exercised without touching the
network.

```bash
python3 -m unittest test_horizons_ui -v
```

The GUI tests need a display; headless machines can use Xvfb:

```bash
xvfb-run -a python3 -m unittest test_horizons_ui -v
```

One test talks to JPL and is opt-in:

```bash
HORIZONS_LIVE=1 python3 -m unittest test_horizons_ui -v
```

## API

https://ssd-api.jpl.nasa.gov/doc/horizons.html
