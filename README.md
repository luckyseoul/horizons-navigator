# JPL Horizons Navigator

<p align="center">
  <img src="logo.png" alt="Horizons orbit navigator" width="168" />
</p>


3D orbit visualization tool for NASA/JPL's Horizons ephemeris system.

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

## Usage

1. Select target (preset or custom like `Apophis;`)
2. Set Center to "Sun (heliocentric)"
3. Set Type to **VECTORS**
4. Click "Get Ephemeris & Plot"
5. Click **▶ Animate** to watch it move

The 3D view draws the Sun and the planetary reference rings only when the
plot origin really is the Sun or the solar-system barycenter. Geocentric,
planet-centred and observatory-centred queries are plotted on their own.

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
