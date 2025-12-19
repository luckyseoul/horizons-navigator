# JPL Horizons Navigator

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

## API

https://ssd-api.jpl.nasa.gov/doc/horizons.html
