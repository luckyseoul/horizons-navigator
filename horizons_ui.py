#!/usr/bin/env python3
"""
JPL Horizons API Navigator with 3D Orbit Visualization
A tkinter UI for querying NASA/JPL's Horizons ephemeris system.

Author: Nicholas Perry
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import urllib.request
import urllib.parse
import json
import re
import math
from datetime import datetime, timedelta
from threading import Thread

# Matplotlib imports for 3D plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


class HorizonsAPI:
    """Interface to JPL Horizons API"""
    
    BASE_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    # Common body codes
    MAJOR_BODIES = {
        "Sun": "10",
        "Mercury": "199",
        "Venus": "299",
        "Earth": "399",
        "Moon": "301",
        "Mars": "499",
        "Phobos": "401",
        "Deimos": "402",
        "Jupiter": "599",
        "Io": "501",
        "Europa": "502",
        "Ganymede": "503",
        "Callisto": "504",
        "Saturn": "699",
        "Titan": "606",
        "Enceladus": "602",
        "Uranus": "799",
        "Neptune": "899",
        "Triton": "801",
        "Pluto": "999",
        "Charon": "901",
        # Barycenters
        "Solar System Barycenter": "0",
        "Mercury Barycenter": "1",
        "Venus Barycenter": "2",
        "Earth-Moon Barycenter": "3",
        "Mars Barycenter": "4",
        "Jupiter Barycenter": "5",
        "Saturn Barycenter": "6",
        "Uranus Barycenter": "7",
        "Neptune Barycenter": "8",
        "Pluto Barycenter": "9",
    }
    
    # Semi-major axes in AU for reference orbits
    PLANET_ORBITS = {
        "Mercury": 0.387,
        "Venus": 0.723,
        "Earth": 1.000,
        "Mars": 1.524,
        "Jupiter": 5.203,
        "Saturn": 9.537,
        "Uranus": 19.19,
        "Neptune": 30.07,
    }
    
    # Planet colors for visualization
    PLANET_COLORS = {
        "Sun": "#FFD700",
        "Mercury": "#B5B5B5",
        "Venus": "#E6C229",
        "Earth": "#6B93D6",
        "Mars": "#C1440E",
        "Jupiter": "#D8CA9D",
        "Saturn": "#F4D59E",
        "Uranus": "#D1E7E7",
        "Neptune": "#5B5DDF",
        "Pluto": "#9C8A7C",
    }
    
    # Observer quantity codes
    QUANTITIES = {
        "1": "Astrometric RA & DEC",
        "2": "Apparent RA & DEC",
        "3": "Rates; RA & DEC",
        "4": "Apparent AZ & EL",
        "7": "Local apparent SOLAR TIME",
        "8": "Airmass & extinction",
        "9": "Visual mag. & Surf Brt",
        "10": "Illuminated fraction",
        "12": "Satellite angle & sep",
        "13": "Target angular diam.",
        "14": "Observer sub-lon & sub-lat",
        "15": "Sun sub-lon & sub-lat",
        "17": "North Pole pos. angle",
        "18": "Helio ecliptic lon & lat",
        "19": "Helio range & range-rate",
        "20": "Observer range & range-rate",
        "21": "One-way down-leg light-time",
        "22": "Speed wrt Sun & observer",
        "23": "Sun-Target-Observer angle",
        "24": "Sun-Observer-Target angle",
        "25": "Target-Observer-Moon angle",
        "29": "Constellation ID",
        "31": "Observer ecliptic lon & lat",
        "32": "North Pole RA & DEC",
        "33": "Galactic latitude",
        "36": "Obs sub-lng & sub-lat (Planet)",
        "37": "Sun sub-lng & sub-lat (Planet)",
        "38": "Phase angle & bisector",
        "39": "Helio range & rng-rate (apparent)",
        "40": "Obs range & rng-rate (apparent)",
        "42": "Phase angle & angle defect",
        "43": "Target-Observer-Earth angle",
        "A": "All available quantities",
    }
    
    CENTERS = {
        "Sun (heliocentric)": "500@10",
        "Solar System Barycenter": "500@0",
        "Geocentric": "500@399",
        "Mars": "500@499",
        "Jupiter": "500@599",
        "Saturn": "500@699",
        # Common observatories
        "Mauna Kea": "568",
        "Paranal": "309",
        "La Silla": "809",
        "Cerro Tololo": "807",
        "Kitt Peak": "695",
        "Palomar Mountain": "675",
        "Arecibo": "251",
        "Goldstone DSN": "-14",
        "Canberra DSN": "-49",
        "Madrid DSN": "-55",
    }
    
    @staticmethod
    def encode_command(cmd: str) -> str:
        """URL-encode the COMMAND parameter"""
        cmd = cmd.replace(";", "%3B")
        cmd = cmd.replace("=", "%3D")
        return cmd
    
    def query(self, params: dict) -> dict:
        """Execute API query and return results"""
        
        # Start with format (no quotes on this one)
        query_parts = ["format=json"]
        
        # Handle COMMAND encoding specially
        if "COMMAND" in params:
            cmd = params.pop("COMMAND")
            cmd_encoded = self.encode_command(cmd)
            query_parts.append(f"COMMAND='{cmd_encoded}'")
        
        # Remove format if accidentally in params
        params.pop("format", None)
        
        # Build query string - properly encode all values
        for key, value in params.items():
            if value:
                encoded_value = urllib.parse.quote(str(value), safe="")
                query_parts.append(f"{key}='{encoded_value}'")
        
        url = f"{self.BASE_URL}?{'&'.join(query_parts)}"
        
        print(f"DEBUG URL: {url}")
        
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "HorizonsUI/1.0")
            
            with urllib.request.urlopen(req, timeout=60) as response:
                data = response.read().decode()
                return json.loads(data)
        except urllib.error.HTTPError as e:
            error_body = e.read().decode() if e.fp else ""
            raise Exception(f"HTTP {e.code}: {e.reason}\n{error_body}")
        except urllib.error.URLError as e:
            raise Exception(f"Connection failed: {e.reason}")
        except json.JSONDecodeError as e:
            raise Exception(f"Invalid JSON response: {e}")


class OrbitPlotter:
    """3D orbit visualization with animation"""
    
    def __init__(self, parent_frame):
        self.frame = parent_frame
        
        # Create figure with dark background
        self.fig = plt.Figure(figsize=(6, 6), dpi=100, facecolor='#1e1e1e')
        self.ax = self.fig.add_subplot(111, projection='3d', facecolor='#1e1e1e')
        
        # Style the 3D axes
        self._style_axes()
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar
        toolbar_frame = ttk.Frame(parent_frame)
        toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()
        
        # Store orbit data
        self.orbits = {}
        
        # Animation state
        self.animation = None
        self.is_animating = False
        self.anim_marker = None
        self.anim_trail = None
        self.anim_frame = 0
        self.anim_speed = 1
        self.trail_length = 20  # Number of points in trail
        
    def _style_axes(self):
        """Apply dark theme to 3D axes"""
        self.ax.set_facecolor('#1e1e1e')
        
        # Axis colors
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        
        self.ax.xaxis.pane.set_edgecolor('#444444')
        self.ax.yaxis.pane.set_edgecolor('#444444')
        self.ax.zaxis.pane.set_edgecolor('#444444')
        
        # Grid
        self.ax.xaxis._axinfo["grid"]['color'] = '#333333'
        self.ax.yaxis._axinfo["grid"]['color'] = '#333333'
        self.ax.zaxis._axinfo["grid"]['color'] = '#333333'
        
        # Labels
        self.ax.set_xlabel('X (AU)', color='#aaaaaa')
        self.ax.set_ylabel('Y (AU)', color='#aaaaaa')
        self.ax.set_zlabel('Z (AU)', color='#aaaaaa')
        
        # Tick colors
        self.ax.tick_params(colors='#888888')
        
    def clear(self):
        """Clear the plot"""
        self.ax.cla()
        self._style_axes()
        self.orbits = {}
        self.canvas.draw()
        
    def plot_sun(self):
        """Plot the Sun at origin"""
        self.ax.scatter([0], [0], [0], c='#FFD700', s=200, marker='o', label='Sun')
        
    def plot_reference_orbits(self, planets=None):
        """Plot circular reference orbits for planets"""
        if planets is None:
            planets = ["Mercury", "Venus", "Earth", "Mars"]
        
        theta = np.linspace(0, 2*np.pi, 100)
        
        for planet in planets:
            if planet in HorizonsAPI.PLANET_ORBITS:
                a = HorizonsAPI.PLANET_ORBITS[planet]
                x = a * np.cos(theta)
                y = a * np.sin(theta)
                z = np.zeros_like(theta)
                
                color = HorizonsAPI.PLANET_COLORS.get(planet, '#666666')
                self.ax.plot(x, y, z, '--', color=color, alpha=0.3, linewidth=1)
    
    def plot_orbit(self, x, y, z, label="Target", color='#00ff00', linewidth=2):
        """Plot an orbit trajectory"""
        self.ax.plot(x, y, z, color=color, linewidth=linewidth, label=label)
        
        # Mark start and end points
        if len(x) > 0:
            self.ax.scatter([x[0]], [y[0]], [z[0]], c='white', s=50, marker='o', zorder=5)
            self.ax.scatter([x[-1]], [y[-1]], [z[-1]], c=color, s=80, marker='*', zorder=5)
        
        self.orbits[label] = (x, y, z)
        
    def plot_current_position(self, x, y, z, label="Current", color='#ff0000'):
        """Plot a single position marker"""
        self.ax.scatter([x], [y], [z], c=color, s=100, marker='o', label=label)
        
    def finalize(self, title="Orbit View"):
        """Finalize the plot with legend and equal aspect"""
        self.ax.set_title(title, color='#ffffff', fontsize=12, pad=10)
        
        # Try to set equal aspect ratio
        self._set_axes_equal()
        
        # Legend
        self.ax.legend(loc='upper left', facecolor='#2e2e2e', edgecolor='#444444',
                      labelcolor='#ffffff', fontsize=8)
        
        self.canvas.draw()
        
    def _set_axes_equal(self):
        """Set equal aspect ratio for 3D plot"""
        limits = []
        for getter in [self.ax.get_xlim3d, self.ax.get_ylim3d, self.ax.get_zlim3d]:
            limits.append(getter())
        
        ranges = [abs(lim[1] - lim[0]) for lim in limits]
        max_range = max(ranges) / 2.0
        
        centers = [(lim[0] + lim[1]) / 2.0 for lim in limits]
        
        self.ax.set_xlim3d([centers[0] - max_range, centers[0] + max_range])
        self.ax.set_ylim3d([centers[1] - max_range, centers[1] + max_range])
        self.ax.set_zlim3d([centers[2] - max_range, centers[2] + max_range])
    
    def start_animation(self, x, y, z, speed=1, on_frame_callback=None):
        """Start orbit animation"""
        self.stop_animation()
        
        self.anim_x = x
        self.anim_y = y
        self.anim_z = z
        self.anim_frame = 0
        self.anim_speed = max(1, int(speed))
        self.is_animating = True
        self.on_frame_callback = on_frame_callback
        
        # Create marker for current position
        self.anim_marker = self.ax.scatter([x[0]], [y[0]], [z[0]], 
                                            c='#ff0000', s=150, marker='o', 
                                            zorder=10, edgecolors='white', linewidths=2)
        
        # Create trail line
        self.anim_trail, = self.ax.plot([], [], [], color='#ff4444', 
                                         linewidth=3, alpha=0.8)
        
        self.canvas.draw()
        
        def animate(frame):
            if not self.is_animating:
                return self.anim_marker, self.anim_trail
            
            # Update frame index
            idx = (self.anim_frame * self.anim_speed) % len(x)
            self.anim_frame += 1
            
            # Update marker position
            self.anim_marker._offsets3d = ([x[idx]], [y[idx]], [z[idx]])
            
            # Update trail
            trail_start = max(0, idx - self.trail_length)
            self.anim_trail.set_data(x[trail_start:idx+1], y[trail_start:idx+1])
            self.anim_trail.set_3d_properties(z[trail_start:idx+1])
            
            # Callback for frame info
            if self.on_frame_callback:
                self.on_frame_callback(idx, len(x))
            
            return self.anim_marker, self.anim_trail
        
        # Create animation
        self.animation = FuncAnimation(
            self.fig, animate, frames=len(x)//self.anim_speed,
            interval=50, blit=False, repeat=True
        )
        
        self.canvas.draw()
    
    def stop_animation(self):
        """Stop the animation"""
        self.is_animating = False
        
        if self.animation:
            self.animation.event_source.stop()
            self.animation = None
        
        # Clean up animation elements
        if self.anim_marker:
            self.anim_marker.remove()
            self.anim_marker = None
        
        if self.anim_trail:
            self.anim_trail.remove()
            self.anim_trail = None
        
        self.canvas.draw()
    
    def pause_animation(self):
        """Pause/resume animation"""
        if self.animation:
            if self.is_animating:
                self.animation.event_source.stop()
                self.is_animating = False
            else:
                self.animation.event_source.start()
                self.is_animating = True
    
    def set_animation_speed(self, speed):
        """Set animation speed (1-10)"""
        self.anim_speed = max(1, min(10, int(speed)))


class HorizonsUI:
    """Main application UI"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("JPL Horizons Navigator - 3D Orbit Viewer")
        self.root.geometry("1400x900")
        
        self.api = HorizonsAPI()
        self.last_vectors = None  # Store last vector query results
        
        # Configure grid
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        
        # Create main paned window
        self.paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        
        # Left panel: controls + output
        self.left_frame = ttk.Frame(self.paned)
        self.paned.add(self.left_frame, weight=1)
        
        # Right panel: 3D view
        self.right_frame = ttk.LabelFrame(self.paned, text="3D Orbit View", padding=5)
        self.paned.add(self.right_frame, weight=1)
        
        # Setup left panel
        self.left_frame.columnconfigure(0, weight=1)
        self.left_frame.rowconfigure(1, weight=1)
        
        self._create_controls()
        self._create_output()
        
        # Setup 3D plotter
        self.plotter = OrbitPlotter(self.right_frame)
        self._init_plot()
        
        # Set defaults
        self._set_defaults()
    
    def _init_plot(self):
        """Initialize the 3D plot with reference orbits"""
        # Stop any running animation
        if hasattr(self, 'plotter') and self.plotter:
            self.plotter.stop_animation()
            if hasattr(self, 'animate_btn'):
                self.animate_btn.config(text="▶ Animate")
                self.frame_var.set("")
        
        self.plotter.clear()
        self.plotter.plot_sun()
        self.plotter.plot_reference_orbits()
        self.plotter.finalize("Solar System - Heliocentric View")
    
    def _create_controls(self):
        """Create the control panel"""
        control_frame = ttk.LabelFrame(self.left_frame, text="Query Parameters", padding=10)
        control_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        control_frame.columnconfigure(1, weight=1)
        control_frame.columnconfigure(3, weight=1)
        
        row = 0
        
        # === Target Selection ===
        ttk.Label(control_frame, text="Target:", font=("", 10, "bold")).grid(
            row=row, column=0, sticky="w", pady=(0, 5))
        row += 1
        
        ttk.Label(control_frame, text="Preset:").grid(row=row, column=0, sticky="w")
        self.preset_var = tk.StringVar()
        preset_combo = ttk.Combobox(control_frame, textvariable=self.preset_var, width=20)
        preset_combo["values"] = list(HorizonsAPI.MAJOR_BODIES.keys())
        preset_combo.grid(row=row, column=1, sticky="w", padx=5)
        preset_combo.bind("<<ComboboxSelected>>", self._on_preset_select)
        
        ttk.Label(control_frame, text="Custom:").grid(row=row, column=2, sticky="w", padx=(10, 0))
        self.command_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.command_var, width=20).grid(
            row=row, column=3, sticky="w", padx=5)
        
        ttk.Button(control_frame, text="Search", command=self._search_body).grid(
            row=row, column=4, padx=5)
        row += 1
        
        ttk.Separator(control_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=5, sticky="ew", pady=8)
        row += 1
        
        # === Observer Location ===
        ttk.Label(control_frame, text="Center:").grid(row=row, column=0, sticky="w")
        self.center_var = tk.StringVar()
        center_combo = ttk.Combobox(control_frame, textvariable=self.center_var, width=20)
        center_combo["values"] = list(HorizonsAPI.CENTERS.keys())
        center_combo.grid(row=row, column=1, sticky="w", padx=5)
        center_combo.bind("<<ComboboxSelected>>", self._on_center_select)
        
        ttk.Label(control_frame, text="Code:").grid(row=row, column=2, sticky="w", padx=(10, 0))
        self.center_code_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.center_code_var, width=15).grid(
            row=row, column=3, sticky="w", padx=5)
        row += 1
        
        ttk.Separator(control_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=5, sticky="ew", pady=8)
        row += 1
        
        # === Time Settings ===
        ttk.Label(control_frame, text="Start:").grid(row=row, column=0, sticky="w")
        self.start_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.start_var, width=15).grid(
            row=row, column=1, sticky="w", padx=5)
        
        ttk.Label(control_frame, text="Stop:").grid(row=row, column=2, sticky="w", padx=(10, 0))
        self.stop_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.stop_var, width=15).grid(
            row=row, column=3, sticky="w", padx=5)
        row += 1
        
        ttk.Label(control_frame, text="Step:").grid(row=row, column=0, sticky="w")
        self.step_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.step_var, width=10).grid(
            row=row, column=1, sticky="w", padx=5)
        ttk.Label(control_frame, text="(1 d, 6 h, 30 m)").grid(
            row=row, column=2, columnspan=2, sticky="w", padx=5)
        row += 1
        
        ttk.Separator(control_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=5, sticky="ew", pady=8)
        row += 1
        
        # === Ephemeris Type ===
        ttk.Label(control_frame, text="Type:").grid(row=row, column=0, sticky="w")
        self.ephem_type_var = tk.StringVar(value="VECTORS")
        type_frame = ttk.Frame(control_frame)
        type_frame.grid(row=row, column=1, columnspan=3, sticky="w", padx=5)
        
        for etype in ["OBSERVER", "VECTORS", "ELEMENTS"]:
            ttk.Radiobutton(type_frame, text=etype, variable=self.ephem_type_var, 
                           value=etype).pack(side="left", padx=8)
        row += 1
        
        # Quantities (for OBSERVER type)
        ttk.Label(control_frame, text="Quantities:").grid(row=row, column=0, sticky="w")
        self.quantities_entry_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.quantities_entry_var, width=20).grid(
            row=row, column=1, sticky="w", padx=5)
        
        # Reference plane
        ttk.Label(control_frame, text="Ref:").grid(row=row, column=2, sticky="w", padx=(10, 0))
        self.ref_plane_var = tk.StringVar(value="ECLIPTIC")
        ref_combo = ttk.Combobox(control_frame, textvariable=self.ref_plane_var, width=12)
        ref_combo["values"] = ["ECLIPTIC", "FRAME", "BODY EQUATOR"]
        ref_combo.grid(row=row, column=3, sticky="w", padx=5)
        row += 1
        
        # Options
        opt_frame = ttk.Frame(control_frame)
        opt_frame.grid(row=row, column=0, columnspan=5, sticky="w", pady=5)
        
        self.csv_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(opt_frame, text="CSV", variable=self.csv_var).pack(side="left", padx=5)
        
        self.obj_data_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(opt_frame, text="Object data", variable=self.obj_data_var).pack(side="left", padx=5)
        
        self.show_ref_orbits_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt_frame, text="Reference orbits", variable=self.show_ref_orbits_var).pack(side="left", padx=5)
        row += 1
        
        # === Execute ===
        ttk.Separator(control_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=5, sticky="ew", pady=8)
        row += 1
        
        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=row, column=0, columnspan=5, pady=5)
        
        self.query_btn = ttk.Button(btn_frame, text="Get Ephemeris & Plot", 
                                     command=self._execute_query)
        self.query_btn.pack(side="left", padx=5)
        
        ttk.Button(btn_frame, text="Clear", command=self._clear_all).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Save", command=self._save_results).pack(side="left", padx=5)
        row += 1
        
        # === Animation Controls ===
        ttk.Separator(control_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=5, sticky="ew", pady=8)
        row += 1
        
        anim_frame = ttk.LabelFrame(control_frame, text="Animation", padding=5)
        anim_frame.grid(row=row, column=0, columnspan=5, sticky="ew", pady=5)
        
        self.animate_btn = ttk.Button(anim_frame, text="▶ Animate", command=self._toggle_animation)
        self.animate_btn.pack(side="left", padx=5)
        
        self.stop_btn = ttk.Button(anim_frame, text="■ Stop", command=self._stop_animation)
        self.stop_btn.pack(side="left", padx=5)
        
        ttk.Button(anim_frame, text="Reset View", command=self._init_plot).pack(side="left", padx=5)
        
        ttk.Label(anim_frame, text="Speed:").pack(side="left", padx=(15, 5))
        self.speed_var = tk.IntVar(value=1)
        speed_scale = ttk.Scale(anim_frame, from_=1, to=10, variable=self.speed_var,
                                orient="horizontal", length=100, command=self._on_speed_change)
        speed_scale.pack(side="left", padx=5)
        self.speed_label = ttk.Label(anim_frame, text="1x")
        self.speed_label.pack(side="left")
        
        # Frame counter
        self.frame_var = tk.StringVar(value="")
        ttk.Label(anim_frame, textvariable=self.frame_var, foreground="gray").pack(side="left", padx=15)
        row += 1
        
        # Status
        status_frame = ttk.Frame(control_frame)
        status_frame.grid(row=row, column=0, columnspan=5, pady=5)
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(status_frame, textvariable=self.status_var, foreground="blue").pack(side="left")
    
    def _create_output(self):
        """Create the output display"""
        output_frame = ttk.LabelFrame(self.left_frame, text="Results", padding=5)
        output_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)
        
        self.output_text = scrolledtext.ScrolledText(
            output_frame, 
            wrap=tk.NONE,
            font=("Courier", 9),
            bg="#1e1e1e",
            fg="#d4d4d4",
            insertbackground="white",
            height=15
        )
        self.output_text.grid(row=0, column=0, sticky="nsew")
        
        h_scroll = ttk.Scrollbar(output_frame, orient="horizontal", command=self.output_text.xview)
        h_scroll.grid(row=1, column=0, sticky="ew")
        self.output_text.configure(xscrollcommand=h_scroll.set)
    
    def _set_defaults(self):
        """Set default values"""
        now = datetime.utcnow()
        self.start_var.set(now.strftime("%Y-%m-%d"))
        self.stop_var.set((now + timedelta(days=365)).strftime("%Y-%m-%d"))
        self.step_var.set("5 d")
        self.center_var.set("Sun (heliocentric)")
        self.center_code_var.set("500@10")
        self.preset_var.set("Mars")
        self.command_var.set("499")
        self.quantities_entry_var.set("1,9,20,23")
    
    def _on_preset_select(self, event=None):
        name = self.preset_var.get()
        if name in HorizonsAPI.MAJOR_BODIES:
            self.command_var.set(HorizonsAPI.MAJOR_BODIES[name])
    
    def _on_center_select(self, event=None):
        name = self.center_var.get()
        if name in HorizonsAPI.CENTERS:
            self.center_code_var.set(HorizonsAPI.CENTERS[name])
    
    def _search_body(self):
        """Search for a body"""
        cmd = self.command_var.get().strip()
        if not cmd:
            messagebox.showwarning("Input Required", "Enter a search term")
            return
        
        if not cmd.endswith(";"):
            cmd = cmd + ";"
        
        self.status_var.set("Searching...")
        self.root.update()
        
        def do_search():
            try:
                result = self.api.query({
                    "COMMAND": cmd,
                    "OBJ_DATA": "YES",
                    "MAKE_EPHEM": "NO"
                })
                self.root.after(0, lambda: self._display_result(result))
                self.root.after(0, lambda: self.status_var.set("Search complete"))
            except Exception as e:
                import traceback
                error_msg = str(e) if str(e) else traceback.format_exc()
                self.root.after(0, lambda msg=error_msg: self._show_error(msg))
        
        Thread(target=do_search, daemon=True).start()
    
    def _build_params(self) -> dict:
        """Build query parameters from UI state"""
        params = {
            "COMMAND": self.command_var.get().strip(),
            "CENTER": self.center_code_var.get().strip(),
            "START_TIME": self.start_var.get().strip(),
            "STOP_TIME": self.stop_var.get().strip(),
            "STEP_SIZE": self.step_var.get().strip(),
            "EPHEM_TYPE": self.ephem_type_var.get(),
            "OBJ_DATA": "YES" if self.obj_data_var.get() else "NO",
            "MAKE_EPHEM": "YES",
            "CSV_FORMAT": "YES" if self.csv_var.get() else "NO",
        }
        
        etype = self.ephem_type_var.get()
        
        if etype == "OBSERVER":
            quantities = self.quantities_entry_var.get().strip()
            if quantities:
                params["QUANTITIES"] = quantities
        elif etype == "VECTORS":
            params["VEC_TABLE"] = "2"  # State vector
            params["REF_PLANE"] = self.ref_plane_var.get()
            params["OUT_UNITS"] = "AU-D"  # AU and days for plotting
        elif etype == "ELEMENTS":
            params["REF_PLANE"] = self.ref_plane_var.get()
        
        return params
    
    def _execute_query(self):
        """Execute the ephemeris query"""
        params = self._build_params()
        
        if not params["COMMAND"]:
            messagebox.showwarning("Input Required", "Please select or enter a target body")
            return
        
        self.status_var.set("Querying Horizons...")
        self.query_btn.state(["disabled"])
        self.root.update()
        
        def do_query():
            try:
                result = self.api.query(params)
                self.root.after(0, lambda: self._display_result(result))
                self.root.after(0, lambda: self._plot_result(result, params))
                self.root.after(0, lambda: self.status_var.set("Query complete"))
            except Exception as e:
                import traceback
                error_msg = str(e) if str(e) else traceback.format_exc()
                self.root.after(0, lambda msg=error_msg: self._show_error(msg))
            finally:
                self.root.after(0, lambda: self.query_btn.state(["!disabled"]))
        
        Thread(target=do_query, daemon=True).start()
    
    def _display_result(self, result: dict):
        """Display API result"""
        self.output_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.output_text.insert(tk.END, f"ERROR: {result['error']}\n")
            return
        
        if "result" in result:
            self.output_text.insert(tk.END, result["result"])
        else:
            self.output_text.insert(tk.END, json.dumps(result, indent=2))
    
    def _parse_vectors(self, result_text: str):
        """Parse vector ephemeris data from result text"""
        x_vals, y_vals, z_vals = [], [], []
        
        # Find data between $$SOE and $$EOE
        match = re.search(r'\$\$SOE\s*(.*?)\s*\$\$EOE', result_text, re.DOTALL)
        if not match:
            return None
        
        data_block = match.group(1)
        lines = data_block.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('*'):
                continue
            
            # Try to extract X, Y, Z values
            # Format varies, but typically: date X= val Y= val Z= val
            x_match = re.search(r'X\s*=\s*([+-]?\d+\.?\d*E?[+-]?\d*)', line, re.IGNORECASE)
            y_match = re.search(r'Y\s*=\s*([+-]?\d+\.?\d*E?[+-]?\d*)', line, re.IGNORECASE)
            z_match = re.search(r'Z\s*=\s*([+-]?\d+\.?\d*E?[+-]?\d*)', line, re.IGNORECASE)
            
            if x_match and y_match and z_match:
                try:
                    x_vals.append(float(x_match.group(1)))
                    y_vals.append(float(y_match.group(1)))
                    z_vals.append(float(z_match.group(1)))
                except ValueError:
                    continue
        
        if x_vals:
            return np.array(x_vals), np.array(y_vals), np.array(z_vals)
        return None
    
    def _plot_result(self, result: dict, params: dict):
        """Plot the result if it's vector data"""
        if "result" not in result:
            return
        
        etype = params.get("EPHEM_TYPE", "")
        
        if etype != "VECTORS":
            self.status_var.set("Note: Select VECTORS type for 3D plot")
            return
        
        vectors = self._parse_vectors(result["result"])
        if vectors is None:
            self.status_var.set("Could not parse vector data")
            return
        
        x, y, z = vectors
        self.last_vectors = (x, y, z)
        
        # Clear and replot
        self.plotter.clear()
        
        # Check if heliocentric
        center = params.get("CENTER", "")
        is_heliocentric = "10" in center or "0" in center
        
        if is_heliocentric:
            self.plotter.plot_sun()
            if self.show_ref_orbits_var.get():
                # Determine which reference orbits to show based on orbit size
                max_dist = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z)))
                if max_dist < 2:
                    planets = ["Mercury", "Venus", "Earth", "Mars"]
                elif max_dist < 6:
                    planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter"]
                elif max_dist < 12:
                    planets = ["Earth", "Mars", "Jupiter", "Saturn"]
                else:
                    planets = ["Jupiter", "Saturn", "Uranus", "Neptune"]
                self.plotter.plot_reference_orbits(planets)
        
        # Get target name
        target_name = self.preset_var.get() or self.command_var.get()
        color = HorizonsAPI.PLANET_COLORS.get(target_name, '#00ff00')
        
        self.plotter.plot_orbit(x, y, z, label=target_name, color=color)
        
        title = f"{target_name} Orbit"
        if is_heliocentric:
            title += " (Heliocentric)"
        
        self.plotter.finalize(title)
        
        self.status_var.set(f"Plotted {len(x)} points")
    
    def _show_error(self, error: str):
        """Display error message"""
        self.status_var.set("Error")
        messagebox.showerror("API Error", error)
    
    def _clear_all(self):
        """Clear output and reset plot"""
        self._stop_animation()
        self.output_text.delete(1.0, tk.END)
        self._init_plot()
        self.status_var.set("Ready")
    
    def _toggle_animation(self):
        """Start or pause animation"""
        if self.last_vectors is None:
            messagebox.showwarning("No Data", "Get ephemeris data first (VECTORS type)")
            return
        
        if self.plotter.is_animating:
            # Pause
            self.plotter.pause_animation()
            self.animate_btn.config(text="▶ Resume")
            self.status_var.set("Animation paused")
        elif self.plotter.animation is not None:
            # Resume
            self.plotter.pause_animation()
            self.animate_btn.config(text="⏸ Pause")
            self.status_var.set("Animation running")
        else:
            # Start new animation
            x, y, z = self.last_vectors
            speed = self.speed_var.get()
            
            def on_frame(idx, total):
                self.frame_var.set(f"Frame {idx+1}/{total}")
            
            self.plotter.start_animation(x, y, z, speed=speed, on_frame_callback=on_frame)
            self.animate_btn.config(text="⏸ Pause")
            self.status_var.set("Animation running")
    
    def _stop_animation(self):
        """Stop animation completely"""
        self.plotter.stop_animation()
        self.animate_btn.config(text="▶ Animate")
        self.frame_var.set("")
        
        # Replot the static orbit
        if self.last_vectors is not None:
            x, y, z = self.last_vectors
            target_name = self.preset_var.get() or self.command_var.get()
            color = HorizonsAPI.PLANET_COLORS.get(target_name, '#00ff00')
            
            # Just redraw the orbit line
            self.plotter.plot_orbit(x, y, z, label=target_name, color=color)
            self.plotter.canvas.draw()
        
        self.status_var.set("Animation stopped")
    
    def _on_speed_change(self, value):
        """Handle speed slider change"""
        speed = int(float(value))
        self.speed_label.config(text=f"{speed}x")
        self.plotter.set_animation_speed(speed)
    
    def _save_results(self):
        """Save results to file"""
        content = self.output_text.get(1.0, tk.END)
        if not content.strip():
            messagebox.showwarning("No Data", "No results to save")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            with open(filename, "w") as f:
                f.write(content)
            self.status_var.set(f"Saved to {filename}")


def main():
    root = tk.Tk()
    
    style = ttk.Style()
    style.theme_use("clam")
    
    app = HorizonsUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
