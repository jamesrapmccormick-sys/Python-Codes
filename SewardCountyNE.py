# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 09:45:08 2025

@author: rapha
"""
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiLineString
from pyproj import Transformer
from matplotlib.patches import Circle

# Paths to shapefiles saved in C:/Data/
county_shapefile = r"C:/Data/cb_2024_us_county_5m/cb_2024_us_county_5m.shp"
state_shapefile  = r"C:/Data/cb_2024_us_state_5m/cb_2024_us_state_5m.shp"
roads_shapefile = r"C:/Data/tl_2024_31_prisecroads/tl_2024_31_prisecroads.shp"

# Read shapefiles
counties_us = gpd.read_file(county_shapefile)
states_us   = gpd.read_file(state_shapefile)
roads_us    = gpd.read_file(roads_shapefile)

# Read shapefiles
counties_us = gpd.read_file(county_shapefile)
states_us   = gpd.read_file(state_shapefile)

# Filter Nebraska (FIPS 31)
ne_counties = counties_us[counties_us["STATEFP"] == "31"].copy()
nebraska      = states_us[states_us["STUSPS"] == "NE"].copy()

# Get Seward County polygon
seward = ne_counties[ne_counties["NAME"].isin(["Seward"])].copy()

# Reproject to a good Nebraska projection
aea = "EPSG:5070"
ne_counties = ne_counties.to_crs(aea)
nebraska      = nebraska.to_crs(aea)
seward  = seward.to_crs(aea)
roads_us    = roads_us.to_crs(aea)

# Keep only Interstate, US, and State highways
major_roads = roads_us[roads_us["RTTYP"].isin(["I", "U", "S"])].copy()

# Clip to Cass County (and surrounding area if you like)
roads_seward = gpd.clip(major_roads, seward)

target_substrings = ["80", "34", "6", "103", "15"]
roads_seward = roads_seward[
    roads_seward["FULLNAME"].str.contains("|".join(target_substrings), na=False)
]
# --- Towns ---
towns = [
    {"name": "Seward",   "lat": 40.9079, "lon": -97.0985},
    {"name": "Beaver Crossing",   "lat": 40.7786, "lon": -97.2823},
    {"name": "Bee",          "lat": 41.0064, "lon": -97.0573},
    {"name": "Cordova",           "lat": 40.7172, "lon": -97.3531},
    {"name": "Garland",         "lat": 40.9447, "lon": -96.9856},
    {"name": "Goehner",         "lat": 40.8322, "lon": -97.2212},
    {"name": "Milford", "lat": 40.7742, "lon": -97.0517},
    {"name": "Pleasant Dale", "lat": 40.7919, "lon": -96.9322},
    {"name": "Staplehurst", "lat": 40.9750, "lon": -97.1725},
    {"name": "Utica", "lat": 40.8944, "lon": -97.3461},
]

town_points = gpd.GeoDataFrame(
    towns,
    geometry=[Point(t["lon"], t["lat"]) for t in towns],
    crs="EPSG:4326"
).to_crs(aea)

# Tornado annotations
ann = {
}

# Per-town label offsets
label_offsets = {
    "Seward":  (500,500),
    "Beaver Crossing": (-2000,1000),
    "Utica":  (1000,-000),
    "Staplehurst":  (-2000,-1000),
    "Pleasant Dale": (-4500,1000),
    "Milford": (-1000,-1000),
    "Goehner": (-1000,1000),
    "Bee": (-1000,1000),
    "Cordova": (-1000,1000),
    "Garland": (-1000,1000),
}


default_dx, default_dy = (5000, 5000)

# --- Plot ---
fig, ax = plt.subplots(figsize=(11, 11))

# State fill
nebraska.plot(ax=ax, facecolor="#f6f6f6", edgecolor="black", linewidth=1.2)

# County outlines — dashed & light gray
ne_counties.boundary.plot(ax=ax, color="lightgray", linewidth=0.8, linestyle="--")

# Plot all towns except Seward as red dots
town_points[town_points["name"] != "Seward"].plot(ax=ax, color="black", markersize=30, zorder=5)

# Plot Seward as a blue star
seward = town_points[town_points["name"] == "Seward"]
seward.plot(ax=ax, color="red", markersize=180, marker="*", zorder=6)

# Add labels
for geom, label in zip(town_points.geometry, town_points["name"]):
    x, y = geom.x, geom.y
    dx, dy = label_offsets.get(label, (default_dx, default_dy))


    ax.text(
        x + dx, y + dy,
        label,
        fontsize=9, color="blue", weight="bold", linespacing=1.25
    )

# Interstates, US highways, and state highways in Cass County
interstates = roads_seward[roads_seward["RTTYP"] == "I"]
us_hwy      = roads_seward[roads_seward["RTTYP"] == "U"]
state_hwy   = roads_seward[roads_seward["RTTYP"] == "S"]

interstates.plot(ax=ax, linewidth=2.5, zorder=3)   # can change linestyle if you like
us_hwy.plot(ax=ax, linewidth=2.0, linestyle="-", zorder=3)
state_hwy.plot(ax=ax, linewidth=1.5, linestyle="--", zorder=3)

def get_main_line(geom):
    """Return a single LineString for labeling (longest segment if MultiLineString)."""
    if isinstance(geom, LineString):
        return geom
    if isinstance(geom, MultiLineString):
        return max(list(geom.geoms), key=lambda g: g.length)
    return None

def plot_shield(ax, x, y, text, route_type, dx=0, dy=0):
    """Draw a shield (circle) with text, offset by dx/dy in map units."""
    x2 = x + dx
    y2 = y + dy

    # Different shield size by highway class (optional)
    radius = 0 if route_type == "I" else 0

    shield = Circle(
        (x2, y2),
        radius,
        facecolor="white",
        edgecolor="black",
        linewidth=1.2,
        zorder=6,
    )
    ax.add_patch(shield)

    ax.text(
        x2, y2,
        text,
        fontsize=9,
        weight="bold",
        ha="center",
        va="center",
        zorder=7,
    )

# (pattern in FULLNAME, route type RTTYP, shield text)
label_specs = [
    ("80",  "I", "I 80"),    # Interstate 29
    ("34",  "U", "US 34"),    # US 75
    ("6",  "U", "US 6"),    # US 34
    ("15",  "S", "NE 15"),    # NE 66
    ("103", "S", "NE 103"),   # NE 50
]

for pattern, rttyp, shield_text in label_specs:
    # roads_cass already clipped to Cass County & filtered by substrings
    subset = roads_seward[
        (roads_seward["RTTYP"] == rttyp) &
        (roads_seward["FULLNAME"].str.contains(pattern, na=False))
    ]

    if subset.empty:
        continue

    # merge segments of that route into one geometry
    dissolved = subset.dissolve()

    geom = dissolved.geometry.iloc[0]
    line = get_main_line(geom)
    if line is None:
        continue

    midpt = line.interpolate(0.5, normalized=True)
    plot_shield(
        ax,
        midpt.x,
        midpt.y,
        shield_text,
        rttyp,
        dx=1700,        # east-west offset (+ right, − left)
        dy=1000      # north-south offset (+ up, − down)
    )

# --- Zoom window ---
lon_west, lon_east = -97.4, -96.80
lat_south, lat_north = 40.7, 41.05
transformer = Transformer.from_crs("EPSG:4326", aea, always_xy=True)
mid_lat = (lat_south + lat_north) / 2
mid_lon = (lon_west + lon_east) / 2
x_west, _ = transformer.transform(lon_west, mid_lat)
x_east, _ = transformer.transform(lon_east, mid_lat)
_, y_south = transformer.transform(mid_lon, lat_south)
_, y_north = transformer.transform(mid_lon, lat_north)
ax.set_xlim(x_west, x_east)
ax.set_ylim(y_south, y_north)


ax.set_title("Seward County: The Communities of SCCDP", fontsize=15)
ax.set_axis_off()
plt.tight_layout()
plt.show()
