import cdsapi
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os

# output folders (ZOOMED OUT VERSION)
era5_folder = r"C:\Users\Jack\OneDrive - Imperial College London\ERA5_Tammy_zoomout"
png_folder  = r"C:\Users\Jack\OneDrive - Imperial College London\Tammy_frames_zoomout"

os.makedirs(era5_folder, exist_ok=True)
os.makedirs(png_folder, exist_ok=True)

# Tammy timestamps 
timestamps = [
    "2023-10-18 18:00","2023-10-19 00:00","2023-10-19 06:00","2023-10-19 12:00",
    "2023-10-19 18:00","2023-10-20 00:00","2023-10-20 06:00","2023-10-20 12:00",
    "2023-10-20 18:00","2023-10-21 00:00","2023-10-21 06:00","2023-10-21 12:00",
    "2023-10-21 18:00","2023-10-22 00:00","2023-10-22 06:00","2023-10-22 12:00",
    "2023-10-22 18:00","2023-10-23 00:00","2023-10-23 06:00","2023-10-23 12:00",
    "2023-10-23 18:00","2023-10-24 00:00","2023-10-24 06:00","2023-10-24 12:00",
    "2023-10-24 18:00","2023-10-25 00:00","2023-10-25 06:00","2023-10-25 12:00",
    "2023-10-25 18:00","2023-10-26 00:00","2023-10-26 06:00","2023-10-26 12:00",
    "2023-10-26 18:00","2023-10-27 00:00","2023-10-27 06:00","2023-10-27 12:00",
    "2023-10-27 18:00","2023-10-28 00:00","2023-10-28 06:00","2023-10-28 12:00",
    "2023-10-28 18:00","2023-10-29 00:00","2023-10-29 06:00","2023-10-29 12:00",
    "2023-10-29 18:00","2023-10-30 00:00","2023-10-30 06:00","2023-10-30 12:00",
    "2023-10-30 18:00","2023-10-31 00:00","2023-10-31 06:00","2023-10-31 12:00",
    "2023-10-31 18:00"
]

storm_lats = [
    12.9,13.0,13.2,13.4,13.5,13.6,13.7,13.9,14.1,14.5,14.9,15.6,16.6,17.5,
    18.0,18.7,19.4,20.2,21.0,21.6,22.1,22.6,23.1,23.6,24.0,24.5,25.1,26.1,
    27.5,29.1,30.0,30.5,31.0,31.5,31.8,32.0,32.3,32.6,32.9,33.1,33.3,32.8,
    32.6,31.8,30.9,29.6,28.5,27.5,26.5,25.6,25.1,24.9,24.6
]

storm_lons = [
    -51.0,-52.5,-54.0,-55.3,-56.4,-57.2,-57.9,-58.4,-58.9,-59.6,-60.3,-60.6,
    -61.0,-61.7,-62.2,-62.8,-63.4,-63.9,-64.0,-64.0,-63.8,-63.5,-63.1,-62.6,
    -62.0,-61.2,-60.4,-59.4,-58.3,-57.8,-58.3,-58.8,-59.5,-60.1,-60.6,-61.0,
    -61.4,-61.5,-61.0,-59.5,-58.0,-56.1,-54.3,-52.5,-50.9,-49.6,-49.0,-48.5,
    -48.1,-48.5,-49.0,-49.9,-50.6
]

# ZOOMED‑OUT domain (wider like Isabel)
north, south, west, east = 50, 5, -85, -35

# ERA5 API client
api = cdsapi.Client()

# loop through all Tammy frames
for i in range(len(timestamps)):

    timestamp = timestamps[i]
    lat = storm_lats[i]
    lon = storm_lons[i]

    # split timestamp into date + hour
    date_str, time_str = timestamp.split(" ")
    year, month, day = date_str.split("-")
    hour = time_str

    # output filename for this frame
    nc_file = os.path.join(era5_folder, f"tammy_{i+1:02d}.nc")

    # remove old file if it exists
    if os.path.exists(nc_file):
        os.remove(nc_file)

    # download ERA5 fields for this timestamp
    api.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': [
                'mean_sea_level_pressure',
                '10m_u_component_of_wind',
                '10m_v_component_of_wind'
            ],
            'year': year,
            'month': month,
            'day': day,
            'time': hour,
            'area': [north, west, south, east],
            'grid': [0.25, 0.25],
            'data_format': 'netcdf'
        },
        nc_file
    )

    # open ERA5 file
    ds = xr.open_dataset(nc_file)

    # convert longitudes to -180 → 180
    ds = ds.assign_coords(
        longitude=((ds.longitude + 180) % 360) - 180
    ).sortby("longitude")

    # extract fields
    mslp = ds["msl"][0] / 100.0     
    u10  = ds["u10"][0]
    v10  = ds["v10"][0]

    lats = ds["latitude"]
    lons = ds["longitude"]

    # create the figure
    plt.figure(figsize=(12, 10))

    # MSLP coloured contours 
    cs = plt.contour(lons, lats, mslp, levels=30, cmap="coolwarm")
    plt.clabel(cs, inline=True, fontsize=8)

    # Wind vectors (quiver plot)
    plt.quiver(lons, lats, u10, v10, scale=700)

    # Storm centre marker 
    plt.plot(lon, lat, "ko", markersize=10, label="Tammy centre")

    plt.legend(loc="upper right")

    plt.title(f"Hurricane Tammy – ERA5 – {timestamp}")
    plt.xlabel("Longitude (°)")
    plt.ylabel("Latitude (°)")

    plt.savefig(
        os.path.join(png_folder, f"frame_{i+1:02d}.png"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
