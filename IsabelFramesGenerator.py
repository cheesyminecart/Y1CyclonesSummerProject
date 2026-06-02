import cdsapi
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os


#Output folders

era5_folder = r"C:\Users\Jack\OneDrive - Imperial College London\ERA5 data isabel"
png_folder  = r"C:\Users\Jack\OneDrive - Imperial College London\animation frames isabel"

os.makedirs(era5_folder, exist_ok=True)
os.makedirs(png_folder, exist_ok=True)

#Track and timestamps

timestamps = [
    "2003-09-06 06:00",
    "2003-09-06 12:00",
    "2003-09-06 18:00",
    "2003-09-07 00:00",
    "2003-09-07 06:00",
    "2003-09-07 12:00",
    "2003-09-07 18:00",
    "2003-09-08 00:00",
    "2003-09-08 06:00",
    "2003-09-08 12:00",
    "2003-09-08 18:00",
    "2003-09-09 00:00",
    "2003-09-09 06:00",
    "2003-09-09 12:00",
    "2003-09-09 18:00",
    "2003-09-10 00:00",
    "2003-09-10 06:00",
    "2003-09-10 12:00",
    "2003-09-10 18:00",
    "2003-09-11 00:00",
    "2003-09-11 06:00",
    "2003-09-11 12:00",
    "2003-09-11 18:00",
    "2003-09-12 00:00",
    "2003-09-12 06:00",
    "2003-09-12 12:00",
    "2003-09-12 18:00",
    "2003-09-13 00:00",
    "2003-09-13 06:00",
    "2003-09-13 12:00",
    "2003-09-13 18:00",
    "2003-09-14 00:00",
    "2003-09-14 06:00",
    "2003-09-14 12:00",
    "2003-09-14 18:00",
    "2003-09-15 00:00",
    "2003-09-15 06:00",
    "2003-09-15 12:00",
    "2003-09-15 18:00",
    "2003-09-16 00:00",
    "2003-09-16 06:00",
    "2003-09-16 12:00",
    "2003-09-16 18:00",
    "2003-09-17 00:00",
    "2003-09-17 06:00",
    "2003-09-17 12:00",
    "2003-09-17 18:00",
    "2003-09-18 00:00",
    "2003-09-18 06:00",
    "2003-09-18 12:00",
    "2003-09-18 17:00",
    "2003-09-18 18:00",
    "2003-09-19 00:00",
    "2003-09-19 06:00",
]

storm_lats = [
    13.9, 13.6, 13.4, 13.5, 13.9, 14.4, 15.2, 15.8, 16.5,
    17.1, 17.6, 18.2, 18.9, 19.4, 20.0, 20.5, 20.9, 21.1,
    21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.6, 21.7,
    21.8, 21.9, 22.1, 22.5, 22.9, 23.2, 23.5, 23.9, 24.3,
    24.5, 24.8, 25.3, 25.7, 26.3, 26.8, 27.4, 28.1, 28.9,
    29.7, 30.6, 31.5, 32.5, 33.7, 34.9, 35.1, 36.7, 38.6
]

storm_lons = [
    -32.7, -33.9, -34.9, -35.8, -36.5, -37.3, -38.5, -39.7, -40.9,
    -42.0, -43.1, -44.1, -45.2, -46.3, -47.3, -48.3, -49.4, -50.4,
    -51.4, -52.3, -53.2, -54.0, -54.8, -55.7, -56.6, -57.4, -58.2,
    -59.1, -60.1, -61.0, -62.1, -63.3, -64.6, -65.8, -67.0, -67.9,
    -68.8, -69.4, -69.8, -70.2, -70.5, -70.9, -71.2, -71.5, -71.9,
    -72.5, -73.0, -73.5, -74.3, -75.2, -76.2, -76.4, -77.7, -78.9
]

n_frames = len(timestamps)

#Domain 

north = 50
south = 10
west  = -90
east  = -20

#ERA5 client

api = cdsapi.Client()


#Loop over frames

for i in range(n_frames):

    timestamp = timestamps[i]
    storm_lat = storm_lats[i]
    storm_lon = storm_lons[i]

    date_str, time_str = timestamp.split(" ")
    year, month, day = date_str.split("-")
    hour = time_str  # already "HH:MM"

    print(f"\n=== Frame {i+1}/{n_frames} — {timestamp} ===")

    nc_file = os.path.join(era5_folder, f"isabel_{i+1:02d}.nc")

#Download ERA5 data
    if os.path.exists(nc_file):
        os.remove(nc_file)

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
            'format': 'netcdf'
        },
        nc_file
    )

    print("Downloaded ERA5.")

    #Plot frame 
    ds = xr.open_dataset(nc_file)

    #Convert longitudes from 0–360 to -180–180
    ds = ds.assign_coords(
        longitude=((ds.longitude + 180) % 360) - 180
    ).sortby("longitude")

    mslp = ds["msl"][0] / 100.0
    u10  = ds["u10"][0]
    v10  = ds["v10"][0]

    lats = ds["latitude"]
    lons = ds["longitude"]

    plt.figure(figsize=(12, 10))

    contours = plt.contour(lons, lats, mslp, levels=30, cmap="coolwarm")
    plt.clabel(contours, inline=True, fontsize=8)

    plt.quiver(lons, lats, u10, v10, scale=700)

    plt.plot(storm_lon, storm_lat, "ko", markersize=10)

    plt.title(f"Hurricane Isabel – ERA5 – {timestamp}", fontsize=14)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    save_path = os.path.join(png_folder, f"frame_{i+1:02d}.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved {save_path}")
