import numpy as np
import scipy as sc
import xarray as xr
import matplotlib.pyplot as plt
import os


# Set working directory to the folder the script is in
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# Reference latitude for the cos correction (in degrees)
ref_lat = 26.1

# Convert degrees to kilometres
latitude_km  = 27.83 * 111.32
longitude_km = 27.83 * 111.32 * np.cos(np.radians(ref_lat))

print(f"Latitude (km):  {latitude_km:.2f}")
print(f"Longitude (km): {longitude_km:.2f}")

# Open and inspect the dataset
ds = xr.open_dataset('C:\\Users\\ypara\\OneDrive\\Desktop\\Cyclone_fitting\\tammy_peakintensity.nc')

u = ds.u10.sel(valid_time="2003-09-07T12:00:00", method='nearest').load()
v = ds.v10.sel(valid_time="2003-09-07T12:00:00", method='nearest').load()

W = np.sqrt(u**2 + v**2)


W = np.sqrt(u**2 + v**2)

# Make a 3 degree box around the rough maximum first, then find the minimum inside, this makes up a rough approximation of the cyclone centre. 
rough_idx = np.unravel_index(W.values.argmax(), W.values.shape)
rough_lat  = float(W.latitude[rough_idx[0]])
rough_lon  = float(W.longitude[rough_idx[1]])

box = W.sel(
    latitude=slice(rough_lat + 3, rough_lat - 3),  
    longitude=slice(rough_lon - 3, rough_lon + 3)
)

# defining centre latitude and longitude 
centre_idx_box = np.unravel_index(box.values.argmin(), box.values.shape)
centre_lat = float(box.latitude[centre_idx_box[0]])
centre_lon = float(box.longitude[centre_idx_box[1]])


# defining the cross sectional velocity profile 
cross_section = W.sel(longitude=centre_lon, method='nearest')
W_profile = cross_section.values

lats = cross_section.latitude.values
dist_km = (lats - centre_lat) * 111.32


# defining a threshold for minimum velocity 
threshold = 1e-2
above_threshold = np.abs(W_profile) > threshold
indices = np.where(above_threshold)[0]
first = indices[0]
last  = indices[-1]

# making a new profile accounting for the threshold 
dist_trimmed    = dist_km[first:last+1]
W_profile_trimmed = W_profile[first:last+1]

# Taking one side (here the negative, as it had more data) and drawing out the cross sectional profile 

dist_negative    = dist_trimmed[dist_trimmed <= 0]
W_profile_negative = W_profile_trimmed[dist_trimmed <= 0]

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(-dist_negative, W_profile_negative, color='steelblue', linewidth=2)
ax.axhline(threshold, color='gray', linestyle=':', linewidth=1, label=f'Threshold ({threshold} m/s)')
ax.set_xlabel('Distance from centre (km)', fontsize=12)
ax.set_ylabel('Wind speed (m/s)', fontsize=12)
ax.set_title(f'Radial wind speed profile', fontsize=13)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

