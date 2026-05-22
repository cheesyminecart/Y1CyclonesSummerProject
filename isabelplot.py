import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


# Load ERA5 data

ds = xr.open_dataset("C:/Users/Jack/OneDrive - Imperial College London/isabel3.nc")

# Extract variables
mslp = ds['msl'][0] / 100        
u10  = ds['u10'][0]
v10  = ds['v10'][0]

lats = ds['latitude']
lons = ds['longitude']


#plot the wind vectors and pressure contours


plt.figure(figsize=(12,10))

# Pressure contours
contours = plt.contour(lons, lats, mslp, levels=30, cmap='coolwarm')
plt.clabel(contours, inline=True, fontsize=8)

# Wind vectors
plt.quiver(lons, lats, u10, v10, scale=700)


# Overlay the strom centre from tcpos.py (thanks Julien)


storm_lat = 35.1
storm_lon = -76.4

plt.plot(storm_lon, storm_lat, 'ko', markersize=10, label="Isabel centre")


#labels and title


plt.title("Hurricane Isabel – ERA5 – 18 Sep 2003 18:00 UTC", fontsize=14)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.legend()

plt.show()
