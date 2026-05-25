import xarray as xr

ds = xr.open_dataset("/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/isabel1.nc")

u = ds["u10"]
v = ds["v10"]

print(u)
print(v)

u_data = ds["u10"].values
v_data = ds["v10"].values

print(u_data.shape)
print(v_data.shape)

u2d = u_data[0]
v2d = v_data[0]

lat = ds["latitude"].values
lon = ds["longitude"].values

import numpy as np

# compute magnitude of wind speed of each grid point
speed = np.sqrt(u2d**2 + v2d**2)
print(speed)

# 2.locate cyclone center
# best approach is to find minimum pressure, but we don't have from data
# our practical approach: min wind speed, estimate center as the point around which winds rotate
import numpy as np

idx = np.unravel_index(np.argmin(speed), speed.shape)

print(idx)

center_lat = lat[idx[0]]
center_lon = lon[idx[1]]
print(f"Estimated cyclone center: Latitude {center_lat}, Longitude {center_lon}")

# 3. compute radius of each grid point from estimated center
import numpy as np

# create 2D coordinate grids
lon2d, lat2d = np.meshgrid(lon, lat)

lat_rad = np.radians(lat2d) # convert to radians
lon_rad = np.radians(lon2d)

clat = np.radians(center_lat)
clon = np.radians(center_lon)

# Earth radius (km)
R = 6371

# differences
dlat = lat_rad - clat
dlon = lon_rad - clon

# approximate distance
r = R * np.sqrt(dlat**2 + (np.cos(clat)*dlon)**2)
# print(r)


# 4. fit with holland model
from scipy.optimize import curve_fit
import numpy as np


r_flat = r.flatten()
v_flat = speed.flatten()

# remove center singularity
mask = r_flat > 1

r_fit = r_flat[mask]
v_fit = v_flat[mask]

def holland(r, vmax, rm, B):
    return vmax * (rm/r)**(B/2) * np.exp(
        0.5 * (1 - (rm/r)**B)
    )


popt, pcov = curve_fit(
    holland,
    r_fit,
    v_fit,
    p0=[50, 30, 1.5]
)

print(popt)
print('max wind speed (vmax):', popt[0], 'm/s')
print('radius of maximum winds (rm):', popt[1], 'km')
print('Holland Parameter (B):', popt[2])

# 6. evaluate error
v_model = holland(r_fit, *popt)

# we assumed that cyclone is symmetric 
rmse = np.sqrt(np.mean((v_fit - v_model)**2))
print('RMSE=', rmse)