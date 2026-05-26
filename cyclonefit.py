import numpy as np
import scipy as sc
import xarray as xr
import matplotlib.pyplot as plt
import os
import cartopy  # Map projections libary
import cartopy.crs as ccrs  # Projections list 
from scipy.optimize import curve_fit
from itertools import product

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
ds = xr.open_dataset(r'C:\\Users\\ypara\\OneDrive\\Desktop\\Cyclone_fitting\\tammy_peakintensity.nc')
#print(ds)

u = ds.u10.sel(valid_time="2003-09-07T12:00:00", method='nearest').load()
v = ds.v10.sel(valid_time="2003-09-07T12:00:00", method='nearest').load()


# let us define a wind speed

W = np.sqrt(u**2 + v**2)

# W here is a 2D array

# The centre of the cyclone is defined as the point with the least pressure, however, the file we have at the moment does not have pressure data.  Therefore, we will be finding the minimum velocity in a given area to define our centre 

# as the maximum is usually at the eyewall, we start by making a 3 degree box around the maximum velocity 

# defining the boundary values 
boundary_idx = np.unravel_index(W.values.argmax(), W.values.shape)
boundary_lat  = float(W.latitude[boundary_idx[0]])
boundary_lon  = float(W.longitude[boundary_idx[1]])

# creating a 2d 'slice' about 3 degrees larger on each side 
box = W.sel(
    latitude=slice(boundary_lat + 3, boundary_lat - 3),   
    longitude=slice(boundary_lon - 3, boundary_lon + 3)
)

# finding the centre by findint he minimum velocity inside the box 
centre_idx_box = np.unravel_index(box.values.argmin(), box.values.shape)
centre_lat = float(box.latitude[centre_idx_box[0]])
centre_lon = float(box.longitude[centre_idx_box[1]])

#defining our cross section to be the central longitude 
cross_section = W.sel(longitude=centre_lon, method='nearest')

# defining distance from the centre, in terms of latitude 

lats = cross_section.latitude.values
dist_km = (lats - centre_lat) * 111.32 

W_profile = cross_section.values

# defining a threshold to end the array, and ordering the distance array 
threshold = 2.0  # m/s — adjust if needed
above_threshold = np.abs(W_profile) > threshold
first = np.argmax(above_threshold)
last  = len(above_threshold) - np.argmax(above_threshold[::-1]) - 1
dist_trimmed    = dist_km[first:last+1]
profile_trimmed = W_profile[first:last+1]

# complete plot (both sides)
# fig, ax = plt.subplots(figsize=(10, 5))

# ax.plot(dist_trimmed, profile_trimmed, color='steelblue', linewidth=2)

# ax.set_xlabel('Distance from centre (km)', fontsize=12)
# ax.set_ylabel('Wind speed (m/s)', fontsize=12)
# ax.set_title(f'Radial wind speed profile at peak intensity',)

# plt.show()

# we only need half of the graph for fitting, so we plot the left half-- 

# Filtering positive latitudes 
pos_mask = dist_trimmed > 0
dist_pos    = dist_trimmed[pos_mask]
profile_pos = profile_trimmed[pos_mask]

# Drawing our plot 
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(dist_pos, profile_pos, color='steelblue', linewidth=2)

ax.set_xlabel('Distance from centre (km)', fontsize=12)
ax.set_ylabel('Wind speed (m/s)', fontsize=12)
ax.set_title('Radial wind speed profile at peak intensity (northern side)', fontsize=13)

ax.set_xlim(left=0)

plt.show()

# fitting to the wood white model


def Wood(r, V_X, R_X, k, n, l):
    eps = 1e-7
    V = (V_X
         * (n ** l)
         * ((r + eps) / R_X) ** k
         / (n - k + k * ((r + eps) / R_X) ** (n / l)) ** l)
    return V

peak_idx = np.argmax(profile_pos)

# define parameter ranges for brute force method (as these factors are qualitative this is sort of the only thing that works )
k_range = np.linspace(1, 20, 50)
n_range = np.linspace(1.0, 15, 50)
l_range = np.linspace(0.5, 20, 50)

# Fix V_X and R_X from the data directly
V_X_fixed = profile_pos.max()
R_X_fixed  = dist_pos[peak_idx]

# initialising the root mean square error (rmse)
best_rmse   = np.inf
best_params = None

# setting up a loop for brute force!!!!

for k, n, l in product(k_range, n_range, l_range): # the number of computations it at LEAST the product of the numbers of k, l and n that are being tested
    if k >= n: 
        continue # if k is greater than or equal to n, the loop skips that combination immediately, as k must be less than n 
    try:
        W_pred = Wood(dist_pos, V_X_fixed, R_X_fixed, k, n, l) # W_pred is the array of predictions based on random values of k, n and l that satisfy the above conditions 
        if np.any(np.isnan(W_pred)):   # if any value in the array is NaN, then it skips that combination of k, n and l 
            continue
        rmse = np.sqrt(np.mean((W_pred - profile_pos) ** 2)) # this part computes the rmse and checks whether its less than the initial/previous rmse.  This is why the initial rmse is deinfed to be infinity. as any value here must be lower than that.
        if rmse < best_rmse:
            best_rmse   = rmse # updates the value of best_rmse 
            best_params = (V_X_fixed, R_X_fixed, k, n, l) # defines a new set of parameters that give said best rmse 
    except Exception:  # if the Wood function crashes/overflows, it skips that combination. 
        continue

# printing our initial guesses 
print(f"Best k    = {best_params[2]}")
print(f"Best n    = {best_params[3]}")
print(f"Best l    = {best_params[4]}")
print(f"Best RMSE = {best_rmse}")


# defining p0 for curve fit 
p0 = list(best_params)
params, cov = curve_fit(
    Wood, dist_pos, profile_pos,
    p0=p0,
    bounds=(
        [0,   0,   0,   0,   0],   # lower: all positive
        [200, 500, 20,  15,  20]   # upper bounds based on the bounds of the brute force calculations
    ),
    maxfev=10000
)


# printing new params from our curve fit 
V_X_fit, R_X_fit, k_fit, n_fit, l_fit = params
print(f"\nRefined V_X = {V_X_fit} m/s")
print(f"Refined R_X = {R_X_fit} km")
print(f"Refined k   = {k_fit}")
print(f"Refined n   = {n_fit}")
print(f"Refined l   = {l_fit}")


# defining a new set of r values for plotting according to the new parameters
r_smooth = np.linspace(dist_pos.min(), dist_pos.max(), 500)
W_fit    = Wood(r_smooth, *params)


# graphing the fit as well as the true velocity profile 
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(dist_pos, profile_pos, color='steelblue', linewidth=2, label='ERA5 data')
ax.plot(r_smooth, W_fit, color='firebrick', linewidth=2, linestyle='--', label='Wood-White fit')
ax.set_xlabel('Distance from centre (km)', fontsize=12)
ax.set_ylabel('Wind speed (m/s)', fontsize=12)
ax.set_title('Wood-White model fit — northern side', fontsize=13)
ax.legend()
plt.tight_layout()
plt.show()

# quiver plot!!
x = np.linspace(-500, 500, 40)   # km
y = np.linspace(-500, 500, 40)   # km
X, Y = np.meshgrid(x, y)
r_grid = np.sqrt(X**2 + Y**2)
theta  = np.arctan2(Y, X)


V_X_fit, R_X_fit, k_fit, n_fit, l_fit = params


eps  = 1e-7
p    = (r_grid + eps) / R_X_fit
phi  = (n_fit ** l_fit * p ** k_fit) / (n_fit - k_fit + k_fit * p ** (n_fit / l_fit)) ** l_fit
V_2D = V_X_fit * phi

U = -V_2D * np.sin(theta)
W =  V_2D * np.cos(theta)


fig, ax = plt.subplots(figsize=(8, 8))

q = ax.quiver(
    X, Y, U, W,
    V_2D,
    scale=0.5,
    scale_units='xy',
    angles='xy',
    cmap='viridis'
)

plt.colorbar(q, ax=ax, label='Wind speed (m/s)')

ax.set_xlabel('Distance from centre (km)', fontsize=12)
ax.set_ylabel('Distance from centre (km)', fontsize=12)
ax.set_title('Wood-White fitted cyclone wind field', fontsize=13)
ax.set_aspect('equal')
ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')
ax.plot(0, 0, 'r+', markersize=12, markeredgewidth=2, label='Centre')
ax.legend()

plt.tight_layout()
plt.show()


