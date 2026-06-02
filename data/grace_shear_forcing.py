# Grace's response timescale to shear forcing
# correlation between shear at T and Vmax at T+dt, where dt is the response timescale
# repeat everything with dt = 0,6,12,18,24,48, r nearest to -1 is best correlated

# obtain grace vws evolution from era5

# first, plot graph of Vmax (y) vs shear (x) for each dt, and then plot correlation coefficient vs dt to find the best dt

import numpy as np
import matplotlib.pyplot as plt

# variable names of evolution .nc fileare valid_time, pressure_level, latitude, longitude(, u, v)

import netCDF4 as nc

# define vws_mag
ds = nc.Dataset('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/grace_vws_evolution.nc') 


u850 = ds.variables['u'][:, 0, :, :]  # shape (time, lat, lon), these are all 3D arrays, like a stack of 2D arrays of u/v over all grids of a certain time
v850 = ds.variables['v'][:, 0, :, :]
u200 = ds.variables['u'][:, 1, :, :]
v200 = ds.variables['v'][:, 1, :, :] 

u_vws = u200 - u850 # shape is (22, 181, 201)
v_vws = v200 - v850
vws_mag = np.sqrt(u_vws**2 + v_vws**2)  # shape is still a 3D array bc each grid has diff value: (22, 181, 201)

# work out annulus-mean VWS (i.e. shear) (scalar at each timestamp)
best_track = [  # manually extracted from Grace (1991) best tracks, every 6h
    ('1991-10-26 00:00', 27.1, -65.2, 30),
    ('1991-10-26 06:00', 27.2, -65.5, 35),
    ('1991-10-26 12:00', 27.3, -66.0, 35),
    ('1991-10-26 18:00', 27.5, -66.5, 40),
    ('1991-10-27 00:00', 28.1, -67.1, 45),
    ('1991-10-27 06:00', 28.9, -66.9, 50),
    ('1991-10-27 12:00', 29.8, -66.4, 55),
    ('1991-10-27 18:00', 30.8, -67.2, 60),
    ('1991-10-28 00:00', 31.6, -68.1, 65),
    ('1991-10-28 06:00', 32.2, -68.5, 65),
    ('1991-10-28 12:00', 32.3, -68.5, 65),
    ('1991-10-28 18:00', 32.4, -67.8, 65),
    ('1991-10-29 00:00', 31.8, -66.8, 70),
    ('1991-10-29 06:00', 31.6, -65.3, 75),
    ('1991-10-29 12:00', 31.5, -63.2, 85),  # last tropical entry
]
# 15 data

R = 6371.0
lat = ds.variables['latitude'][:]
lon = ds.variables['longitude'][:]
lon2d, lat2d = np.meshgrid(lon, lat)

VWS_list = []

best_track_tropical = best_track

for i, (dt_str, lat0, lon0, Vmax_tracks) in enumerate(best_track_tropical):
    # distance mask (distance from storm centre)
    dlat_km = (lat2d - lat0) * (np.pi/180) * R
    dlon_km = (lon2d - lon0) * (np.pi/180) * R * np.cos(np.radians(lat0))
    dist_km = np.sqrt(dlat_km**2 + dlon_km**2)
    annulus = (dist_km >= 200) & (dist_km <= 800)

    vws_t = vws_mag[i, :, :]  # 2D slice at time i
    vws_mean = np.mean(vws_t[annulus])  # single value
    VWS_list.append(vws_mean)

VWS_array = np.array(VWS_list)  # shape = 10 (index 0-9) - one value per timestamp


# define Vmax as the max V of annulus wind speed at 10m above sea level
ds = nc.Dataset('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/grace_surface_V_evolution.nc')  # 28 data sets, contains u10 and v10, 2023-09-06 to 2023-09-19, 00:00 and 12:00

u10 = ds.variables['u10'][:, :, :]  # shape (time, lat, lon)
v10 = ds.variables['v10'][:, :, :]

Vmax_list = []
for i, (dt_str, lat0, lon0, Vmax_tracks) in enumerate(best_track_tropical): 

    # for Vmax, use inner disk instead of annulus
    inner_disk = dist_km <= 200

    wspd10 = np.sqrt(u10[i, :, :]**2 + v10[i, :, :]**2)
    Vmax = np.max(wspd10[inner_disk]) # maximum value within the annulus
    Vmax_list.append(Vmax)

Vmax_array = np.array(Vmax_list)  # shape = 15 (index 0-14) (1D array of Vmax at each timestamp)
# print(Vmax_array[1])
# print(VWS_array)


# correlation between VWS and Vmax for different dt 
# expect r coef to be negative, since higher shear should cause lower Vmax
# before importing detrend, all r values are positive, which may be explained by the fact that both VWS and Vmax are increasing over time. Thus, import detrend to analyze fluctuations around trend
from scipy.signal import detrend

# for the following, choose either detrend or differencing
vws_detrended = detrend(VWS_array)
vmax_detrended = detrend(Vmax_array)
Vmax_tracks_array = np.array([vmax for _, _, _, vmax in best_track_tropical]) * 0.514  # convert knots to m/s

vws_diff  = np.diff(VWS_array)         # change in VWS between timesteps
vmax_diff = np.diff(Vmax_tracks_array) # change in Vmax between timesteps


fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()  # makes it easier to index as axes[0], axes[1], ...

# choose either detrend or differencing

# using detrend
for idx, dt in enumerate([0, 12, 24, 36, 48]):  # 5 values → axes[0] to axes[4]
    steps = dt // 12  # convert hours to timesteps (each step = 12h)

    vws = vws_detrended[:15 - steps] if steps > 0 else vws_detrended     
    # vmax = vmax_detrended[steps:]      # vmax from era5 data
    vmax = Vmax_tracks_array[steps:]     # vmax from best track data [either choose this or era5]

    
    r = np.corrcoef(vws, vmax)[0, 1]
    print(f'dt={dt}h, r={r:.3f}')

    axes[idx].scatter(vws, vmax, label='Data points')
    axes[idx].set_xlabel("Annulus-mean VWS (m/s)")
    axes[idx].set_ylabel("Vmax (m/s)")
    axes[idx].set_title(f"dt = {dt}h, r = {r:.3f}")
    axes[idx].grid()

axes[5].set_visible(False)  # hide the 6th (empty) subplot

plt.suptitle("Vmax vs Annulus-mean VWS for Ida (Detrended Data)", fontsize=14)
plt.tight_layout(pad=3.0)
plt.savefig('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/VWS_shear_forcing_plots/Grace_Vmax_vs_VWS_detrended.png')
plt.show()

# # using differencing
# for idx, dt in enumerate([0, 12, 24, 36, 48]):  # 5 values → axes[0] to axes[4]
#     steps = dt // 12  # convert hours to timesteps (each step = 12h)

#     vws = vws_diff[:14 - steps] if steps > 0 else vws_diff  
#     vmax = vmax_diff[steps:]
    
#     r = np.corrcoef(vws, vmax)[0, 1]
#     print(f'dt={dt}h, r={r:.3f}')

#     axes[idx].scatter(vws, vmax, label='Data points')
#     axes[idx].set_xlabel("Annulus-mean VWS (m/s)")
#     axes[idx].set_ylabel("Vmax (m/s)")
#     axes[idx].set_title(f"dt = {dt}h, r = {r:.3f}")
#     axes[idx].grid()

# axes[5].set_visible(False)  # hide the 6th (empty) subplot

# plt.suptitle("Vmax vs Annulus-mean VWS for Grace (Differenced Data)", fontsize=14)
# plt.tight_layout(pad=3.0)
# plt.savefig('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/VWS_shear_forcing_plots/Grace_Vmax_vs_VWS_differenced.png')
# plt.show()
