
# fitting of tammy starting at line 135
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
# An attempt to model viscous drag force on a single circulation loop

# setting our variables
x = np.linspace(-10, 10, 21)
y = np.linspace(-10, 10, 21)
X, Y = np.meshgrid(x, y)


r = np.sqrt(X**2 + Y**2)

theta = np.arctan2(Y, X)


condition = r == 5

V = np.where(condition,5, 0)

U = -V * np.sin(theta)
W =  V * np.cos(theta)




# .quiver(
#       X, Y,
#       U, W,
#       V,                  # colour by magnitude
#      scale= 3,            # smaller = bigger arrows
#      scale_units='xy',
#      angles='xy'
#     )

# plt.xlim(-10, 10)
# plt.ylim(-10,10)

# plt.xlabel("x")

# plt.ylabel("y")

# plt.title("Single Circulation field")

# plt.colorbar(label="Velocity Magnitude")

# plt.axis("equal")
# plt.show()

# defining initial circulation: 

C_0 = 2 * np.pi * 5 * 6

print (C_0)

# defining dynamic viscosity...

# The formula requires a reference temperature, which we define with the subscript 0 
u_0 = 100
T_0 = 500 
# It also requires a gas-specific constant called sutherland's constant, which is given for air as follows: 
S = 110.4
T = 800

u = u_0 * ((T/T_0)**1.5) * (T_0 + S)/(T + S)

# Now, kinematic viscosity (v) is defined as follows: 
p = 1.3
v = u/p
# where p is the mass density

# As shown on the one note, the laplacian of the velocity vector field needs to computed.  This turns out to be surprisingly easy for a constant vector field in polar coordinates.  The azimuthal compnent of this is as follows:  

V_lap = V * 1/(r**2 + 1e-10)

# As our vector field is circular, the complete line integral, whcih is also the rate of change of circulation, becomes
 
C_d = -v * V_lap * 2 * np.pi * r

# Now, for euler's method, we define a span of time dt, 

dt = 0.05

# testing for the first 3 seconds: 

C_n = C_0 + dt * C_d


idx = np.unravel_index(np.argmin(np.abs(r - 5)), r.shape)
C_n_scalar = C_n[idx]
C_d_scalar = C_d[idx]

print(C_n_scalar)
# defining an array to store our data

V_change= []
time_steps = []

#trying to set up a loop hope this works 
C_n_scalar = C_0 
time = 0
while time < 3:  # Run for 3 seconds
   
    V_n_scalar = C_n_scalar / (2 * np.pi * 5)
    
    V = np.where(np.abs(r - 5) < 0.5, V_n_scalar, 0)

    V_lap = V_n_scalar / (r**2 + 1e-10)  
    C_d = -v * V_lap * 2 * np.pi * r
    C_d_scalar = C_d[idx]
    
    # Step forward
    C_n_scalar = C_n_scalar + dt * C_d_scalar
       
    V_change.append(V_n_scalar)
    time_steps.append(time)
   
    
    time += dt
    print(V_n_scalar)

# plt.figure()
# plt.plot(time_steps, V_change)
# plt.xlabel("Time (s)")
# plt.ylabel("V (m/s)")
# plt.title("V_mag vs Time")
# plt.grid()
# plt.show()

# shear force model ends here

#---------
# trying to fit tammy peak intensity to model
#---------
import os
import glob
import xarray as xr
import pandas as pd

files = sorted(glob.glob("/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/tammy*.nc"))

# sort by actual timestamp
file_times = []
for f in files:
    ds = xr.open_dataset(f)
    t = pd.Timestamp(ds["valid_time"].values[0])
    file_times.append((t, f))

file_times.sort(key=lambda x: x[0])

# extract V_tan at r=5 from each file
times_out = []
V_tammy   = []

target_radius_km = 200  # somewhere noted that radius of max is 5km

for t, f in file_times:
    ds       = xr.open_dataset(f)
    u_data   = ds["u10"].values[0]
    v_data   = ds["v10"].values[0]
    lat      = ds["latitude"].values
    lon      = ds["longitude"].values

    # find eye coor by min wind speed
    speed   = np.sqrt(u_data**2 + v_data**2)
    min_idx = np.unravel_index(np.argmin(speed), speed.shape)
    lat_c   = lat[min_idx[0]]
    lon_c   = lon[min_idx[1]]

    dlat     = np.radians(lat - lat_c) * 6371.0
    dlon     = np.radians((lon - lon_c) * np.cos(np.radians(lat_c))) * 6371.0
    LON_G, LAT_G = np.meshgrid(dlon, dlat)
    R_grid   = np.sqrt(LAT_G**2 + LON_G**2)

    theta_grid = np.arctan2(LAT_G, LON_G)
    V_tan    = -u_data * np.sin(theta_grid) + v_data * np.cos(theta_grid)

    mask     = (R_grid >= target_radius_km - 25) & (R_grid < target_radius_km + 25)
    if mask.sum() > 0:
        times_out.append(t)
        V_tammy.append(np.mean(V_tan[mask]))

# convert time to hours since first timestamp
t0         = times_out[0]
hours      = [(t - t0).total_seconds() / 3600 for t in times_out]

# plot on top of your existing model
plt.figure()
plt.plot(time_steps, V_change, label="Viscous decay model")
plt.plot(hours, V_tammy, 'o--', color='orange', label=f"Tammy V at r={target_radius_km}km")
plt.xlabel("Time (hours from start)")
plt.ylabel("V (m/s)")
plt.title("V_mag vs Time — Model vs Tammy")
plt.legend()
plt.grid()
plt.show()


