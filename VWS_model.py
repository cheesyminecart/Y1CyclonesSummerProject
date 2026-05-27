# model for VWS, using tammy (23 Oct, 2023, 12pm, 200hPa and 850hPa data)
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

ds = nc.Dataset('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/tammy_vws.nc')

# Coordinates
lat = ds.variables['latitude'][:]   
lon = ds.variables['longitude'][:]  
levs = ds.variables['pressure_level'][:]  # index 0=850hPa, index 1=200hPa

# [valid_time, pressure_level, latitude, longitude]
u850 = ds.variables['u'][0, 0, :, :]   
v850 = ds.variables['v'][0, 0, :, :]
u200 = ds.variables['u'][0, 1, :, :]
v200 = ds.variables['v'][0, 1, :, :]


# Storm centre at this timestep (Oct 23 12Z from best track)
lat0      = 21.6
lon0_plot = -64.0 # ERA5 uses 0-360, Tammy is at 64°W
lon0_era5 = 296.0  

lat_mask = (lat >= 5)   & (lat <= 45)
lon_mask = (lon >= 270) & (lon <= 360)
lat_sub  = lat[lat_mask]
lon_sub  = lon[lon_mask]

u850_sub = u850[np.ix_(lat_mask, lon_mask)]
v850_sub = v850[np.ix_(lat_mask, lon_mask)]
u200_sub = u200[np.ix_(lat_mask, lon_mask)]
v200_sub = v200[np.ix_(lat_mask, lon_mask)]

# meshgrid on the subsetted coords
lon2d, lat2d = np.meshgrid(lon_sub, lat_sub)
lon2d_plot   = lon2d - 360      

# Distance from storm centre
lon2d, lat2d = np.meshgrid(lon, lat)
R = 6371.0
dlat_km = (lat2d - lat0) * (np.pi/180) * R
dlon_km = (lon2d - lon0_era5) * (np.pi/180) * R * np.cos(np.radians(lat0))
dist_km = np.sqrt(dlat_km**2 + dlon_km**2)


# meshgrid, convert ERA5 lons (0-360) to (-180 to 180)
lon2d, lat2d = np.meshgrid(lon_sub, lat_sub)
lon2d_plot   = lon2d - 360

# VWS
du = u200_sub - u850_sub
dv = v200_sub - v850_sub
vws_mag = np.sqrt(du**2 + dv**2)
vws_dir = (270 - np.degrees(np.arctan2(dv, du))) % 360

# area-averaged VWS in 200-800 km annulus
R       = 6371.0
dlat_km = (lat2d - lat0) * (np.pi/180) * R
dlon_km = (lon2d - lon0_era5) * (np.pi/180) * R * np.cos(np.radians(lat0))
dist_km = np.sqrt(dlat_km**2 + dlon_km**2)
annulus = (dist_km >= 200) & (dist_km <= 800)

vws_mean = np.mean(vws_mag[annulus])
du_mean  = np.mean(du[annulus])
dv_mean  = np.mean(dv[annulus])
vws_dir  = (270 - np.degrees(np.arctan2(dv_mean, du_mean))) % 360

# quiver subsample (every 8th point ~ 2 deg spacing)
step  = 8
lon_q = lon2d_plot[::step, ::step]
lat_q = lat2d[::step, ::step]
du_q  = du[::step, ::step]
dv_q  = dv[::step, ::step]

# --- plot ---
fig, ax = plt.subplots(figsize=(10, 7))

cf = ax.contourf(lon2d_plot, lat2d, vws_mag,
                 levels=np.arange(0, 42, 2),
                 cmap='jet', extend='max')

cs = ax.contour(lon2d_plot, lat2d, vws_mag,
                levels=[5, 10, 15, 20, 25, 30],
                colors='k', linewidths=0.5, alpha=0.5)
ax.clabel(cs, fmt='%d m/s', fontsize=7)

ax.quiver(lon_q, lat_q, du_q, dv_q,
          color='k', scale=600, width=0.0015,
          headwidth=4, headlength=4, alpha=0.6)

# annulus rings
for r_km, ls in [(200, '--'), (800, '-')]:
    r_lat = r_km / 111.0
    r_lon = r_km / (111.0 * np.cos(np.radians(lat0)))
    theta  = np.linspace(0, 2*np.pi, 300)
    ax.plot(lon0_plot + r_lon*np.cos(theta),
            lat0      + r_lat*np.sin(theta),
            'k', linewidth=1.2, linestyle=ls)

# storm centre and mean shear arrow
ax.plot(lon0_plot, lat0, 'k*', markersize=14)
ax.annotate('',
    xy=(lon0_plot + du_mean*0.7/10, lat0 + dv_mean*0.7/10),
    xytext=(lon0_plot, lat0),
    arrowprops=dict(arrowstyle='->', color='k', lw=2))

cbar = plt.colorbar(cf, ax=ax, pad=0.02)
cbar.set_label('VWS (m/s)', fontsize=11)

ax.set_xlim(-90, 0)
ax.set_ylim(5, 45)
ax.set_xlabel('Longitude', fontsize=11)
ax.set_ylabel('Latitude', fontsize=11)
ax.set_title('Hurricane Tammy - VWS (200-850 hPa)\n2023-10-23 12:00 UTC', fontsize=12)
ax.text(0.02, 0.97,
        f'Annulus VWS: {vws_mean:.1f} m/s @ {vws_dir:.0f}°',
        transform=ax.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(facecolor='white', edgecolor='grey', boxstyle='round,pad=0.4'))
ax.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
# plt.savefig('tammy_vws.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"VWS magnitude : {vws_mag:.2f} m/s")
print(f"VWS direction : {vws_dir:.1f}°")