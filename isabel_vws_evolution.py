import netCDF4 as nc

ds = nc.Dataset('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/ida_vws_evolution.nc') 

# variable names are valid_time, pressure_level, latitude, longitude(, u, v)

# print(ds)
# print(ds.variables['valid_time'][:])        # check all timesteps
# print(ds.variables['pressure_level'][:])  # should be [850, 200]
# print(ds.variables['u'].shape)        # should be (20, 2, lat, lon) from claude, but it is (44, 2, 181, 201)
###!!!! because there are only 20 timestamps, but we have 4 times per day, and requested 11 days... seems like some timestamps are duplicated, need checking!!!###

# (datetime, lat, lon) -- every 12h from best track
best_track = [  # manually extracted from isabel best tracks
    ('2003-09-06 00:00', 13.80, -31.40, 30),  # (datetime, lat, lon, Vmax)
    ('2003-09-06 12:00', 13.60, -33.90, 40),
    ('2003-09-07 00:00', 13.50, -35.80, 55),
    ('2003-09-07 12:00', 14.40, -37.30, 65),
    ('2003-09-08 00:00', 15.80, -39.70, 80),
    ('2003-09-08 12:00', 17.10, -42.00, 110),
    ('2003-09-09 00:00', 18.20, -44.10, 115),
    ('2003-09-09 12:00', 19.40, -46.30, 115),
    ('2003-09-10 00:00', 20.50, -48.30, 110),
    ('2003-09-10 12:00', 21.10, -50.40, 115),
    ('2003-09-11 00:00', 21.20, -52.30, 125),
    ('2003-09-11 12:00', 21.40, -54.00, 135),
    ('2003-09-12 00:00', 21.60, -55.70, 140),
    ('2003-09-12 12:00', 21.60, -57.40, 140),
    ('2003-09-13 00:00', 21.80, -59.10, 135),
    ('2003-09-13 12:00', 22.10, -61.00, 135),
    ('2003-09-14 00:00', 22.90, -63.30, 135),
    ('2003-09-14 12:00', 23.50, -65.80, 135),
    ('2003-09-15 00:00', 24.30, -67.90, 130),
    ('2003-09-15 12:00', 24.80, -69.40, 120),
    ('2003-09-16 00:00', 25.70, -70.20, 105),
    ('2003-09-16 12:00', 26.80, -70.90, 95),
    ('2003-09-17 00:00', 28.10, -71.50, 95),
    ('2003-09-17 12:00', 29.70, -72.50, 90),
    ('2003-09-18 00:00', 31.50, -73.50, 90),
    ('2003-09-18 12:00', 33.70, -75.20, 90),
    ('2003-09-19 00:00', 36.70, -77.70, 65),
    ('2003-09-19 12:00', 40.90, -80.30, 35),  # from here onwards is extratropical (i.e. transitioned into another weather system)
    ('2003-09-20 00:00', 48.00, -81.00, 25),
] 
best_track = best_track[:28]  
# VWS model, which is the same as the one ins VWS_model.py, but this time we run it over all timesteps
import numpy as np

def compute_vws(ds, time_idx, lat0, lon0_plot):
    lat = ds.variables['latitude'][:]
    lon = ds.variables['longitude'][:]

    # subset region
    lat_mask = (lat >= 5)  & (lat <= 60)
    lon_mask = (lon >= -90) & (lon <= -20)  # extend east to -20 to catch early track
    lat_sub  = lat[lat_mask]
    lon_sub  = lon[lon_mask]

    u850 = ds.variables['u'][time_idx, 0, :, :][np.ix_(lat_mask, lon_mask)]
    v850 = ds.variables['v'][time_idx, 0, :, :][np.ix_(lat_mask, lon_mask)]
    u200 = ds.variables['u'][time_idx, 1, :, :][np.ix_(lat_mask, lon_mask)]
    v200 = ds.variables['v'][time_idx, 1, :, :][np.ix_(lat_mask, lon_mask)]

    lon2d, lat2d = np.meshgrid(lon_sub, lat_sub)
    # lon2d_plot   = lon2d - 360 # shouldn't blindly subtract 360 from everything
    lon2d_plot = np.where(lon2d > 180, lon2d - 360, lon2d)  # correct: convert 0-360 to -180 to 180     

    du      = u200 - u850
    dv      = v200 - v850
    vws_mag = np.sqrt(du**2 + dv**2)

    # annulus mean
    R        = 6371.0
    # lon0_era5 = lon0_plot + 360

    dlat_km = (lat2d - lat0)  * (np.pi/180) * R
    dlon_km = (lon2d - lon0_plot) * (np.pi/180) * R * np.cos(np.radians(lat0))
    # dlat_km  = (lat2d - lat0) * (np.pi/180) * R
    # dlon_km  = (lon2d - lon0_era5) * (np.pi/180) * R * np.cos(np.radians(lat0))
    dist_km  = np.sqrt(dlat_km**2 + dlon_km**2)
    annulus  = (dist_km >= 200) & (dist_km <= 800)

    vws_mean = np.mean(vws_mag[annulus])
    du_mean  = np.mean(du[annulus])
    dv_mean  = np.mean(dv[annulus])
    vws_dir  = (270 - np.degrees(np.arctan2(dv_mean, du_mean))) % 360

    return lon2d_plot, lat2d, vws_mag, du, dv, vws_mean, du_mean, dv_mean, vws_dir

# plot of one frame
import matplotlib.pyplot as plt

def plot_frame(lon2d_plot, lat2d, vws_mag, du, dv,
               vws_mean, du_mean, dv_mean, vws_dir,
               lat0, lon0_plot, title):

    fig, ax = plt.subplots(figsize=(10, 7))

    cf = ax.contourf(lon2d_plot, lat2d, vws_mag,
                     levels=np.arange(0, 42, 2),
                     cmap='jet', extend='max')

    cs = ax.contour(lon2d_plot, lat2d, vws_mag,
                    levels=[5, 10, 15, 20, 25, 30],
                    colors='k', linewidths=0.5, alpha=0.5)
    ax.clabel(cs, fmt='%d m/s', fontsize=7)

    step = 8
    ax.quiver(lon2d_plot[::step, ::step], lat2d[::step, ::step],
              du[::step, ::step], dv[::step, ::step],
              color='k', scale=600, width=0.0015,
              headwidth=4, headlength=4, alpha=0.6)

    for r_km, ls in [(200, '--'), (800, '-')]:
        r_lat = r_km / 111.0
        r_lon = r_km / (111.0 * np.cos(np.radians(lat0)))
        theta = np.linspace(0, 2*np.pi, 300)
        ax.plot(lon0_plot + r_lon*np.cos(theta),
                lat0      + r_lat*np.sin(theta),
                'k', linewidth=1.2, linestyle=ls)

    ax.plot(lon0_plot, lat0, 'k*', markersize=14)
    ax.annotate('',
        xy=(lon0_plot + du_mean*0.7/10, lat0 + dv_mean*0.7/10),
        xytext=(lon0_plot, lat0),
        arrowprops=dict(arrowstyle='->', color='k', lw=2))

    cbar = plt.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label('VWS (m/s)', fontsize=11)

    ax.set_xlim(-90, -20)
    ax.set_ylim(5, 55)
    ax.set_xlabel('Longitude', fontsize=11)
    ax.set_ylabel('Latitude', fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.text(0.02, 0.97,
            f'Annulus VWS: {vws_mean:.1f} m/s @ {vws_dir:.0f}°',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(facecolor='white', edgecolor='grey', boxstyle='round,pad=0.4'))
    ax.grid(True, linestyle='--', alpha=0.4)

    return fig

# loop over all frames 
import os

os.makedirs('frames', exist_ok=True)

ds = nc.Dataset('/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/data/isabel_vws_evolution.nc')

for i, (dt_str, lat0, lon0_plot, vmax) in enumerate(best_track):
    print(f'Processing {dt_str}...')

    lon2d_plot, lat2d, vws_mag, du, dv, \
    vws_mean, du_mean, dv_mean, vws_dir = compute_vws(ds, i, lat0, lon0_plot)

    fig = plot_frame(lon2d_plot, lat2d, vws_mag, du, dv,
                     vws_mean, du_mean, dv_mean, vws_dir,
                     lat0, lon0_plot,
                     title=f'Hurricane Isabel - VWS (200-850 hPa)\n{dt_str} UTC')

    fig.savefig(f'frames/frame_{i:02d}.png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved frame {i}')

    # make gif yay
    import imageio

frames = []
for i in range(len(best_track)):
    frames.append(imageio.imread(f'frames/frame_{i:02d}.png'))

imageio.mimsave('isabel_vws_animation.gif', frames, fps=2)
print('Animation saved.')

# debugging prints
# print(lon2d_plot.min(), lon2d_plot.max())  # should be ~-90 to -40
# print(lat2d.min(), lat2d.max())            # should be ~5 to 50
# print(vws_mag.min(), vws_mag.max())        # should have real values
# print(np.isnan(vws_mag).mean())  # if close to 1.0, this is your problem
# print(lon2d_plot.shape, vws_mag.shape, du.shape)  # all must match