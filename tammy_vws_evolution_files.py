# produce nc file that contains wind data at two pressure levels at each timestep
# 28 May: the animation is missing all the contours and quivers, need refining again
import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            'u_component_of_wind',
            'v_component_of_wind',
        ],
        'pressure_level': ['200', '850'],
        'year': '2023',
        'month': '10',
        'day': [
            '18', '19', '20', '21', '22', '23',
            '24', '25', '26', '27', '28'
        ],
        'time': ['00:00', '06:00', '12:00', '18:00'],

        # bounding box around Tammy's full track
        'area': [50, -90, 5, -40],  # N, W, S, E

        'format': 'netcdf',
    },
    'tammy_vws_evolution.nc'
)

# checking
import netCDF4 as nc

ds = nc.Dataset('tammy_vws_evolution.nc')

# variable names are valid_time, pressure_level, latitude, longitude(, u, v)

# print(ds)
# print(ds.variables['valid_time'][:])        # check all timesteps
# print(ds.variables['pressure_level'][:])  # should be [850, 200]
# print(ds.variables['u'].shape)        # should be (20, 2, lat, lon) from claude, but it is (44, 2, 181, 201)
###!!!! because there are only 20 timestamps, but we have 4 times per day, and requested 11 days... seems like some timestamps are duplicated, need checking!!!###

# (datetime, lat, lon) -- every 12h from best track
best_track = [
    ('2023-10-18 00:00', 12.9, -51.0),
    ('2023-10-18 12:00', 13.2, -54.0),
    ('2023-10-19 00:00', 13.0, -52.5),
    ('2023-10-19 12:00', 13.4, -55.3),
    ('2023-10-20 00:00', 13.6, -57.2),
    ('2023-10-20 12:00', 13.9, -58.4),
    ('2023-10-21 00:00', 14.5, -59.6),
    ('2023-10-21 12:00', 15.6, -60.6),
    ('2023-10-22 00:00', 17.5, -61.7),
    ('2023-10-22 12:00', 18.7, -62.8),
    ('2023-10-23 00:00', 20.2, -63.9),
    ('2023-10-23 12:00', 21.6, -64.0),
    ('2023-10-24 00:00', 22.6, -63.5),
    ('2023-10-24 12:00', 23.6, -62.6),
    ('2023-10-25 00:00', 24.5, -61.2),
    ('2023-10-25 12:00', 26.1, -59.4),
    ('2023-10-26 00:00', 29.1, -57.8),
    ('2023-10-26 12:00', 30.5, -58.8),
    ('2023-10-27 00:00', 31.5, -60.1),
    ('2023-10-27 12:00', 32.0, -61.0),
]

# VWS model, which is the same as the one ins VWS_model.py, but this time we run it over all timesteps
import numpy as np

def compute_vws(ds, time_idx, lat0, lon0_plot):
    lat = ds.variables['latitude'][:]
    lon = ds.variables['longitude'][:]

    # subset region
    lat_mask = (lat >= 5)   & (lat <= 55)
    # lon_mask = (lon >= 270) & (lon <= 360)
    lon_mask = (lon >= -90) & (lon <= -40)
    lat_sub  = lat[lat_mask]
    lon_sub  = lon[lon_mask]

    u850 = ds.variables['u'][time_idx, 0, :, :][np.ix_(lat_mask, lon_mask)]
    v850 = ds.variables['v'][time_idx, 0, :, :][np.ix_(lat_mask, lon_mask)]
    u200 = ds.variables['u'][time_idx, 1, :, :][np.ix_(lat_mask, lon_mask)]
    v200 = ds.variables['v'][time_idx, 1, :, :][np.ix_(lat_mask, lon_mask)]

    lon2d, lat2d = np.meshgrid(lon_sub, lat_sub)
    lon2d_plot   = lon2d - 360

    du      = u200 - u850
    dv      = v200 - v850
    vws_mag = np.sqrt(du**2 + dv**2)

    # annulus mean
    R        = 6371.0
    lon0_era5 = lon0_plot + 360
    dlat_km  = (lat2d - lat0) * (np.pi/180) * R
    dlon_km  = (lon2d - lon0_era5) * (np.pi/180) * R * np.cos(np.radians(lat0))
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

    ax.set_xlim(-90, -40)
    ax.set_ylim(5, 50)
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

ds = nc.Dataset('tammy_vws_evolution.nc')

for i, (dt_str, lat0, lon0_plot) in enumerate(best_track):
    print(f'Processing {dt_str}...')

    lon2d_plot, lat2d, vws_mag, du, dv, \
    vws_mean, du_mean, dv_mean, vws_dir = compute_vws(ds, i, lat0, lon0_plot)

    fig = plot_frame(lon2d_plot, lat2d, vws_mag, du, dv,
                     vws_mean, du_mean, dv_mean, vws_dir,
                     lat0, lon0_plot,
                     title=f'Hurricane Tammy - VWS (200-850 hPa)\n{dt_str} UTC')

    fig.savefig(f'frames/frame_{i:02d}.png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved frame {i}')

    # make gif yay
    import imageio

frames = []
for i in range(len(best_track)):
    frames.append(imageio.imread(f'frames/frame_{i:02d}.png'))

imageio.mimsave('tammy_vws_animation.gif', frames, fps=2)
print('Animation saved.')