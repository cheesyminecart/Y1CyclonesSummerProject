import numpy as np
import matplotlib.pyplot as plt

N = 200
x = np.linspace(-100, 100, N)  # km
y = np.linspace(-100, 100, N)

X, Y = np.meshgrid(x, y)

# changing alpha ratios
V_max = 50  # max wind speed (i.e. peak intensity?)
R_max = 20  # radius of max wind speed 
B = 1.2  # Holland B-parameter: increases with storm intensity

def wind_field(X, Y):  # right now currently only depend on r
    r = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)

    V = V_max * np.sqrt((R_max/r)**B * np.exp(1 - (R_max/r)**B))

    u = -V * np.sin(theta)
    v = V * np.cos(theta)
    return (u, v)
#V[r==0] = 0

(u, v) = wind_field(X, Y)
speed = np.sqrt(u**2 + v**2)

# plt.figure(figsize=(8,8))

# plt.contourf(X, Y, speed, levels=50)  # Wind speed background

# plt.streamplot(X, Y, u, v, color='white', density=2)  # Streamlines

# plt.xlabel("x (km)")
# plt.ylabel("y (km)")
# plt.title("Cyclone Wind Field (Holland Model)")

# plt.colorbar(label="Wind Speed")

# plt.axis('equal')
# plt.show()

# plt.quiver(X,Y,u,v)
# plt.savefig('holland_wind_field.png')
# plt.show()

import numpy as np, matplotlib.pyplot as plt, imageio, os

def category(v):
    return ('TD' if v < 17 else 'TS' if v < 33 else 'Cat 1' if v < 43
            else 'Cat 2' if v < 50 else 'Cat 3' if v < 59
            else 'Cat 4' if v < 70 else 'Cat 5')

RMAX, B = 20, 1.2
vmaxs = np.linspace(18, 70, 60)

os.makedirs('vmax_frames', exist_ok=True)

for i, vmax in enumerate(vmaxs):

    # u, v = holland_wind_field(
    #     lon2d, lat2d, lat0, lon0_plot,
    #     vmax=vmax, Rmax=RMAX, B=B
    # )

    ws = np.sqrt(u**2 + v**2)

    fig, ax = plt.subplots(figsize=(8,6))

    cf = ax.contourf(
        0, 0, ws,
        levels=np.arange(0, 76, 2),
        cmap='turbo', extend='max'
    )

    ax.streamplot(0, 0, u, v,
                  density=2, color='k', linewidth=0.8)

    ax.plot(0, 0, 'k*', ms=12)

    ax.text(
        0.02, 0.98,
        f'Rmax={RMAX} km\nB={B}',
        transform=ax.transAxes, va='top',
        bbox=dict(fc='white', alpha=0.9)
    )

    ax.set_title(
        f'Vmax={vmax:.1f} m/s ({category(vmax)})'
    )

    plt.colorbar(cf, ax=ax, label='Wind speed (m/s)')

    fig.savefig(f'vmax_frames/f_{i:03d}.png',
                dpi=120, bbox_inches='tight')
    plt.close()

imageio.mimsave(
    'vmax_sweep.gif',
    [imageio.imread(f'vmax_frames/f_{i:03d}.png')
     for i in range(len(vmaxs))],
    fps=10
)