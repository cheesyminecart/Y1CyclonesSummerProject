import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import imageio.v2 as imageio
import os

# -----------------------------
# Grid and physics parameters
# -----------------------------
N = 200
x = np.linspace(-100, 100, N)  # km
y = np.linspace(-100, 100, N)
X, Y = np.meshgrid(x, y)

RMAX, B = 20.0, 1.2
vmaxs = np.linspace(18, 70, 60)  # m/s

# -----------------------------
# Category label (Saffir–Simpson-like)
# -----------------------------
def category_label(v):
    if v < 17:
        return "Tropical Depression"
    elif v < 33:
        return "Tropical Storm"
    elif v < 43:
        return "Category 1"
    elif v < 50:
        return "Category 2"
    elif v < 59:
        return "Category 3"
    elif v < 70:
        return "Category 4"
    else:
        return "Category 5"

# -----------------------------
# Holland-like wind field
# -----------------------------
def wind_field(X, Y, Vmax, Rmax, B):
    r = np.hypot(X, Y)
    with np.errstate(divide='ignore', invalid='ignore'):
        term = (Rmax / r)**B
        V = Vmax * np.sqrt(term * np.exp(1.0 - term))
        V = np.where(np.isfinite(V), V, 0.0)  # set r=0 (and any inf/nan) to 0

    theta = np.arctan2(Y, X)
    u = -V * np.sin(theta)  # tangential, cyclonic
    v =  V * np.cos(theta)
    u = np.where(r == 0, 0.0, u)
    v = np.where(r == 0, 0.0, v)
    return u, v, V

# -----------------------------
# Streamline seed points (rings)
# -----------------------------
def circular_seed_points(radii, n_per_ring=64):
    pts = []
    for r in radii:
        ang = np.linspace(0, 2*np.pi, n_per_ring, endpoint=False)
        pts.append(np.column_stack([r*np.cos(ang), r*np.sin(ang)]))
    return np.vstack(pts)

seed_radii = np.arange(5, 100, 5)  # km
seeds = circular_seed_points(seed_radii, n_per_ring=48)

# -----------------------------
# Plot/GIF settings
# -----------------------------
os.makedirs("vmax_frames", exist_ok=True)
vmin, vmax = 0.0, float(vmaxs.max())  # fix color range across frames
norm = Normalize(vmin=vmin, vmax=vmax)
levels = np.linspace(vmin, vmax, 30)  # smooth color gradations
cmap = "viridis"                      # close to your reference look

# -----------------------------
# Generate frames
# -----------------------------
for i, vmax_frame in enumerate(vmaxs):
    u, v, ws = wind_field(X, Y, Vmax=float(vmax_frame), Rmax=RMAX, B=B)

    fig, ax = plt.subplots(figsize=(7.5, 7.5))  # square plot

    # Filled contours of wind speed
    cf = ax.contourf(X, Y, ws, levels=levels, cmap=cmap, norm=norm, extend="max")

    # White streamlines, circular look via seeded rings
    ax.streamplot(
        X, Y, u, v,
        color="white",
        density=1.5,       # base density
        linewidth=0.9,
        arrowsize=0.9,
        start_points=seeds
    )

    # Center marker
    ax.plot(0, 0, marker="o", ms=8, color="purple", markeredgecolor="white", markeredgewidth=0.8)

    # Axes styling
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")

    # Title with category
    ax.set_title(f"Vmax = {vmax_frame:.1f} m/s  ({category_label(vmax_frame)})")

    # Colorbar
    cbar = plt.colorbar(cf, ax=ax, label="Wind Speed (m/s)")
    cbar.ax.tick_params(length=0)

    fig.savefig(f"vmax_frames/f_{i:03d}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
# Write GIF
# -----------------------------
imageio.mimsave(
    "vmax_sweep.gif",
    [imageio.imread(f"vmax_frames/f_{i:03d}.png") for i in range(len(vmaxs))],
    fps=10,
)

print("Saved vmax_sweep.gif")