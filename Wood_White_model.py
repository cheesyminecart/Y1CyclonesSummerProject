import numpy as np 
import matplotlib.pyplot as plt 
# This is a new model that uses 5 parameters, V_X, R_X, k, n and l
# V_X is the maximum velocity, or the velocity at the boundary (either that of the eye/eye wall)
# R_X is the radial distance at which said maximum velocity is attained 
# k is ..... (please write about this later)
# n is ..... (please write about this later)
# l is ..... (please write abou this later)
# f is the coriolis parameter 

# These parameters are defined for each region.  For now, we consider 2 regions, the eye and the eye wall.  In this model, we will attempt to visualise these 2 regions as a vector field 

# for conversion into cartesian later 

x = np.linspace(-10, 10, 40)
y = np.linspace(-10, 10, 40)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)


theta = np.arctan2(Y, X)


# defining paramters for region 1, the eye 
R_1X = 3
n_1 = 2
l_1 = 0.2
f = 0.0001 
k_1 = 1.5 
V_1X = 2
p_1 = r/R_1X


R_1CX = V_1X/(f * R_1X)

l_1 = 0.5 * ((1/R_1CX) + 2 * np.sqrt(1 + 1/(4 * (R_1CX ** 2))))

phi_1 = ((n_1 ** l_1) * (p_1 ** k_1)) * (1/(n_1 - k_1 + k_1 * p_1 ** (n_1/l_1))) ** l_1 

V_1 = (V_1X * l_1) * phi_1 

# defining parameters for region 2 
R_2X = 5
n_2 = 2
l_2 = 0.2
f = 0.0001 
k_2 = 1.5
V_2X = 2
p_2 = r/R_2X


R_2CX = V_2X/f * R_2X

l_2 = 0.5 * ((1/R_2CX) + 2 * np.sqrt(1 + 1/(4 * R_2CX ** 2)))

phi_2 = ((n_2 ** l_2) * (p_2 ** k_2)) * (1/(n_2 - k_2 + k_2 * p_2 ** (n_2/l_2))) ** l_2

V_2 = (V_2X * l_2) * phi_2 

# combining the two...


V = np.where(r <= R_1X, V_1, V_2)

U = -V * np.sin(theta)
W =  V * np.cos(theta)

plt.quiver(
    X, Y,
    U, W,
    V,                  # colour by magnitude
    scale= 2,            # smaller = bigger arrows
    scale_units='xy',
    angles='xy'
)

plt.xlim(-10, 10)
plt.ylim(-10,)

plt.xlabel("x")

plt.ylabel("y")

plt.title("Cyclone Vector Field")

plt.colorbar(label="Velocity Magnitude")

plt.axis("equal")

plt.show()
