# Pressure holland Model, gradient wind equation from 'Holland 1998'

import numpy as np
import matplotlib.pyplot as plt

N = 200
x = np.linspace(-100, 100, N)  # km
y = np.linspace(-100, 100, N)

X, Y = np.meshgrid(x, y)

# Parameters
p_n = 101325   # ambient pressure (Pa) - standard atmosphere
p_c = 95000    # central pressure (Pa) - adjust for storm intensity
rho = 1.15     # air density (kg/m^3)
f   = 1e-4     # Coriolis parameter (s^-1) - typical midlatitude value
A   = 20**1.2  # scaling parameter (related to R_max and B)
B   = 1.2      # Holland B-parameter

def V_g(r):
    r = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)
    pressure_term = (A * B * (p_n - p_c) * np.exp(-A / r**B)) / (rho * r**B)
    coriolis_term = (r**2 * f**2) / 4
    V = np.sqrt(pressure_term + coriolis_term) - (r * f) / 2

    u = -V * np.sin(theta)
    v = V * np.cos(theta)
    return (u, v)
#V[r==0] = 0

(u, v) = wind_field(X, Y)
speed = np.sqrt(u**2 + v**2)

plt.figure(figsize=(8,8))

plt.contourf(X, Y, speed, levels=50)  # Wind speed background

plt.streamplot(X, Y, u, v, color='white', density=2)  # Streamlines

plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Gradient wind velocity (Holland Model)")

plt.colorbar(label="Wind Speed")

plt.axis('equal')
plt.show()

plt.quiver(X,Y,u,v)
plt.show()
print(A**(1/B))