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

plt.figure(figsize=(8,8))

plt.contourf(X, Y, speed, levels=50)  # Wind speed background

plt.streamplot(X, Y, u, v, color='white', density=2)  # Streamlines

plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Cyclone Wind Field (Holland Model)")

plt.colorbar(label="Wind Speed")

plt.axis('equal')
plt.show()

plt.quiver(X,Y,u,v)
plt.show()
