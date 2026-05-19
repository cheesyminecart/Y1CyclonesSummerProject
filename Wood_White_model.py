import numpy as np 
import matplotlib.pyplot as plt 
# This is a new model that uses 5 parameters, V_X, R_X, k, n and l
# V_X is the maximum velocity, or the velocity at the boundary (either that of the eye/eye wall)
# R_X is the radial distance at which said maximum velocity is attained 
# find out about k, n and l in the onenote!!!!
# f is the coriolis parameter 


# for conversion into cartesian later 

x = np.linspace(-10, 10, 40)
y = np.linspace(-10, 10, 40)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)


theta = np.arctan2(Y, X)


# defining paramters ()
R_1X = 3
n_1 = 4
l_1 = 0.01
k_1 = 2
f = 2. * 7.292115e-5 * np.sin(np.deg2rad(21.5))
V_1X = 9
p_1 = r/R_1X


R_1CX = V_1X/(f * R_1X)

phi_1 = ((n_1 ** l_1) * (p_1 ** k_1)) * (1/(n_1 - k_1 + k_1 * p_1 ** (n_1/l_1))) ** l_1 

V = (V_1X) * phi_1 

# we can define an equation for the coriolis parameter, see https://fabienmaussion.info/climate_system/week_04/01_Lesson_Wind-Derivatives-Integrals.html
# this will let us get more accurate values this parameter, since i have a program that gets us the HOUR ON HOUR COORDINATE POSITIONS OF THE CENTER OF CIRCULATION! (we remove 1 more source of error!)

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

plt.xlim(0, 9)
plt.ylim(0,8)

plt.xlabel("x")

plt.ylabel("y")

plt.title("Cyclone Vector Field")

plt.colorbar(label="Velocity Magnitude")

plt.axis("equal")
plt.show()

#plt.plot(r, V)
#plt.show()



