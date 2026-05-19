
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


plt.figure()
plt.plot(time_steps, V_change)
plt.xlabel("Time (s)")
plt.ylabel("V (m/s)")
plt.title("V_mag vs Time")
plt.grid()
plt.show()
