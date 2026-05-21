import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as ani

# This is a new model that uses 5 parameters, V_X, R_X, k, n and l
# V_X is the maximum velocity, or the velocity at the boundary (either that of the eye/eye wall)
# R_X is the radial distance at which said maximum velocity is attained 
# find out about k, n and l in the onenote!!!!
# f is the coriolis parameter 


# for conversion into cartesian later 

x = np.linspace(-10, 10, 20)
y = np.linspace(-10, 10, 20)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)


theta = np.arctan2(Y, X)


# defining paramters ()
R_1X = 3
n_1 = 4
l_1 = 0.01

# we can define an equation for the coriolis parameter, see https://fabienmaussion.info/climate_system/week_04/01_Lesson_Wind-Derivatives-Integrals.html
# this will let us get more accurate values this parameter, since i have a program that gets us the HOUR ON HOUR COORDINATE POSITIONS OF THE CENTER OF CIRCULATION! (we remove 1 more source of error!)
# JC 
#f = 2. * 7.292115e-5 * np.sin(np.deg2rad(21.5))

f = 2
k_1 = 2
V_1X = 2

# to avoid division by 0, let us define a small quantity epsilon... 

eps = 1e-10

# defining the velocity profile function... 


def V_profile(r):
    
    return(V_1X *((n_1 ** l_1) * (((r + eps)/R_1X) ** k_1)) * (1/(n_1 - k_1 + k_1 * ((r + eps)/R_1X) ** (n_1/l_1))) ** l_1 
)

# V at every point of the grid 

V_field = V_profile(r)

# Circulation at every r 
C_field = 2 * np.pi * r * V_field 

# defining initial circulation

C_0 = 2 * np.pi * R_1X * V_1X 

# Defining viscosity with sutherland's law 


u_0 = 0.5
T_0 = 300
S   = 110.4
T   = 310
p = 1.3

u = u_0 * ((T / T_0)**1.5) * (T_0 + S) / (T + S)
v = u / p  # kinematic viscosity

# defining the timesteps and putting the frames into an array
dt = 0.05
t_end = 4
time_steps = np.arange(0, t_end, dt)

# defining the derivative of C 

C_d = -v * V_profile(r) * 2 * np.pi * r 

C_current = C_field.copy()
frames_U = []
frames_W = []
frames_V = []
all_time = []
time = 0.0

while np.max(C_current) > 1e-6:
    V_current = C_current / (2 * np.pi * (r + eps))
    
    U = -V_current * np.sin(theta)
    W =  V_current * np.cos(theta)
    
    frames_U.append(U.copy())
    frames_W.append(W.copy())
    frames_V.append(V_current.copy())
    all_time.append(time)
    
    C_d = -v * V_current * 2 * np.pi * r
    C_current = C_current + dt * C_d
    time += dt

# Initial plot parameters
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Decay of a cyclone due to viscous torque ')

# Setting up the initial frame with colour grading 
quiv = ax.quiver(X, Y, frames_U[0], frames_W[0], frames_V[0], cmap='viridis', clim=[0, V_1X])
fig.colorbar(quiv, ax=ax, label='V')


# time counter (adapted from claude)
time_text = ax.text(0.03, 0.95, f't = {all_time[0]:.2f} s',
                    transform=ax.transAxes, fontsize=10,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


# frame update function 
def update(frame):
    quiv.set_UVC(frames_U[frame], frames_W[frame], frames_V[frame])
    time_text.set_text(f't = {all_time[frame]:.2f} ')
    return quiv, time_text



animation = ani.FuncAnimation(fig, update, frames=len(frames_U), interval=50)
plt.show()

V_max = [np.max(V) for V in frames_V]

fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.plot(all_time, V_max)
ax2.set_xlabel('Time')
ax2.set_ylabel('Max Velocity')
ax2.set_title('V_max over time ')
ax2.grid(True, alpha=0.3)
plt.show()