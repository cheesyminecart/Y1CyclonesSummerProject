import numpy as np 
import matplotlib.pyplot as plt 


R = 5 # R, the radius of the eye, is constant 
V = 2 # V is the maximum velocity

# defining a buncha points in  space
x = np.linspace(-10, 10, 20)
y = np.linspace(-10, 10, 20)
X, Y = np.meshgrid(x, y)

# converting to polar coordinates 

r = np.sqrt(X**2 + Y**2)
phi = np.arctan2(Y, X)

# defining the magnitude of velocity
v_a = np.where(r <= R, V * (r / R), V * (R / r))

# converting BACK to cartesian (there should be an easier way to do this, please give suggestions!!)
U = -v_a * np.sin(phi)
W =  v_a * np.cos(phi)


plt.quiver(X, Y, U, W)
plt.show()






