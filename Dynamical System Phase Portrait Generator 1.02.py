# Dynamical System Phase Portrait Generator
# This is a Generator for 2D dynamical systems
# The user inputs the function for dy/dt and dx/dt and the phase portrait is displaed

# Necessary Imports
import numpy as np
import matplotlib.pyplot as plt

# Initial Variables

Lx = 5        #Length in x direction
Nx = 41        #Number of spacial points in x direction
Ny = 41        #Number of spacial points in y direction

# Generate Line
x_coord = np.linspace(-Lx,Lx,Nx)
y_coord = np.linspace(-Lx,Lx,Nx)

X, Y = np.meshgrid(x_coord, y_coord)


# dx/dt
def derivative(x, y):
    # Here you input dx/dt and dy/dt
    dx = x*(4-x-y)
    dy = y*(3-x-2*y)
    return dx, dy


# calculate the derivative of every point on the grid
U, V = derivative(X,Y)



#scale all the vectors so they have unit length
scalar = np.sqrt(U**2 + V**2)

# avoid division by zero
scalar[scalar == 0] = 1

U_unit = U / scalar
V_unit = V / scalar


plt.quiver(X,Y,U_unit, V_unit, color = 'red')
plt.axhline(0,color = 'black', linewidth=1)
plt.axvline(0,color = 'black', linewidth=1)
plt.show()


    

