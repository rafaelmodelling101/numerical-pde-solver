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
def derivative_x(x, y):
    # Here you input what dx/dt is
    dx = x*(4-x-y)
    return dx

def derivative_y(x, y):
    # Here you input what dx/dt is
    dy = y*(3-x-2*y)
    return dy

# calculate the derivative of every point on the grid
x_vector = derivative_x(X,Y)
y_vector = derivative_y(X,Y)



plt.quiver(X,Y, x_vector, y_vector, color = 'red')
plt.axhline(0,color = 'black', linewidth=1)
plt.axvline(0,color = 'black', linewidth=1)
plt.show()


    

