# Dynamical System Phase Portrait Generator
# This is a Generator for 2D dynamical systems
# The user inputs the function for dy/dt and dx/dt and the phase portrait is displaed

# Necessary Imports
import numpy as np
import matplotlib.pyplot as plt

# Initial Variables

Lx = 10        #Length in x direction
Nx = 11        #Number of spacial points on dimension x

# Generate Line
x_line = np.linspace(0,Lx,Nx)
zero_vector = np.zeros(Nx)

# dx/dt
def derivative_x(x_coord):
    # Here you input what dx/dt is
    dx = x_coord*(x_coord -2)
    return dx

# calculate the derivative of every point on the line
def derivative_vector(vector):
    for i in range(Nx):
        zero_vector[i] = derivative_x(vector[i])
    return zero_vector

gradient_vector = derivative_vector(x_line)
plt.plot(x_line, gradient_vector)
plt.show()
#print(derivative_vector(x_line))
#print(x_line)

    

