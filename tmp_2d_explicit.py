import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0            # length of square
Nx = 50            # number of spatial points in x direction
Ny = 50            # number of spatial points in y direction
alpha = 0.1       # thermal diffusivity
dx = L / (Nx - 1)  # spatial step in x
dy = L / (Ny - 1)  # spatial step in y

# Stability Condition:
dt = dx**2 / (4*alpha)  # safe dt

t_target = 1
Nt = int(t_target/dt)   #number of time steps defined by how much time you want to pass

# Establish r: value to represent rate of temperature spread

r = alpha * dt/(dx**2)


# Initialize grid function for 2D: 4 different starting states

def initialize_grid_2d(X, Y, initial_type="gaussian"):
    
    # initialize u(x,y,t) for t = 0

    if initial_type == "constant":
        u = np.ones_like(X)

    elif initial_type == "hot_center":
        u = np.zeros_like(X)
        cx, cy = 0.5, 0.5
        u[(X - cx)**2 + (Y - cy)**2 < 0.05] = 1.0

    elif initial_type == "gaussian":
        u = np.exp(-50 * ((X - 0.5)**2 + (Y - 0.5)**2))

    elif initial_type == "sine":
        u = np.sin(np.pi * X) * np.sin(np.pi * Y)

    else:
        raise ValueError("Unknown initial condition")

    return u

# Bounday conditions; 2 types

def apply_boundary(u,boundary_type):
    u_copy = u.copy()

    # Dirichlet: temperature = 0 at boundary
    if boundary_type == "dirichlet":
        for i in range(Nx):
            u_copy[i,0] = 0
            u_copy[i,Ny-1] = 0
        for j in range(Ny):
            u_copy[0,j] = 0
            u_copy[Nx-1,j] = 0


    # Neumann: temperature is defined by nearest points to the boundary
    elif boundary_type == "neumann":
        for i in range(0,Nx):
            u_copy[i,0] = u[i,1]
            u_copy[i,Ny-1] = u[i,Ny-2]
        for j in range(0,Ny):
            u_copy[0,j] = u[1,j]
            u_copy[Nx-1,j] = u[Nx-2,j]
    return(u_copy)

# Time step function

def time_step(u, r):
    u_new = u.copy()
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            u_new[i, j] = u[i, j] + r * (
                    u[i + 1, j] + u[i - 1, j] +
                    u[i, j + 1] + u[i, j - 1] -
                    4 * u[i, j]
            )
    return u_new

# Run simulation function
# takes numb number of timestamps and outputs them

def run_simulation(u,numb):
    u_new = u.copy()
    time_stamps = []

    for i in range(Nt):
        u_new = time_step(u_new, r)
        u_new = apply_boundary(u_new, "dirichlet")
        if (i % snapshot(numb)) == 0:
            time_stamps.append(u_new)
    for plt in time_stamps:
        plot_2d(plt)


# Plot function
def plot_2d(u_new):
    plt.imshow(u_new, origin='lower', cmap='hot')
    plt.colorbar()
    plt.show()

# Function to form the timestamps

def snapshot(numsnaps):

    numdiv = Nt // numsnaps
    return numdiv

# Full code execution: initialisaion and run simulation
x = np.linspace(0, L, Nx)
y = np.linspace(0, L, Ny)

X, Y = np.meshgrid(x, y, indexing='ij')
u = initialize_grid_2d(X, Y, initial_type="constant")


u_new = u.copy()
run_simulation(u_new, 10)


