import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import spsolve


# Parameters

Lx = 1.0                # Length of x component
Ly = 1.0                # Length of y component
Nx = 50                 # Number of spatial points x
Ny = 50                 # Number of spatial points y
alpha = 0.1             # Thermal diffusivity constant

dx = Lx / (Nx - 1)      # Spacial step in x
dy = Ly / (Ny - 1)      # Spacial step in y

dt = dx**2 / (4 * alpha)        # Suitable timestep
t_target = 1.0
Nt = int(t_target / dt)         # Number of iterations calculated using the target time passed


# Initialize Grid

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')


# Initial conditions function: 4 states

def initialize_grid_2d(X, Y, initial_type="constant"):
    if initial_type == "constant":
        return np.ones_like(X)

    elif initial_type == "hot_center":
        u = np.zeros_like(X)
        u[(X - 0.5)**2 + (Y - 0.5)**2 < 0.05] = 1.0
        return u

    elif initial_type == "gaussian":
        return np.exp(-50 * ((X - 0.5)**2 + (Y - 0.5)**2))

    elif initial_type == "sine":
        return np.sin(np.pi * X) * np.sin(np.pi * Y)

    else:
        raise ValueError("Unknown initial condition")

u0 = initialize_grid_2d(X, Y, "constant")


# Boundary conditions (2D) Dirichlet

def apply_dirichlet(u):
    u[0, :] = 0
    u[-1, :] = 0
    u[:, 0] = 0
    u[:, -1] = 0
    return u


# Laplacian construction: Unique to 2D implicit solution


# Function to build the matrix representation of the second derivative d^2/dx^2

def second_derivative_matrix(N, d):
    main = -2 * np.ones(N)
    off = np.ones(N - 1)
    return diags([off, main, off], [-1, 0, 1]) / d**2


# Using the constructed second derivatives to create the full operation to calculate the Laplacian ∇^2
def laplacian_2d(Nx, Ny, dx, dy):
    Dxx = second_derivative_matrix(Nx, dx)
    Dyy = second_derivative_matrix(Ny, dy)
    Ix = identity(Nx)
    Iy = identity(Ny)
    return kron(Iy, Dxx) + kron(Dyy, Ix)

Lmat = laplacian_2d(Nx, Ny, dx, dy)
I = identity(Nx * Ny)

# Creating Crank Nicolson Matrices A and B

A = I - 0.5 * alpha * dt * Lmat
B = I + 0.5 * alpha * dt * Lmat


# Crank–Nicolson step (VECTOR): Solves a linear system to find the field at t = n+1

def crank_nicolson_step(u_vec):
    rhs = B @ u_vec
    return spsolve(A, rhs)


# Plotting function

def plot_2d(u):
    plt.imshow(u, origin="lower", cmap="hot")
    plt.colorbar()
    plt.show()


# Run simulation function:
# Includes the time steps in as well as num_snaps to take a number of snapshots of the simulation

def run_simulation(u0, num_snaps=10):
    u = u0.copy()
    u_vec = u.flatten()                 # flattens the matrix into a vector to perform linear algebra
    snap_every = Nt // num_snaps

    for n in range(Nt):
        u_vec = crank_nicolson_step(u_vec)      #perform a crank nicolson step on this vector flattened form

        # Reshape, apply boundary conditions and flatten again
        u = u_vec.reshape((Nx, Ny))
        u = apply_dirichlet(u)
        u_vec = u.flatten()

        if n % snap_every == 0:
            plot_2d(u)


# Run the simulation

run_simulation(u0, num_snaps=10)

