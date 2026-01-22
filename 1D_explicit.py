import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0            # Length of rod
Nx = 50            # Number of spatial points
alpha = 0.1       # Thermal diffusivity
dx = L / (Nx - 1)  # Spatial step

# Stability condition for Forward Euler
dt = dx**2 / (2*alpha)  # Safe dt

t_target = 1
Nt = int(t_target/dt)   # Number of time steps defined by how much time you want to pass

# Functions utilised

# Initialize grid: 3 different starting states

def initialize_grid(x, initial_type):

    if initial_type == "sine":
        u = np.sin(np.pi * x)
    elif initial_type == "gaussian":
        u = np.exp(-100*(x-0.5)**2)
    elif initial_type == "step":
        u = np.where((x>= 0.4) & (x<= 0.6), 1.0, 0.0)
    else:
        raise ValueError(f"initial_type {initial_type} not recognized")
    return u

# Boundary conditions: two possible types

def apply_boundary(u, boundary_type):
    # Dirichlet both ends 0
    if boundary_type == "dirichlet":
        u[0] = 0
        u[-1]= 0

    # Neumman (fixed flux)
    elif boundary_type == "neumann":
        u[0] = u[1]
        u[-1] = u[-2]
    else:
        raise ValueError(f"boundary_type {boundary_type} not recognized")

    return u

# Time step for grid

def time_step(u,dx,dt, alpha):
    u_new = u.copy()
    for i in range(1, Nx-1):
        u_new[i] = u[i] + alpha * dt / dx ** 2 * (u[i + 1] - 2 * u[i] + u[i - 1])
    return u_new


# Run simulation function
def run_simulation():
    temp = u.copy()
    for i in range(Nt):
        temp = time_step(temp,dx,dt, alpha)
        temp = apply_boundary(temp, boundary_type = "dirichlet")
    return temp


# Plot the final temperature distribution
def plot_final(v):
    plt.plot(x, v, label=f't={Nt * dt:.2f}')
    plt.xlabel('x')
    plt.ylabel('Temperature u(x)')
    plt.title('1D Heat Equation - Forward Euler')
    plt.legend()
    plt.show()


# Initial conditions: u(x,0) = sin(pi*x)
# Initial conditions can be adjusted

x = np.linspace(0, L, Nx)
u = initialize_grid(x, initial_type = "sine")


# Run simulation
final_u = run_simulation()
plot_final(final_u)
