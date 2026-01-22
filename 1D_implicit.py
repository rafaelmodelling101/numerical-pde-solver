import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0            # length of rod
Nx = 50            # number of spatial points
alpha = 0.1       # thermal diffusivity
dx = L / (Nx - 1)  # spatial step
cfl = 0.4          #
if alpha <= 0 :
    raise ValueError("alpha must be > 0")


# Stability condition for Forward Euler
dt = cfl * dx**2 / alpha  # safe dt

t_target = 1
Nt = int(t_target/dt)   # Number of time steps defined by how much time you want to pass

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

# Boundary conditions: 2 types
def apply_boundary(u, boundary_type):
    #Dirichlet: both ends 0
    if boundary_type == "dirichlet":
        u[0] = 0
        u[-1]= 0

    #Neumman (fixed flux)
    elif boundary_type == "neumann":
        u[0] = u[1]
        u[-1] = u[-2]
    else:
        raise ValueError(f"boundary_type {boundary_type} not recognized")

    return u


# Run simulation function
def run_simulation():
    temp = u.copy()
    for i in range(Nt):
        temp = crank_nicolson_step(temp,A,B)
        temp = apply_boundary(temp, boundary_type = "dirichlet")
    return temp

# Crank nicolson adaptation

# Create tridiagonal matrix
def crank_nicolson_matrices(Nx, r):
    A = np.zeros((Nx-2, Nx-2))
    B = np.zeros((Nx-2, Nx-2))
    for i in range(Nx-2):
        A[i,i] = 1 + r
        B[i,i] = 1 - r
        if i > 0:
            A[i,i-1] = -r/2
            B[i,i-1] = r/2
        if i < Nx-3:
            A[i,i+1] = -r/2
            B[i,i+1] = r/2
    return A, B

# Perform a crank nicolson step
def crank_nicolson_step(u,A,B):
    u_inner = u[1:-1]
    rhs = B @ u_inner

    u_new_inner = np.linalg.solve(A,rhs)
    u_new = u.copy()
    u_new[1:-1] = u_new_inner
    return u_new

# Final plot function

def final_plot(v):
    plt.plot(x, v, label=f't={Nt * dt:.2f}')
    plt.xlabel('x')
    plt.ylabel('Temperature u(x)')
    plt.title('1D Heat Equation - Forward Euler')
    plt.legend()
    plt.show()

# Initialise line and grid
x = np.linspace(0, L, Nx)
u = initialize_grid(x, initial_type = "sine")

# Establish Crank Nicolson Matrix
r = alpha * dt / dx**2
A,B = crank_nicolson_matrices(Nx,r)

# Run simulation
final_u = run_simulation()
final_plot(final_u)
