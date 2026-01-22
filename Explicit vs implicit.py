import numpy as np
import matplotlib.pyplot as plt
from fontTools.misc.cython import returns

# Parameters
L = 1.0            # length of rod
Nx = 50            # number of spatial points
alpha = 0.1       # thermal diffusivity
dx = L / (Nx - 1)  # spatial step
cfl = 0.6
# Note: cfl = 0.6 intentionally violates the explicit stability condition
# to demonstrate instability of Forward Euler



# Stability condition for Forward Euler
dt = cfl* dx**2 / alpha  # safe dt

#coefficient for crank nicolson
r = alpha * dt / dx**2

t_target = 1
Nt = int(t_target/dt)   #number of time steps defined by how much time you want to pass

if alpha <= 0 :
    raise ValueError("alpha must be > 0")

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

#boundary conditions:

def apply_boundary(u, boundary_type):
    #dirichlet both ends 0
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


def time_step_explicit(u,dx,dt, alpha):
    u_new = u.copy()
    for i in range(1, Nx-1):
        u_new[i] = u[i] + alpha * dt / dx ** 2 * (u[i + 1] - 2 * u[i] + u[i - 1])
    return u_new


#run simulation function
def run_simulation_explicit():
    temp = u.copy()
    for i in range(Nt):
        temp = time_step_explicit(temp,dx,dt, alpha)
        temp = apply_boundary(temp, boundary_type = "dirichlet")
    return temp



#run simulation function
def run_simulation_implicit():
    temp = u.copy()
    for i in range(Nt):
        temp = crank_nicolson_step(temp,A,B)
        temp = apply_boundary(temp, boundary_type = "dirichlet")
    return temp

#crank nicolson adaptation

#Create tridiagonal matrix
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

#Perform a crank nicolson step
def crank_nicolson_step(u,A,B):
    u_inner = u[1:-1]
    rhs = B @ u_inner

    u_new_inner = np.linalg.solve(A,rhs)
    u_new = u.copy()
    u_new[1:-1] = u_new_inner
    return u_new

def final_plot(v,condit):
    plt.plot(x, v, label=f't={Nt * dt:.2f}')
    plt.xlabel('x')
    plt.ylabel('Temperature u(x)')
    if condit == "F":
        plt.title('1D Heat Equation - Forward Euler')
    elif condit == "CN":
        plt.title('1D Heat Equation - Crank Nicolson')
    plt.legend()
    plt.show()

# Initial conditions: u(x,0) = sin(pi*x)
x = np.linspace(0, L, Nx)
u = initialize_grid(x, initial_type = "sine")

A,B = crank_nicolson_matrices(Nx,r)

# Run using functions created:

implicit_u = run_simulation_implicit()
final_plot(implicit_u,"CN")


explicit_u = run_simulation_explicit()
final_plot(explicit_u,"F")