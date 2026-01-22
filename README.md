# Numerical PDE Solvers — Heat Equation (1D & 2D)

This project implements numerical solutions of the heat equation using finite difference methods in one and two spatial dimensions, with both explicit and implicit time-stepping schemes.

The repository is designed to demonstrate a solid understanding of:

- Discretisation of partial differential equations  
- Explicit vs implicit finite difference methods  
- Stability and convergence (CFL condition)  
- Matrix-based formulations of PDE solvers  
- Scientific computing and visualisation in Python  

This work was completed as part of my preparation for MSc programmes in Applied Mathematics / Scientific Computing.

---

## Mathematical Model

The heat equation solved throughout the project is

∂u/∂t = α ∇²u

where:
- u(x,t) or u(x,y,t) is the temperature field
- α is the thermal diffusivity constant

Both one-dimensional and two-dimensional versions are considered, with Dirichlet and Neumann boundary conditions.

---

## Implemented Solvers

### 1D Heat Equation

- **Forward Euler (Explicit)**
  - Finite difference discretisation in space and time
  - Conditionally stable (CFL restriction)
  - Numerical instability demonstrated when stability conditions are violated

- **Crank–Nicolson (Implicit)**
  - Second-order accurate in time
  - Unconditionally stable for the heat equation
  - Implemented by solving a tridiagonal linear system

- **Explicit vs Implicit Comparison**
  - Same grid, timestep, and initial condition
  - CFL condition intentionally violated to show:
    - Instability of the explicit scheme
    - Stability of the implicit scheme

---

### 2D Heat Equation

- **Explicit Forward Euler (2D)**
  - Five-point stencil approximation of the Laplacian
  - Stability-restricted timestep
  - Time evolution visualised via 2D heat maps

- **Implicit Crank–Nicolson (2D)**
  - Matrix formulation of the 2D Laplacian
  - Sparse matrices constructed using Kronecker products
  - Linear systems solved using sparse solvers from SciPy

---

## Numerical Techniques Used

- Finite difference discretisation of spatial derivatives  
- CFL stability analysis  
- Sparse matrix construction (`scipy.sparse`)  
- Kronecker products for multidimensional operators  
- Direct linear solvers for implicit schemes  
- Snapshot-based visualisation of time evolution  

---

## Initial Conditions

The solvers support multiple initial temperature distributions, including:

- Sine profiles  
- Gaussian distributions  
- Step functions  
- Hot-spot / hot-centre configurations  
- Constant initial fields  

---

## Boundary Conditions

- **Dirichlet**: fixed temperature at the boundaries  
- **Neumann**: zero-flux (insulated) boundaries  

Boundary conditions are applied consistently after each timestep.

---

## Visualisation

- 1D temperature profiles plotted over space  
- 2D temperature fields visualised using heat maps  
- Configurable snapshots to track time evolution  

---

## Project Structure

Each solver is implemented as a standalone, readable script:

- 1D explicit solver  
- 1D implicit (Crank–Nicolson) solver  
- 2D explicit solver  
- 2D implicit (Crank–Nicolson) solver  
- Explicit vs implicit comparison script  

All code is written in pure Python using:
- NumPy  
- Matplotlib  
- SciPy (sparse linear algebra)

---

## Key Learning Outcomes

Through this project I developed:

- A strong practical understanding of numerical PDE solvers  
- Insight into the stability limitations of explicit schemes  
- Experience translating mathematical formulations into efficient code  
- Familiarity with matrix-based and sparse linear algebra methods  
- Good scientific coding and modular design practices  

---

## Possible Extensions

- Adaptive time-stepping  
- Non-uniform spatial grids  
- Alternative implicit methods (e.g. ADI schemes)  
- Extension to other PDEs such as the wave equation or reaction–diffusion systems  

---

## Author

**Rafael Marquez**  
Undergraduate Mathematics  
Birkbeck, University of London

