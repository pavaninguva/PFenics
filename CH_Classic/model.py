import ufl
from dolfinx import log, plot
from dolfinx.fem import Function, FunctionSpace
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.io import XDMFFile
from dolfinx.mesh import CellType, create_unit_square
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import numpy as np

from mpi4py import MPI
from petsc4py import PETSc

import pyvistaqt as pvqt
import pyvista as pv

"""
Introduction

The classic CH equation is given as

dc/dt = \nabla \dot (M \nabla \mu)

where:
c: Concentration / composition 
t: Time
M: Mobility, which is assumed constant in this case
mu: Chemical potential

To identify mu, first consider the homogeneous free energy of the system, f:

f(c) = f_mix + \sum_{i} (c_i f*_i),

where f_mix is the free energy of mixing and f*_i is the free energy of the pure component i

For the purposes of the CH equation, it does not see f*_i as it eliminated in the derivation. 
So we just have to consider f_mix, which in the case of the classic CH equation is

f_mix = A*c^2*(1-c)^2,

where A is a parameter that controls the "depth" of the free energy well

We can abstract this a little bit more to account for different equilibrium compositions of the phases,

f_mix = A*(c-\alpha)^2*(c-\beta)^2,

where alpha and beta are the compositions of the two phases.

To get the total free energy of the system, we include the nonhomogeneous term

F = \int_v f_homo + f_non-homo dv = \int_v f_mix + (k/2)*(\nabla c)^2 dv,

where k is the gradient energy parameter. Computing the variation derivative:

mu_i = dF/dc_i - \nabla \dot dF/d(\nabla c_i) 

The actual mu in the CH equation is the difference between mu_1 and mu_2 

mu = df_mix/dc + k*\nabla^2 c

We thus end up with our final equation system which is a 4th order nonlinear PDE

It is commonly expressed as a set of coupled 2nd order PDEs:

dc/dt = M*\nabla^2 \mu
\mu = df_mix/dc + k*\nabla^2 c 

References:

1. Lee, Dongsun, et al. "Physical, mathematical, and numerical derivations of the Cahn–Hilliard equation." Computational Materials Science 81 (2014): 216-225.

2. Inguva, Pavan K., et al. "Continuum-scale modelling of polymer blends using the Cahn–Hilliard equation: transport and thermodynamics." Soft Matter 17.23 (2021): 5645-5665.

3. Nauman, E. Bruce, and David Qiwei He. "Nonlinear diffusion and phase separation." Chemical Engineering Science 56.6 (2001): 1999-2018.

"""

"""
Define the model
"""
#Model parameters
k = 1.0e-2  # gradient energy parameter
alpha = 0.0 # Compositions of phase 1
beta = 1.0 #Composition of phase 2

#Timestepping parameters
dt = 5.0e-06  # time step
theta = 0.5  # time stepping family, e.g. theta=1 -> backward Euler, theta=0.5 -> Crank-Nicholson

#Define mesh and finite element functions
msh = create_unit_square(MPI.COMM_WORLD, 96, 96, CellType.triangle)
P1 = ufl.FiniteElement("Lagrange", msh.ufl_cell(), 1)
ME = FunctionSpace(msh, P1 * P1)

q, v = ufl.TestFunctions(ME)

u = Function(ME)  # current solution
u0 = Function(ME)  # solution from previous converged step

# Split mixed functions
c, mu = ufl.split(u)
c0, mu0 = ufl.split(u0)

# Zero u
u.x.array[:] = 0.0

# Interpolate initial condition
u.sub(0).interpolate(lambda x: 0.63 + 0.02 * (0.5 - np.random.rand(x.shape[1])))
u.x.scatter_forward()

# Compute the chemical potential df/dc
c = ufl.variable(c)
f = 100 * (c-alpha)**2 * (c-beta)**2
dfdc = ufl.diff(f, c)

# mu_(n+theta)
mu_mid = (1.0 - theta) * mu0 + theta * mu

# Weak statement of the equations
F0 = inner(c, q) * dx - inner(c0, q) * dx + dt * inner(grad(mu_mid), grad(q)) * dx
F1 = inner(mu, v) * dx - inner(dfdc, v) * dx - k * inner(grad(c), grad(v)) * dx
#Couple the equations
F = F0 + F1

"""
Set up the solver
"""

# Create nonlinear problem and Newton solver
problem = NonlinearProblem(F, u)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-6

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()


"""
Timestep
"""
T = 10*dt
t=0.0
while t < T:
    t += dt
    res = solver.solve(u)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    u0.x.array[:] = u.x.array


"""
Plot
"""

V0, dofs = ME.sub(0).collapse()
topology, cell_types, x = plot.create_vtk_mesh(V0)
#Generate mesh object for pyvista
grid = pv.UnstructuredGrid(topology,cell_types,x)

#Get last c values from last timestep
grid.point_data["c"] = u.x.array[dofs].real
grid.set_active_scalars("c")

pv.OFF_SCREEN = True
pv.plot(grid,show_edges=True, screenshot = "c.png")




