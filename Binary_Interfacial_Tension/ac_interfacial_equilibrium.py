import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import ufl

import matplotlib.pyplot as plt

#Create mesh
nx = 50
domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(-2.0,2.0), nx=100)

#Specify basis functions
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
V = fem.FunctionSpace(domain,P1)
v = ufl.TestFunction(V)

#Define problem
c = fem.Function(V)
c.name="$c$" 

#Compute f'
f = 0.25*(c**2 -1)**2
dfdc = ufl.diff(f,c)
kappa = 0.005

F = kappa*inner(grad(c),grad(v))*dx + inner(dfdc,v)*dx

def ic(x):
    return np.tanh(x[0]/kappa)

#Apply IC
c.interpolate(ic)

#Setup solver

problem = NonlinearProblem(F, c)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.atol = 1e-8

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly"
opts[f"{option_prefix}pc_type"] = "lu"
ksp.setFromOptions()

sol = solver.solve(c)
print(sol[0])


#Plot
cells, types, x = plot.create_vtk_mesh(V)
print(x[:,0])
fig1, ax1 = plt.subplots(1,1,num=1)
#Plt numerical solution
ax1.plot(x[:,0],c.x.array.real)
ax1.plot(x[:,0], np.tanh(x[:,0]/kappa))

plt.show()