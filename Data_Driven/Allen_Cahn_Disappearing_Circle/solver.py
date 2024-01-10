import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import ufl

#Define time parameters
t = 0.0
tend = 200
nsteps = 100
dt = (tend-t)/nsteps

#Define mesh
nx = ny = 100

domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0, 0],[2*np.pi, 2*np.pi]], [nx,ny])

# domain = mesh.create_unit_square(MPI.COMM_WORLD,nx,ny,mesh.CellType.triangle)
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),2)
V = fem.FunctionSpace(domain,P1)
v = ufl.TestFunction(V)

#Define initial condition

def initial_condition(x,a=0.1,r0=2.0):
    return np.tanh((np.sqrt(( x[0]-np.pi)**2  + (x[1]-np.pi)**2)-r0)/(np.sqrt(2)*a))


#no flux boundary conditions are the default

#Create output file
xdmf = io.XDMFFile(domain.comm,"ac.xdmf","w")
xdmf.write_mesh(domain)

#Define solution variable and interpolate initial solution
#Solution variable from previous timestep
c0 = fem.Function(V)
c0.name="c0"
c0.interpolate(initial_condition)

c = fem.Function(V)
c.name="$c$" 
c.interpolate(initial_condition)
xdmf.write_function(c,t)


#Compute f'
f = 0.25*(c**2 -1)**2
dfdc = ufl.diff(f,c)
kappa = 0.01

F = inner(c,v)*dx - inner(c0,v)*dx + kappa*dt*inner(grad(c),grad(v))*dx + dt*inner(dfdc,v)*dx

#Setup solver
problem = NonlinearProblem(F, c)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.atol = 1e-8

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "gmres"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()

#Timestepping

while t < tend:
    #Update t
    t += dt
    #Solve
    res = solver.solve(c)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    #Update c0
    c0.x.array[:] = c.x.array
    #Write to xdmf
    xdmf.write_function(c,t)
xdmf.close()