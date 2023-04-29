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
tend = 0.5
nsteps = 50
dt = (tend-t)/nsteps

#Define mesh
nx,ny = 50,50
domain = mesh.create_unit_square(MPI.COMM_WORLD,nx,ny,mesh.CellType.triangle)
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
V = fem.FunctionSpace(domain,P1)
v = ufl.TestFunction(V)

#Define initial condition

def initial_condition(x,a=0.1,r0=0.25):
    return np.tanh((r0 - np.sqrt(( x[0]-0.5)**2  + (x[1]-0.5)**2))/(np.sqrt(2)*a))

# u_n = fem.Function(V)
# u_n.name="u_n"
# u_n.interpolate(initial_condition)

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
kappa = 0.1

F = inner(c,v)*dx - inner(c0,v)*dx + kappa*dt*inner(grad(c),grad(v))*dx + inner(dfdc,v)*dx

#Setup solver
problem = NonlinearProblem(F, c)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.rtol = 1e-6

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



