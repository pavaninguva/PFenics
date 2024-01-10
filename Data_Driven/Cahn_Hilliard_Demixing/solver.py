import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, ln
import ufl


#Define time parameters
t = 0.0
tend = 1000
dt = 1.0

#Parameters
Nchi = 3
kappa = (2/3)*Nchi

#Define mesh
nx = ny = 128
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0.0,0.0],[128.0, 128.0]], [nx,ny])

#Create output file
xdmf = io.XDMFFile(domain.comm,"spinodal.xdmf","w")
xdmf.write_mesh(domain)

#Define test functions
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
ME = fem.FunctionSpace(domain,P1*P1)
q,v = ufl.TestFunctions(ME)

#Define solution variables and split 
u = fem.Function(ME)
u0 = fem.Function(ME)

c,mu = ufl.split(u)
c0,mu0 = ufl.split(u0)

#Define chemical potential
c = ufl.variable(c)
f = c*ln(c) + (1-c)*ln(1-c) + Nchi*c*(1-c)
dfdc = ufl.diff(f,c)

F0 = inner(c,q)*dx - inner(c0,q)*dx + (c*(1-c))*dt*inner(grad(mu),grad(q))*dx
F1 = inner(mu,v)*dx - inner(dfdc,v)*dx - kappa*inner(grad(c),grad(v))*dx
F = F0 + F1

def initial_condition(x):
    values = 0.5 + 0.02*(0.5-np.random.rand(x.shape[1]))
    return values

#Apply IC
u.x.array[:] = 0.0
u.sub(0).interpolate(initial_condition)
u.x.scatter_forward()
c = u.sub(0)
u0.x.array[:] = u.x.array
xdmf.write_function(c,t)

#Setup solver
problem = NonlinearProblem(F, u)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.atol = 1e-8
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly"
opts[f"{option_prefix}pc_type"] = "lu"
# opts[f"{option_prefix}pc_hypre_type"] = "euclid"
# opts[f"{option_prefix}ksp_max_it"] = "1000"
ksp.setFromOptions()

#Introduce skipping for output to output only every nth step
stride = 1
counter = 0

#Timestepping
while t < tend:
    #Update t
    t += dt
    #Solve
    res = solver.solve(u)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    #Update u0
    u0.x.array[:] = u.x.array
    counter = counter +1
    if counter % stride == 0:
        xdmf.write_function(c,t)
xdmf.close()