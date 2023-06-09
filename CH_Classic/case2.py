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
tend = 20.0
nsteps = 2e3
dt = (tend-t)/nsteps

#Model parameters
kappa = 1e-2
M = 1.0

#Define mesh
nx = ny = 200
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0.0,0.0],[2*np.pi, 2*np.pi]], [nx,ny])

#Create output file
xdmf = io.XDMFFile(domain.comm,"CH2_2.xdmf","w")
xdmf.write_mesh(domain)

#Define model
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
ME = fem.FunctionSpace(domain,P1*P1)
q,v = ufl.TestFunctions(ME)

#Define solution variables and split 
u = fem.Function(ME)
u0 = fem.Function(ME)

c,mu = ufl.split(u)
c0,mu0 = ufl.split(u0)

#Define Equation
# c0 = ufl.variable(c0)
# f0 = 0.25*(1-c0**2)**2
# dfdc = ufl.diff(f0,c0)
c = ufl.variable(c)
f = 0.25*(1-c**2)**2
dfdc = ufl.diff(f,c)

F0 = inner(c,q)*dx - inner(c0,q)*dx + M*dt*inner(grad(mu),grad(q))*dx
F1 = inner(mu,v)*dx - inner(dfdc,v)*dx - kappa*inner(grad(c),grad(v))*dx
F = F0 + F1

#Initial conditions
# def initial_condition(x):
#     #Unpack x and y 
#     x_ = x[0]
#     y_ = x[1]

#     f = 0.05*(
#         np.cos(3*x_)*np.cos(4*y_) + (np.cos(3*x_)*np.cos(4*y_))**2 +
#         np.cos(x_ - 5*y_)*np.cos(2*x_ - y_))
#     return f

def initial_condition(x):
    # values = np.zeros((2,x.shape[1]))
    values = 0.0 + 0.02*(0.5-np.random.rand(x.shape[1]))
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
# opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()


#Introduce skipping for output to output only every nth step
stride = 50
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