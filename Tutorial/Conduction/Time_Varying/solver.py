import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io
import ufl

#Timestepping parameters
t = 0 # Start time
T = 50.0 # Final time
num_steps = 500     
dt = T / num_steps # time step size
theta = 0.5

#Define mesh
nx, ny = 50, 50
domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([-2, -2]), np.array([2, 2])], 
                               [nx, ny], mesh.CellType.triangle)

#Define function space
V = fem.FunctionSpace(domain, ("CG", 1))


#Create homogeneous Dirichlet BCs
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)

#Create and interpolate initial condition
def initial_condition(x, a=1):
    return 10*np.exp(-a*(x[0]**2+x[1]**2))

u0 = fem.Function(V)
u0.interpolate(initial_condition)

#Define solution variable 
uh = fem.Function(V)
uh.interpolate(initial_condition)


#Define model equations
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, PETSc.ScalarType(0))
D=1.0
D_ = fem.Constant(domain, PETSc.ScalarType(D))


a = ufl.inner(u,v)*ufl.dx + dt*theta*D_*ufl.dot(ufl.grad(u), ufl.grad(v))*ufl.dx
L = u0*v*ufl.dx +  dt*f*v*ufl.dx - dt*(1-theta)*D_*ufl.dot(ufl.grad(u0),ufl.grad(v))*ufl.dx

#Build solver
problem = fem.petsc.LinearProblem(a,L, bcs=[bc],
                                  petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

#Create output file
xdmf = io.XDMFFile(domain.comm,"heat.xdmf","w")
xdmf.write_mesh(domain)
xdmf.write_function(uh,t)

#Timestep
while t < T - 1e-8:
    #Update t
    t += dt
    #Update D
    D_.value = D*np.exp(-2*t)
    #Solve
    uh = problem.solve()
    #Update u0
    u0.x.array[:] = uh.x.array
    #Output to xdmf
    xdmf.write_function(uh,t)