import ufl
import numpy
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem, nls, log, io

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with CG Element
V = fem.FunctionSpace(domain, ("P", 1))

#Define boundary conditions
def bc_fun(x):
    return 1 + x[0] + 2*x[1]

u_D = fem.Function(V)
u_D.interpolate(bc_fun)
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: numpy.full(x.shape[1], True, dtype=bool))
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets))

#Define problem
def q(u):
    return 1 + u**2

x = ufl.SpatialCoordinate(domain)
u_ufl = 1 + x[0] + 2*x[1]
f = - ufl.div(q(u_ufl)*ufl.grad(u_ufl))

uh = fem.Function(V)
v = ufl.TestFunction(V)
F = q(uh)*ufl.dot(ufl.grad(uh), ufl.grad(v))*ufl.dx - f*v*ufl.dx

problem = fem.petsc.NonlinearProblem(F, uh, bcs=[bc])

#Define solver
solver = nls.petsc.NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-6
solver.report = True

#Set linear solver settings
ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "cg"
opts[f"{option_prefix}pc_type"] = "gamg"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()


#Solve
log.set_log_level(log.LogLevel.INFO)
n, converged = solver.solve(uh)

#Output XDMF
xdmf = io.XDMFFile(domain.comm,"nonlin_poisson_ana.xdmf","w")
xdmf.write_mesh(domain)
xdmf.write_function(u_D)