from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from petsc4py.PETSc import ScalarType


#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with Lagrange element
V = fem.FunctionSpace(domain, ("P", 1))

#Create BC
uD = ScalarType(0.0)

#Identify locations of BCs
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)

boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs, V=V)

#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define source term
f = fem.Constant(domain, ScalarType(1))

#Declare equation
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

#Solve with direct solver
problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

#Output to xdmf
xdmf = io.XDMFFile(domain.comm,"poisson1.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh,0)

