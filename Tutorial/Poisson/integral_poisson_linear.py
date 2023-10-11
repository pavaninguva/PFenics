from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from petsc4py.PETSc import ScalarType

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with Lagrange element
V = fem.FunctionSpace(domain, ("P", 1))
R = fem.FunctionSpace(domain, ("R", 0))

#Create BC
uD = ScalarType(0.0)

#Identify locations of BCs
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)

boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs, V=V)