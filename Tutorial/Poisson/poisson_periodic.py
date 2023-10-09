from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from petsc4py.PETSc import ScalarType
# import dolfinx_mpc.utils
from dolfinx_mpc import LinearProblem, MultiPointConstraint
import numpy as np


"""
Implementing the example from
https://fenicsproject.org/olddocs/dolfin/1.4.0/python/demo/documented/periodic/python/documentation.html
"""

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with Lagrange element
V = fem.FunctionSpace(domain, ("P", 1))

#Define Dirichlet BC
def dirichlet_bc(x):
    return np.logical_or(np.isclose(x[1],0),np.isclose(x[1],1))

#Create BC
uD = ScalarType(0.0)

#Create Dirichlet BC
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.locate_entities_boundary(domain, 1, dirichlet_bc)
boundary_dofs = fem.locate_dofs_topological(V,fdim,boundary_facets)
bc = fem.dirichletbc(uD,boundary_dofs,V)

#Define periodic BC
def periodic_bc(x):
    return np.isclose(x[0],1)

def periodic_relation(x):
    out_x = np.zeros(x.shape)
    out_x[0] = 1 - x[0]
    out_x[1] = x[1]
    out_x[2] = x[2]
    return out_x

mpc = MultiPointConstraint(V)
mpc.create_periodic_constraint_geometrical(V, periodic_bc, periodic_relation, [bc])
mpc.finalize()

#Define equation
#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define source term
x = ufl.SpatialCoordinate(domain)
f = x[0]*ufl.sin(5*np.pi*x[1]) + ufl.exp(-((x[0]-0.5)**2 + (x[1]-0.5)**2)/0.02)

#Declare equation
a = ufl.dot(ufl.grad(u),ufl.grad(v))*ufl.dx
L = f*v*ufl.dx

problem=LinearProblem(a,L,mpc,bcs=[bc])

#Solve
uh = problem.solve()

#Output
xdmf = io.XDMFFile(domain.comm,"periodic_poisson1.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh)

"""
Test with Neumann BCs instead of periodic
"""
# problem_ = fem.petsc.LinearProblem(a,L,bcs=[bc])
# uh_ = problem_.solve()
# #Output
# xdmf = io.XDMFFile(domain.comm,"periodic_poisson1_.xdmf","w")
# #Output mesh
# xdmf.write_mesh(domain)
# #Output solution
# xdmf.write_function(uh_)

