from mpi4py import MPI
from dolfinx import mesh, fem, io
import numpy as np
import ufl 
from ufl import sin, exp
from petsc4py.PETSc import ScalarType

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with CG Element
V = fem.FunctionSpace(domain, ("P", 2))

#Define Dirichlet BC
facets = mesh.locate_entities_boundary(domain, domain.topology.dim-1, 
                                        marker = lambda x: np.logical_or(np.isclose(x[0],0.0),
                                                                         np.isclose(x[0],1.0),
                                        ))

dof = fem.locate_dofs_topological(V=V, entity_dim=domain.topology.dim-1,entities=facets)
bc = fem.dirichletbc(value=ScalarType(0.0),dofs=dof, V=V)

#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define Equation
x = ufl.SpatialCoordinate(domain)
f = 10 * exp(-((x[0] - 0.5)**2 + (x[1] - 0.5)**2) / 0.02)
g = sin(5 * x[0])

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = ufl.inner(f, v) * ufl.dx + ufl.inner(g, v) * ufl.ds

#Solve with direct solver
problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

#Output to xdmf
xdmf = io.XDMFFile(domain.comm,"poisson4.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh)