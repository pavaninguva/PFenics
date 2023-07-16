from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with CG Element
V = fem.FunctionSpace(domain, ("P", 1))

#Create BCs

#Bottom, 
facets = mesh.locate_entities_boundary(domain, domain.topology.dim-1, 
                                        marker = lambda x: np.isclose(x[1],0.0)
                                        )
dof = fem.locate_dofs_topological(V=V, entity_dim=domain.topology.dim-1,entities=facets)
bc = fem.dirichletbc(value=ScalarType(0.0),dofs=dof, V=V)

# Left, Right
facets1 = mesh.locate_entities_boundary(domain, domain.topology.dim-1, 
                                        marker = lambda x: np.logical_or(np.isclose(x[0],0.0),
                                                                         np.isclose(x[0],1.0),
                                        ))

dof1 = fem.locate_dofs_topological(V=V, entity_dim=domain.topology.dim-1,entities=facets1)
bc1 = fem.dirichletbc(value=ScalarType(0.0),dofs=dof1, V=V)

#Top
facets2 = mesh.locate_entities_boundary(domain, domain.topology.dim-1, 
                                        marker = lambda x: np.isclose(x[1],1.0)
                                        )
dof2 = fem.locate_dofs_topological(V=V, entity_dim=domain.topology.dim-1,entities=facets2)
bc2 = fem.dirichletbc(value=ScalarType(1.0),dofs=dof2, V=V)

#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define source term
f = fem.Constant(domain, ScalarType(0))

#Declare equation
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

#Solve with direct solver
problem = fem.petsc.LinearProblem(a, L, bcs=[bc,bc1,bc2], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

#Output to xdmf
xdmf = io.XDMFFile(domain.comm,"poisson2.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh)

