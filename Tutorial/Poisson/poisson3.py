from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from ufl import cos
import numpy as np
from petsc4py.PETSc import ScalarType


#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 50, 50, mesh.CellType.quadrilateral)

#Create function space with CG Element
V = fem.FunctionSpace(domain, ("P", 2))

#No BCs are required

#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define source term
x = ufl.SpatialCoordinate(domain)
f = (3 + 2*np.pi**2)*cos(np.pi*x[0])*cos(np.pi*x[1])

#Declare equation
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx + 3*u*v*ufl.dx  
L = f * v * ufl.dx

#Solve with direct solver
problem = fem.petsc.LinearProblem(a, L, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

#Output to xdmf
xdmf = io.XDMFFile(domain.comm,"poisson3.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh)