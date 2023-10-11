from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from petsc4py.PETSc import ScalarType
from dolfinx_mpc import LinearProblem, MultiPointConstraint
import numpy as np

"""
Implementing a modified case from
https://ammar-hakim.org/sj/je/je1/je1-periodic-poisson.html
"""

#Create mesh
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0.0,0.0],[np.pi, np.pi]], [50,50])

#Create function space with Lagrange element
V = fem.FunctionSpace(domain, ("P", 1))

#Create Periodic BC for whole domain
#Define periodic BC
def periodic_bc(x):
    return np.logical_or(np.isclose(x[0],np.pi),np.isclose(x[1],np.pi))

def periodic_relation(x):
    out_x = np.zeros(x.shape)
    out_x[0] = np.pi - x[0]
    out_x[1] = np.pi - x[1]
    out_x[2] = x[2]
    return out_x

mpc = MultiPointConstraint(V)
mpc.create_periodic_constraint_geometrical(V, periodic_bc, periodic_relation,[])
mpc.finalize()


#Define equation
#Define test function and solution variable
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Define source term
x = ufl.SpatialCoordinate(domain)
f = ufl.exp(-((x[0]-(np.pi/2))**2 + (x[1]-(np.pi/2))**2)/0.02) 

#Declare equation
a = ufl.dot(ufl.grad(u),ufl.grad(v))*ufl.dx
L = f*v*ufl.dx

problem=LinearProblem(a,L,mpc)

#Solve
uh = problem.solve()

#Output
xdmf = io.XDMFFile(domain.comm,"periodic_poisson2.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
xdmf.write_function(uh)
