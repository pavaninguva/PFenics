from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl

"""
https://fenicsproject.org/olddocs/dolfin/latest/python/demos/neumann-poisson/demo_neumann-poisson.py.html
"""

#Create mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10, mesh.CellType.triangle)

#Create function space with Lagrange multiplier
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
R0 = ufl.FiniteElement("Real",domain.ufl_cell(),0)
V = ufl.MixedElement([P1,R0])
ME = fem.FunctionSpace(domain,V)

#Define Equation. Boundary condition is captured within weak form
v,d = ufl.TestFunctions(ME)
u,c = ufl.TrialFunctions(ME)

x = ufl.SpatialCoordinate(domain)
f = 10 * ufl.exp(-((x[0] - 0.5)**2 + (x[1] - 0.5)**2) / 0.02)
g = -ufl.sin(5 * x[0])

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx + ufl.inner(c,v)*ufl.dx + ufl.inner(u,d)*ufl.dx
L = ufl.inner(f, v) * ufl.dx + ufl.inner(g, v) * ufl.ds 

#Solve with direct solver
problem = fem.petsc.LinearProblem(a, L, petsc_options={"ksp_type": "preonly",
                                                        "pc_type": "lu"
                                                        })
uh = problem.solve()

#Output to xdmf
xdmf = io.XDMFFile(domain.comm,"poisson5.xdmf","w")
#Output mesh
xdmf.write_mesh(domain)
#Output solution
u_,c_ = uh.split()
xdmf.write_function(u_)
xdmf.write_function(c_)


