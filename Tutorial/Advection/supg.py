import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
import dolfinx
import ufl
import matplotlib.pyplot as plt

"""
https://www.firedrakeproject.org/demos/DG_advection.py.html

https://bitbucket.org/fenics-project/dolfin/src/master/python/demo/undocumented/advection-diffusion/demo_advection-diffusion.py
"""

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

#Create mesh
nx = 100
Lx = 1.0
domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(0.0,Lx), nx=nx)

#Timestepping
CFL = 1.0
dt = CFL*(Lx/nx)/(1.0)

#Define function space
V = fem.FunctionSpace(domain, ("DG", 1))

#No treatment of BCs are needed

#Define initial condition
def initial_condition(x):
    return (1/(np.sqrt(2*np.pi)))*np.exp(-((x[0]-0.5)**2)/0.00025)

u0 = fem.Function(V)
u0.interpolate(initial_condition)

#Store u0
u0_vals = u0.x.array.real

#Define solution variable 
uh = fem.Function(V)
uh.interpolate(initial_condition)

#Define model equations
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Velocity vector
vel = fem.Constant(domain, PETSc.ScalarType((1.0,)))

#Midpoint solution
u_mid = 0.5*(u0 + u)


#Define equation
F = v*(u - u0)*ufl.dx + dt*v*ufl.dot(vel, ufl.grad(u_mid))*ufl.dx

#SUPG Stabilization Term
vnorm = ufl.sqrt(ufl.dot(vel,vel))
tdim = domain.topology.dim
numcells = domain.topology.index_map(tdim).size_local
h = dolfinx.cpp.mesh.h(domain,tdim,range(numcells))

res = u - u0 + dt*ufl.dot(vel, ufl.grad(u_mid))
stab = (h[0]/(2.0*vnorm))*ufl.dot(vel, ufl.grad(v))*res*ufl.dx

F += stab

a = ufl.lhs(F)
L = ufl.rhs(F)

#Setup solver
problem = fem.petsc.LinearProblem(a,L,
                                  petsc_options={"ksp_type": "preonly", "pc_type": "lu"})




#Solve
Tend = 1*dt
t = 0.0

while t < Tend - 1e-8:
    t += dt
    uh = problem.solve()
    u0.x.array[:] = uh.x.array

uh_vals = uh.x.array

#Plot
cells, types, x = plot.create_vtk_mesh(V)

#Plot IC
def ic(x):
    return (1/(np.sqrt(2*np.pi)))*np.exp(-((x-0.5)**2)/0.00025)

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(x[:,0],ic(x[:,0]),label="Initial Conditions")
ax1.plot(x[:,0],uh_vals,"--k",label="u(t=%s)"%t)
ax1.set_ylabel(r"$u$")
ax1.set_xlabel(r"$x$")
ax1.legend()
fig1.tight_layout()

plt.show()

