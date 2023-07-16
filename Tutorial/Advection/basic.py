import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
import ufl
import matplotlib.pyplot as plt

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
V = fem.FunctionSpace(domain, ("CG", 1))

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


a = ufl.inner(u,v)*ufl.dx 
L = ufl.inner(u0,v)*ufl.dx - dt*ufl.dot(vel,ufl.grad(u0))*v*ufl.dx

#Setup solver
problem = fem.petsc.LinearProblem(a,L,
                                  petsc_options={"ksp_type": "preonly", "pc_type": "lu"})




#Solve
print(dt)
Tend = 0.01
t = 0.0

while t < Tend - 1e-8:
    t += dt
    uh = problem.solve()
    u0.x.array[:] = uh.x.array

print(t)
uh_vals = uh.x.array

#Plot
cells, types, x = plot.create_vtk_mesh(V)

#Plot IC
def ic(x):
    return (1/(np.sqrt(2*np.pi)))*np.exp(-((x-0.5)**2)/0.00025)

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(x[:,0],ic(x[:,0]),label="Initial Conditions")
ax1.plot(x[:,0],uh_vals,"--k",label=r"$u(t=0.01)$")
ax1.set_ylabel(r"$u$")
ax1.set_xlabel(r"$x$")
ax1.legend()
fig1.tight_layout()
plt.savefig("naive.png",dpi=300)

plt.show()