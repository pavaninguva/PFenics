import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, ln, atan, tan, dot
import ufl

#Define mesh
nx = ny = 50
Lx = Ly = 5.0
domain = mesh.create_rectangle(MPI.COMM_WORLD,[[0.0,0.0],[Lx,Ly]], [nx,ny])

#Create output file
xdmf = io.XDMFFile(domain.comm,"Anisotropic.xdmf","w")
xdmf.write_mesh(domain)

P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
ME = fem.FunctionSpace(domain,P1*P1)
q,v = ufl.TestFunctions(ME)

#Define solution variables and split 
u = fem.Function(ME)
u0 = fem.Function(ME)

DT,c = ufl.split(u)
DT0,c0 = ufl.split(u0)

#model parameters
D_T = 2.25
N = 6.0
k3 = 0.02
k1 = 0.9
k2 = 20.0
tau_c = 3e-4

#Timestepping
dt = 0.001
t = 0.0
tend = 0.1

#Heat balance
F0 = inner(DT,q)*dx - inner(DT0, q)*dx + D_T*dt*inner(grad(DT),grad(q))*dx - (inner(c,q) - inner(c0,q))*dx

#Phase Field
#var.dx(i) is dc/x_i :O
psi = np.pi/8 + atan(c.dx(1)/c.dx(0))
phi = tan((N/2)*psi)
beta = (1-phi**2)/(1 + phi**2)
dbetadpsi = ufl.diff(beta,ufl.variable(psi))
D_tensor = ufl.as_tensor([[1 + k3*beta , -k3*dbetadpsi], [k3*dbetadpsi,1+k3*beta]])

m = c - 0.5 - (k1/np.pi)*atan(k2*DT)
F1 = tau_c*(inner(c,v) - inner(c0,v))*dx + dt*inner(dot(D_tensor,grad(c)),grad(v))*dx - dt*inner(c*(1-c)*m, v)*dx

#Couple
F = F0 + F1

#Define ICs
def ic_c(x):
    radius = 0.5
    C = (Lx/2, Ly/2)
    
    x_ = x[0]
    y_ = x[1]

    s_vals = np.sqrt((x_-C[0])**2 + (y_-C[1])**2)
    f_vals = np.zeros_like(s_vals)
    for idx, val in np.ndenumerate(s_vals):
        if val <= radius:
            f_vals[idx] = 1.0

    return f_vals

def ic_T(x):
    return -0.5*np.ones_like(x[0])

# Zero u
u.x.array[:] = 0.0
#Apply initial condition to c
u.sub(1).interpolate(ic_c)
#Apply IC to DT
u.sub(0).interpolate(ic_T)
u.x.scatter_forward()
c = u.sub(1)
DT = u.sub(0)
#Write IC to u0
u0.x.array[:] = u.x.array
#Output IC to file
xdmf.write_function(c,t)
xdmf.write_function(DT,t)


#Setup solver
problem = NonlinearProblem(F, u)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.atol = 1e-8
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly"
opts[f"{option_prefix}pc_type"] = "lu"
ksp.setFromOptions()


#Timestepping
while t < tend - 1e-8:
    #Update t
    t += dt
    #Solve
    res = solver.solve(u)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    #Update u0
    u0.x.array[:] = u.x.array
    xdmf.write_function(c,t)
    xdmf.writa_function(DT,t)
xdmf.close()


