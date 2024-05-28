import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, ln
import ufl

#Define time parameters
t = 0.0
tend = 400
dt = 1.0
theta = 1

#Parameters
Nchi_AB = 3
Nchi_AC = 3
Nchi_BC = 3

kappa_AA = (2.0 / 3.0) * Nchi_AC
kappa_BB = (2.0 / 3.0) * Nchi_BC
kappa_AB = (1.0 / 3.0) * (Nchi_AC +  Nchi_BC - 1.0 * Nchi_AB)

A0 = 0.2
B0 = 0.2

#Define mesh
nx = ny = 100
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0.0,0.0],[100.0, 100.0]], [nx,ny])

#Create output file
xdmf1 = io.XDMFFile(domain.comm,"a.xdmf","w")
xdmf1.write_mesh(domain)
xdmf2 = io.XDMFFile(domain.comm,"b.xdmf","w")
xdmf2.write_mesh(domain)
xdmf3 = io.XDMFFile(domain.comm, "c.xdmf","w")
xdmf3.write_mesh(domain)


#Define test functions
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
ME = ufl.MixedElement([P1,P1,P1,P1,P1])
VP = fem.FunctionSpace(domain,ME)
h_1, h_2, j_1, j_2, j_3 = ufl.TestFunctions(VP)

#Define solution variables and split 
u = fem.Function(VP)
u0 = fem.Function(VP)

a, b, mu_AB, mu_AC, mu_BC = ufl.split(u)
a0, b0, mu_AB0, mu_AC0, mu_BC0 = ufl.split(u0)

#Formulate equations
a = ufl.variable(a)
b = ufl.variable(b)
c = 1 - a - b
c = ufl.variable(c)
f = a*ln(a) + b*ln(b) + (1-a-b)*ln(1-a-b) + a*b*Nchi_AB + a*c*Nchi_AC + b*c*Nchi_BC

dfda = ufl.diff(f,a)
dfdb = ufl.diff(f,b)
dfdc = ufl.diff(f,c)

mu_AB_mid = (1-theta)*mu_AB0 + theta*mu_AB
mu_AC_mid = (1-theta)*mu_AC0 + theta*mu_AC
mu_BC_mid = (1-theta)*mu_BC0 + theta*mu_BC

F_a = (
    inner(a,h_1)*dx 
    - inner(a0,h_1)*dx
    + dt*a*b*inner(grad(mu_AB_mid),grad(h_1))*dx
    + dt*a*c*inner(grad(mu_AC_mid),grad(h_1))*dx
)

F_b = (
    inner(b,h_2)*dx
    - inner(b0,h_2)*dx
    - dt*a*b*inner(grad(mu_AB_mid),grad(h_2))*dx
    + dt*b*c*inner(grad(mu_BC_mid),grad(h_2))*dx
)

F_mu_AB = (
    inner(mu_AB,j_1)*dx
    -inner(dfda,j_1)*dx
    +inner(dfdb,j_1)*dx
    -(kappa_AA-kappa_AB)*inner(grad(a),grad(j_1))*dx
    +(kappa_BB-kappa_AB)*inner(grad(b),grad(j_1))*dx
)

F_mu_AC = (
    inner(mu_AC,j_2)*dx
    - inner(dfda,j_2)*dx
    + inner(dfdc,j_2)*dx
    - kappa_AA*inner(grad(a),grad(j_2))*dx
    - kappa_AB*inner(grad(b),grad(j_2))*dx
)

F_mu_BC = (
    inner(mu_BC,j_3)*dx
    -inner(dfdb,j_3)*dx
    +inner(dfdc,j_3)*dx
    -kappa_BB*inner(grad(a),grad(j_3))*dx
    -kappa_AB*inner(grad(b),grad(j_3))*dx
)

#Couple everything
F = F_a + F_b + F_mu_AB + F_mu_AC + F_mu_BC

#Apply initial conditions
def initial_conditions_a(x):
    values = A0 + 0.02*(0.5-np.random.rand(x.shape[1]))
    return values

def initial_conditions_b(x):
    values = B0 + 0.02*(0.5-np.random.rand(x.shape[1]))
    return values

u.x.array[:] = 0.0
u.sub(0).interpolate(initial_conditions_a)
u.sub(1).interpolate(initial_conditions_b)
u.x.scatter_forward()
u0.x.array[:] = u.x.array

#Write a0 and b0
a = u.sub(0)
b = u.sub(1)
xdmf1.write_function(a,t)
xdmf2.write_function(b,t)

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

#Introduce skipping for output to output only every nth step
stride = 10
counter = 0

#Timestepping
while t < tend:
    #Update t
    t += dt
    #Solve
    res = solver.solve(u)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    #Update u0
    u0.x.array[:] = u.x.array
    counter = counter +1
    if counter % stride == 0:
        xdmf1.write_function(a,t)
        xdmf2.write_function(b,t)
xdmf1.close()
xdmf2.close()