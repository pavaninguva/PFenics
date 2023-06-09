import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import ufl
import matplotlib.pyplot as plt
#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)


def model(nx, kappa, profile=True):

    domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(-4.0,4.0), nx=nx)

    #Specify basis functions
    P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
    V = fem.FunctionSpace(domain,P1)
    v = ufl.TestFunction(V)

    #Define problem
    c = fem.Function(V)
    c.name="$c$" 

    #Compute f'
    f = 0.25*(c**2 -1)**2
    dfdc = ufl.diff(f,c)

    F = kappa*inner(grad(c),grad(v))*dx + inner(dfdc,v)*dx

    def ic(x):
        return np.tanh(x[0]/kappa)

    #Apply IC
    c.interpolate(ic)

    #Setup solver

    problem = NonlinearProblem(F, c)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.convergence_criterion = "residual"
    solver.atol = 1e-8

    ksp = solver.krylov_solver
    opts = PETSc.Options()
    option_prefix = ksp.getOptionsPrefix()
    opts[f"{option_prefix}ksp_type"] = "preonly"
    opts[f"{option_prefix}pc_type"] = "lu"
    ksp.setFromOptions()

    sol = solver.solve(c)
    cells, types, x = plot.create_vtk_mesh(V)

    if profile == True:
        return x[:,0], c.x.array.real
    else:
        dc = np.gradient(c.x.array.real, x[:,0])
        dcsq = dc**2
        sigma = np.trapz(kappa*dcsq, x=x[:,0])
        return sigma
    

#Run plots
sigma_25 = []
sigma_50 = []
sigma_100 = []
sigma_150 = []
sigma_200 = []
sigma_250 = []

kappa_vals = np.linspace(0.001, 0.5, 50)

for kappa in kappa_vals:
    sigma_25.append(model(nx=25,kappa=kappa,profile=False))
    sigma_50.append(model(nx=50,kappa=kappa,profile=False))
    sigma_100.append(model(nx=100,kappa=kappa,profile=False))
    sigma_150.append(model(nx=150,kappa=kappa,profile=False))
    sigma_200.append(model(nx=200,kappa=kappa,profile=False))
    sigma_250.append(model(nx=250,kappa=kappa,profile=False))
    print("beep")

#Plot 
fig1, ax1 = plt.subplots(1,1,num=1)
ax1.loglog(kappa_vals,sigma_25,label=r"$nx = 25$")
ax1.plot(kappa_vals,sigma_50, label=r"$nx = 50$")
ax1.plot(kappa_vals,sigma_100, label=r"$nx = 100$")
ax1.plot(kappa_vals,sigma_150, label=r"$nx = 150$")
ax1.plot(kappa_vals,sigma_200, label=r"$nx = 200$")
ax1.plot(kappa_vals,sigma_250, label=r"$nx = 250$")
ax1.plot(kappa_vals, (2/3)*(2**0.5)*np.sqrt(kappa_vals),"rx",label="Analytical")
ax1.legend()
ax1.set_ylabel(r"$\sigma$")
ax1.set_xlabel(r"$\kappa$")
fig1.tight_layout()
plt.savefig("ac_ift_equilibrium_vals.png",dpi=300)

#Plot profiles
x1, c1 = model(nx=250,kappa=0.001)
x2, c2 = model(nx=250,kappa=0.005)
x3, c3 = model(nx=250,kappa=0.01)
x4, c4 = model(nx=250,kappa=0.05)
x5, c5 = model(nx=250, kappa=0.1)
x6, c6 = model(nx=250, kappa=0.5)

fig2, ax2 = plt.subplots(1,1,num=2)
ax2.plot(x1,c1, label=r"$\kappa = 0.001$")
ax2.plot(x2,c2, label=r"$\kappa = 0.005$")
ax2.plot(x3,c3, label=r"$\kappa = 0.01$")
ax2.plot(x4,c4, label=r"$\kappa = 0.05$")
ax2.plot(x5,c5, label=r"$\kappa = 0.1$")
ax2.plot(x6,c6, label=r"$\kappa = 0.5$")
ax2.legend()
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$c$")
fig2.tight_layout()
plt.savefig("ac_profile_equilibrium_vals.png",dpi=300)


plt.show()
