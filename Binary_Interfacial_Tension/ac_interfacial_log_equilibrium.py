import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, ln
import ufl
import matplotlib.pyplot as plt

from scipy.optimize import fsolve

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

def fh_minima(Nchi):
    f = lambda x: -2*Nchi*x + Nchi - np.log(1-x) + np.log(x)
    x_vals = fsolve(f, [0.001, 0.999])
    return x_vals



def model(nx, Nchi, profile=True):
    domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(-8.0,8.0), nx=nx)

    #Specify basis functions
    P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
    V = fem.FunctionSpace(domain,P1)
    v = ufl.TestFunction(V)

    #Define problem
    c = fem.Function(V)
    c.name="$c$" 

    #Define chemical potential
    f = c*ln(c) + (1-c)*ln(1-c) + Nchi*c*(1-c)
    dfdc = ufl.diff(f,c)

    kappa = (2/3)*Nchi
    F = kappa*inner(grad(c),grad(v))*dx + inner(dfdc,v)*dx

    def ic(x):
        xvals = fh_minima(Nchi)
        f = (xvals[0] - xvals[1])/(1+ np.exp(x[0]*Nchi)) + xvals[1]
        return f
    
    #Apply IC
    c.interpolate(ic)

    problem = NonlinearProblem(F, c)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.convergence_criterion = "residual"
    solver.atol = 1e-6

    ksp = solver.krylov_solver
    opts = PETSc.Options()
    option_prefix = ksp.getOptionsPrefix()
    opts[f"{option_prefix}ksp_type"] = "preonly"
    opts[f"{option_prefix}pc_type"] = "lu"
    opts[f"{option_prefix}ksp_max_it"] = "1000"
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


#Plot convergence
sigma_50 = []
sigma_100 = []
sigma_150 = []
sigma_200 = []
sigma_250 = []

Nchi_vals = np.linspace(2.5,7.0,25)

for Nchi in Nchi_vals:
    sigma_50.append(model(nx=50,Nchi=Nchi,profile=False))
    sigma_100.append(model(nx=100,Nchi=Nchi,profile=False))
    sigma_150.append(model(nx=150,Nchi=Nchi,profile=False))
    sigma_200.append(model(nx=200,Nchi=Nchi,profile=False))
    sigma_250.append(model(nx=250,Nchi=Nchi,profile=False))
    print("beep")

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(Nchi_vals,sigma_50,label=r"$nx = 50$")
ax1.plot(Nchi_vals,sigma_100,label=r"$nx = 100$")
ax1.plot(Nchi_vals,sigma_150,label=r"$nx = 150$")
ax1.plot(Nchi_vals,sigma_200,label=r"$nx = 200$")
ax1.plot(Nchi_vals,sigma_250,label=r"$nx = 250$")
ax1.legend()
ax1.set_ylabel(r"$\sigma$")
ax1.set_xlabel(r"$N\chi$")
fig1.tight_layout()
plt.savefig("ac_ift_equilibrium_log.png",dpi=300)

#Plot profiles
x1, c1 = model(nx=100,Nchi=2.5)
x2, c2 = model(nx=100,Nchi=3.5)
x3, c3 = model(nx=100,Nchi=4.5)
x4, c4 = model(nx=100,Nchi=5.5)
x5, c5 = model(nx=100,Nchi=6.5)

fig2, ax2 = plt.subplots(1,1,num=2)
ax2.plot(x1,c1, label=r"$N\chi = 2.5$")
ax2.plot(x2,c2, label=r"$N\chi = 3.5$")
ax2.plot(x3,c3, label=r"$N\chi = 4.5$")
ax2.plot(x4,c4, label=r"$N\chi = 5.5$")
ax2.plot(x5,c5, label=r"$N\chi = 6.5$")
ax2.legend()
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$c$")
fig2.tight_layout()
plt.savefig("ac_profile_equilibrium_log.png",dpi=300)

plt.show()