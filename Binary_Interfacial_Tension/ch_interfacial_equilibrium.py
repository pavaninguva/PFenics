import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, SpatialCoordinate
import ufl
import matplotlib.pyplot as plt
#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

def model(nx, kappa, profile=True):

    #Create mesh
    domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(-4.0,4.0), nx=nx)

    #Specify basis functions
    P1 = ufl.FiniteElement("CG",domain.ufl_cell(),1)
    ME = fem.FunctionSpace(domain,P1*P1)
    q,v = ufl.TestFunction(ME)

    V0, dofs = ME.sub(0).collapse()
    v1, dofs1 = ME.sub(1).collapse()
    cells, types, x = plot.create_vtk_mesh(V0)

    #Define solution variables and split
    u = fem.Function(ME)
    c,mu = ufl.split(u)

    #Define Equation
    f = 0.25*(1-c**2)**2
    def f(c):
        return 0.25*(1-c**2)**2
    
    def dfdc_(c):
        return -c*(1-c**2)
    
    dfdc = ufl.diff(f(c),ufl.variable(c))

    F0 = inner(grad(mu),grad(q))*dx
    F1 = inner(mu,v)*dx - inner(dfdc,v)*dx - kappa*inner(grad(c),grad(v))*dx

    F = F0 + F1

    #Apply IC
    def ic_c(x):
        return np.tanh(x[0]/np.sqrt(kappa))
    
    # def ic_mu(x):
    #     val = dfdc_(ic_c(x)) + np.tanh(x[0]/np.sqrt(2*kappa))*(1/(np.cosh(x[0]/np.sqrt(2*kappa))**2))
    #     return val
    
    u.x.array[:] = 0.0
    u.sub(0).interpolate(ic_c)
    # u.sub(1).interpolate(ic_mu)
    u.x.scatter_forward()

    problem = NonlinearProblem(F, u)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.convergence_criterion = "residual"
    solver.atol = 1e-8

    ksp = solver.krylov_solver
    opts = PETSc.Options()
    option_prefix = ksp.getOptionsPrefix()
    opts[f"{option_prefix}ksp_type"] = "preonly"
    opts[f"{option_prefix}pc_type"] = "lu"
    ksp.setFromOptions()

    sol = solver.solve(u)

    if profile == True:
        return x[:,0], u.x.array[dofs].real
    else:
        dc = np.gradient(u.x.array[dofs].real, x[:,0])
        dcsq = dc**2
        sigma = np.trapz(kappa*dcsq, x=x[:,0])
        return sigma



#Plot profiles
x1, c1, = model(nx=500,kappa=0.01)

fig2, ax2 = plt.subplots(1,1,num=2)
ax2.plot(x1,c1, label=r"$\kappa = 0.001$")
ax2.legend()
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$c$")
# ax2.set_ylim(bottom=-1.2, top=1.2)
fig2.tight_layout()

plt.show()
