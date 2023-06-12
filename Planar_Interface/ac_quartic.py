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

def model(nx, dt, tend, kappa, profile=True, animation=False):

    #Timestepping parameters
    t = 0.0
    n_steps = int((tend - t)/dt)
    n_output = 25
    stride = int(n_steps/n_output)

    #Generate mesh
    domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(-6.0,6.0), nx=nx)

    #Specify basis functions
    P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),1)
    V = fem.FunctionSpace(domain,P1)
    v = ufl.TestFunction(V)

    #Initial conditions
    def ic(x, kappa = kappa):
        return 0.5*(1-np.tanh((x[0])/(2*np.sqrt(2*kappa))))
    
    #Define solution variable and interpolate initial solution
    #Solution variable from previous timestep
    c0 = fem.Function(V)
    c0.name="c0"
    c0.interpolate(ic)

    c = fem.Function(V)
    c.name="$c$" 
    c.interpolate(ic)
    
    #Compute f'
    f = 0.25*(1- c**2)**2
    dfdc = ufl.diff(f,c)

    F = inner(c,v)*dx - inner(c0,v)*dx + kappa*dt*inner(grad(c),grad(v))*dx + dt*inner(dfdc,v)*dx
    
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

    while t < tend - 1e-8:
        t += dt
        #Solve
        res = solver.solve(c)
        print(f"Step {int(t/dt)}: num iterations: {res[0]}")
        #Update c0
        c0.x.array[:] = c.x.array[:]
    cells, types, x = plot.create_vtk_mesh(V)
    
    if profile == True:
        return x[:,0], c.x.array.real
    

#Plot profile
def analytical (x,t,kappa):
    u = (3*np.sqrt(kappa))/(np.sqrt(2))
    return 0.5*(1-np.tanh((x-u*t)/(2*np.sqrt(2*kappa))))


x1, c1 = model(nx=100, dt=0.01, tend=5.0, kappa=0.05)
x2, c2 = model(nx=200, dt = 0.05, tend=5.0, kappa = 0.05)
x3, c3 = model(nx=500, dt = 0.05, tend=5.0, kappa=0.05)
ic_vals = analytical(x3,0.0,0.05)
ana_vals = analytical(x3,5.0,0.05)


fig1, ax1 = plt.subplots(1,1,num=1)
ax1.plot(x3,ic_vals,"--k",label=r"Initial Conditions")
ax1.plot(x3,ana_vals,"--r" ,label=r"Analytical")
ax1.plot(x1,c1, label=r"nx = 100, dt = 0.01")
ax1.plot(x2,c2, label=r"nx = 200, dt = 0.05")
ax1.plot(x3,c3, label=r"nx = 500, dt = 0.05")
ax1.set_ylabel(r"$c$")
ax1.set_xlabel(r"$x$")
ax1.legend()
fig1.tight_layout()
plt.savefig("ac_travelingwave_profile.png",dpi=300)

#Plot RMSE
# nlist = [50,100,150,200,250,300,350,400,450,500]
dt_list = [0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.002, 0.001]
rmse_vals = []
for n in dt_list:
    xval, cval = model(nx = 200, dt=n, tend=5.0, kappa=0.05)
    ana_val = analytical(xval, 5.0,kappa =0.05)
    rmse_vals.append(np.sqrt(np.square(np.subtract(cval,ana_val)).mean()))

fig2, ax2 = plt.subplots(1,1,num=2)
ax2.loglog(dt_list,rmse_vals)
ax2.set_xlabel(r"$dt$")
ax2.set_ylabel(r"RMSE")
fig2.tight_layout()
plt.savefig("ac_travelingwave_error.png",dpi=300)


plt.show()

        

