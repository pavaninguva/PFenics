import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner, ln
import ufl
import matplotlib.pyplot as plt
import pandas as pd

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

def model(nx, dt, tend, Nchi):
    #Parameters
    kappa = (2/3)*Nchi

    t = 0.0
    n_steps = int((tend - t)/dt)
    n_output = 50
    stride = int(n_steps/n_output)


    #Create mesh
    domain = mesh.create_interval(comm=MPI.COMM_WORLD, points=(0.0,4.0), nx=nx)

    #Specify basis functions
    P1 = ufl.FiniteElement("CG",domain.ufl_cell(),1)
    ME = fem.FunctionSpace(domain,P1*P1)
    q,v = ufl.TestFunction(ME)

    V0, dofs = ME.sub(0).collapse()
    v1, dofs1 = ME.sub(1).collapse()
    cells, types, x = plot.create_vtk_mesh(V0)

    #Define solution variables and split 
    u = fem.Function(ME)
    u0 = fem.Function(ME)

    c,mu = ufl.split(u)
    c0,mu0 = ufl.split(u0)

    c = ufl.variable(c)
    f = c*ln(c) + (1-c)*ln(1-c) + Nchi*c*(1-c)
    dfdc = ufl.diff(f,c)

    F0 = inner(c,q)*dx - inner(c0,q)*dx + (c*(1-c))*dt*inner(grad(mu),grad(q))*dx
    F1 = inner(mu,v)*dx - inner(dfdc,v)*dx - kappa*inner(grad(c),grad(v))*dx
    F = F0 + F1

    def initial_condition(x):
        values = 0.5 + 0.02*(0.5-np.random.rand(x.shape[1]))
        return values
    
    #Apply IC
    u.x.array[:] = 0.0
    u.sub(0).interpolate(initial_condition)
    u.x.scatter_forward()
    c = u.sub(0)
    u0.x.array[:] = u.x.array

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

    #Set output
    counter = 0

    output_df = pd.DataFrame()
    #Append xvals 
    output_df["x"] = x[:,0]
    output_df["0.0"] = u.x.array[dofs].real
    # local_vol = fem.assemble_scalar(fem.form(fem.Constant(domain, 1)*dx))
    # vol = MPI.reduce(local_vol, MPI.SUM, root=0)
    # local_u = fem.form(u*dx)
    while t < tend - 1e-8:
        try:
            #Update t
            t += dt
            #Solve
            res = solver.solve(u)
            print(f"Step {int(t/dt)}: num iterations: {res[0]}")
            #Update u0
            # l_c = fem.assemble_scalar(local_u)
            # lc_u = MPI.reduce(l_c, MPI.SUM, root=0)
            # print(lc_u/vol)
            u0.x.array[:] = u.x.array
        
            counter += 1
            if counter % stride == 0:
                output_df["%s"%t] = u.x.array[dofs].real
        except:
            #Output last converged value for
            output_df["%s"%t] = u0.x.array[dofs].real
            break
    return output_df
    
#Plot and final
df1 = model(50,0.000005,2.0, 30.0)
print(df1)

counter = 0
for time in list(df1.columns)[1:]:
    fig1,ax1 = plt.subplots(1,1)
    ax1.plot(df1["x"],df1[time])
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$c$")
    ax1.set_ybound(lower=0.0, upper=1.0)
    ax1.set_title("t = %s"%time)
    fig1.tight_layout()
    plt.savefig("plot_%s.png"%counter,dpi=300)
    plt.clf()
    plt.close()
    counter += 1


