import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import ufl
import multiprocessing

#Define time parameters

def cahn_hilliard(r0,stride=8):

    t = 0.0
    tend = 40
    nsteps = 800 
    dt = (tend-t)/nsteps

    #Model parameters
    kappa = 1e-2
    M = 1.0

    #Define mesh
    nx = ny = 100
    domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0.0,0.0],[2*np.pi, 2*np.pi]], [nx,ny])

    #Create output file
    filename = str(r0) +".xdmf"
    xdmf = io.XDMFFile(domain.comm,filename,"w")
    xdmf.write_mesh(domain)

    #Define model
    P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),2)
    ME = fem.FunctionSpace(domain,P1*P1)
    q,v = ufl.TestFunctions(ME)

    #Define solution variables and split 
    u = fem.Function(ME)
    u0 = fem.Function(ME)

    c,mu = ufl.split(u)
    c0,mu0 = ufl.split(u0)

    #Define Equation
    c = ufl.variable(c)
    f = 0.25*(1-c**2)**2
    dfdc = ufl.diff(f,c)


    F0 = inner(c,q)*dx - inner(c0,q)*dx + M*dt*inner(grad(mu),grad(q))*dx
    F1 = inner(mu,v)*dx - inner(dfdc,v)*dx - kappa*inner(grad(c),grad(v))*dx
    F = F0 + F1

    x_list = [x*np.pi for x in [0.5,1.5]]
    y_list = [x*np.pi for x in [1.0, 1.0]]
    r_list = [x*np.pi for x in [r0,0.25]]
    loc_list = [list(x) for x in zip(x_list,y_list,r_list)]


    def initial_condition(x, kappa = kappa, loclist = loc_list):
        #Unpack x and y 
        x_ = x[0]
        y_ = x[1]

        def f(x_,y_,loc, kappa):
            xval,yval,rval = loc
            s_vals = np.sqrt((x_-xval)**2 + (y_-yval)**2) -rval
            f_vals = np.zeros_like(s_vals)
            for idx,s in np.ndenumerate(s_vals):
                if s < 0:
                    f_vals[idx] = 2*np.exp(-kappa/(s**2))
                else:
                    f_vals[idx] = 0.0
            return f_vals
        
        sum_ = 0.0
        for loc_ in loclist:
            sum_ = sum_ + f(x_,y_,loc_,kappa)

        return -1.0 + sum_

    # Zero u
    u.x.array[:] = 0.0
    #Apply initial condition
    u.sub(0).interpolate(initial_condition)
    u.x.scatter_forward()
    c = u.sub(0)
    #Write IC to u0
    u0.x.array[:] = u.x.array
    #Output IC to file
    xdmf.write_function(c,t)

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
    stride = stride
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
            xdmf.write_function(c,t)
    xdmf.close()

#Run cases
rvals = np.linspace(0.2,0.3,11).tolist()
it = iter(rvals)
for val in it:
    val1 = round(val,2)
    val2 = round(next(it,np.nan),2)
    if np.isnan(val2) == True:
        process1 = multiprocessing.Process(target=cahn_hilliard, args=(val1,))
        process1.start()
        process1.join()
    else:
        process1 = multiprocessing.Process(target=cahn_hilliard, args=(val1,))
        process2 = multiprocessing.Process(target=cahn_hilliard,args=(val2,))
        process1.start()
        process2.start()
        print("waiting...")
        process1.join()
        process2.join()
