import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import dx, grad, inner
import ufl

#Define time parameters
t = 0.0
tend = 100
nsteps = 100
dt = (tend-t)/nsteps

#Define mesh
nx = ny = 100
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[-50, -50],[50, 50]], [nx,ny])

# domain = mesh.create_unit_square(MPI.COMM_WORLD,nx,ny,mesh.CellType.triangle)
P1 = ufl.FiniteElement("Lagrange",domain.ufl_cell(),2)
V = fem.FunctionSpace(domain,P1)
v = ufl.TestFunction(V)

def initial_condition(x,r0=1.0*5):
    return 0.5*(1-np.tanh((np.sqrt((x[0])**2 +(x[1])**2)-r0)/(np.sqrt(2))))

#Create output file
xdmf = io.XDMFFile(domain.comm,"nb_100.xdmf","w")
xdmf.write_mesh(domain)

#Define solution variable and interpolate initial solution
#Solution variable from previous timestep
c0 = fem.Function(V)
c0.name="c0"
c0.interpolate(initial_condition)

c = fem.Function(V)
c.name="$c$" 
c.interpolate(initial_condition)
xdmf.write_function(c,t)

#Compute g' and p'
g = (c**2)*(1-c)**2
dgdc = ufl.diff(g,c)

DeltaF = np.sqrt(2)/30.0
p = (c**3)*(10 - 15*c + 6*c**2)
dpdc = ufl.diff(p,c)

#Define problem
F = (inner(c,v)*dx - inner(c0,v)*dx + 
    dt*inner(grad(c),grad(v))*dx + 
    dt*inner(dgdc,v)*dx -
    DeltaF*dt*inner(dpdc,v)*dx
    )

#Setup solver
problem = NonlinearProblem(F, c)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "residual"
solver.atol = 1e-8
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "gmres"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()

#Timestepping

#Introduce skipping for output to output only every nth step
stride = 2
counter = 0

while t < tend:
    #Update t
    t += dt
    #Solve
    res = solver.solve(c)
    print(f"Step {int(t/dt)}: num iterations: {res[0]}")
    #Update c0
    c0.x.array[:] = c.x.array
    #Write to xdmf
    counter = counter +1
    if counter % stride == 0:
        xdmf.write_function(c,t)

xdmf.close()