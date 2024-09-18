import numpy as np
import ufl
from dolfinx import fem, mesh, plot
from mpi4py import MPI
from petsc4py import PETSc
import matplotlib.pyplot as plt

# Create mesh and function space
nx = 50
domain = mesh.create_interval(MPI.COMM_WORLD, nx, [0.0, 1.0])
V = fem.FunctionSpace(domain, ("CG", 1))

# Define boundary condition
u_D = fem.Function(V)
u_D.interpolate(lambda x: np.sin(2 * np.pi * x[0]))

def boundary(x):
    return np.isclose(x[0], 0.0)

bc = fem.dirichletbc(u_D, fem.locate_dofs_geometrical(V, boundary))

# Define the advection velocity
a = 1.0

# Define trial and test functions
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Time-stepping parameters
dt = 0.01
T = 2.0
num_steps = int(T / dt)

# Define the stabilization parameter tau
h = ufl.CellDiameter(domain)
tau = h / (2 * a)

# Create functions for current and previous time step
u_n = fem.Function(V)
u_n.interpolate(lambda x: np.sin(2 * np.pi * x[0]))  # Initial condition

u = fem.Function(V)

# Define the SUPG-stabilized weak form
F = ((u - u_n) / dt + a * ufl.grad(u)[0]) * (v + tau * a * ufl.grad(v)[0]) * ufl.dx
a_form = ufl.lhs(F)
L_form = ufl.rhs(F)

# Assemble the matrix (bilinear form)
A = fem.petsc.assemble_matrix(fem.form(a_form), bcs=[bc])
A.assemble()

# Time-stepping loop
solver = PETSc.KSP().create(MPI.COMM_WORLD)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.GMRES)
solver.getPC().setType(PETSc.PC.Type.JACOBI)

for n in range(num_steps):
    t = (n + 1) * dt
    print(f"Time step {n+1}/{num_steps}, Time {t:.2f}")

    # Update the previous solution
    u_n.x.array[:] = u.x.array[:]

    # Assemble the right-hand side (linear form)
    b = fem.petsc.create_vector(fem.form(L_form))
    with b.localForm() as loc_b:
        loc_b.set(0)
    fem.petsc.assemble_vector(b, fem.form(L_form))
    fem.petsc.apply_lifting(b, [fem.form(a_form)], [bc])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, [bc])

    # Solve the linear system
    solver.solve(b, u.vector)
    u.x.scatter_forward()

    # Plot the solution
    if n % 10 == 0 or n == num_steps - 1:
        plot.plot(u)
        plt.title(f"Time: {t:.2f}")
        plt.show()

print("Simulation completed.")