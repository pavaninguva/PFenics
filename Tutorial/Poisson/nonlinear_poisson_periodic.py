from mpi4py import MPI
from dolfinx import mesh, fem, io
import ufl
from petsc4py.PETSc import ScalarType
from dolfinx_mpc import LinearProblem, MultiPointConstraint
import numpy as np
from Nonlinear_MPC import NonlinearMPCProblem,NewtonSolverMPC

"""
Implementing the modified example from
https://fenicsproject.org/olddocs/dolfin/1.3.0/python/demo/documented/nonlinear-poisson/python/documentation.html
"""
