using Pkg
Pkg.activate("CahnHilliard")

using Symbolics
using SparseArrays
using Printf
using Plots

# Parameters
nx = 100
L = 10.0
dx = L / (nx - 1)
xvals = LinRange(0, L, nx)

# Create symbolic variables
@variables t
@variables u[1:nx]

# Define the nonlinear heat equation symbolically
function k(u)
    1 + u^2
end

eqs = Vector{Equation}(undef, nx)

# Interior points
for i in 2:nx-1
    eqs[i] = Differential(t)(u[i]) ~ (k(u[i+1]) * (u[i+1] - u[i]) - k(u[i-1]) * (u[i] - u[i-1])) / dx^2
end

# Boundary conditions
eqs[1] = u[1] ~ 0.0
eqs[nx] = u[nx] ~ 0.0

# Substitute variables with symbolic placeholders to allow Symbolics.jacobian to work
@variables usym[1:nx]

# eqs_sym = [substitute(e, Dict(u .=> usym)) for e in eqs]

# Compute the Jacobian of the system
J = Symbolics.jacobian(eqs, usym)

# Generate the sparsity pattern
sparsity_pattern = Symbolics.sparsity_pattern(J)

# Convert to sparse matrix for visualization
sparse_matrix = SparseMatrixCSC(sparsity_pattern)

# Plot the sparsity pattern
plot()
spy(sparse_matrix, markersize=5, title="Jacobian Sparsity Pattern", grid=false)
savefig("jacobian_sparsity_pattern.png")