using Pkg
Pkg.activate("CahnHilliard")
using NonlinearSolve

f(u, p) = u .* u .- p
u0 = [1.0, 1.0]
p = 2.0
prob = NonlinearProblem(f, u0, p)
sol = solve(prob)