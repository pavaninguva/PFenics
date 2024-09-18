using Pkg
Pkg.activate("CahnHilliard")
using DifferentialEquations
using Plots
using Sundials
using LaTeXStrings




# Parameters
L = 4.0    # Length of the domain
nx = 201    # Number of spatial points
dx = L / (nx - 1)
xvals = LinRange(0, L, nx)
tf = 20.0    # Final time

χ = 25
κ = (2/3)*χ


dfdphi(ϕ) = -2 .* χ.*ϕ .+ χ - log.(1-ϕ) + log.(ϕ)

mobility(ϕ) = ϕ*(1-ϕ)

function M_func_half(ϕ₁,ϕ₂,option=1)
    if option == 1
        M_func = 0.5 .*(mobility.(ϕ₁) .+ mobility.(ϕ₂))
    elseif option == 2
        M_func = mobility.(0.5 .* (ϕ₁ .+ ϕ₂))
    elseif option == 3
        M_func = (2 .* mobility.(ϕ₁) .* mobility.(ϕ₂)) ./ (mobility.(ϕ₁) .+ mobility.(ϕ₂))
    end
    return M_func
end


# Define CH equations with spatial discretization
function CH(ϕ,dx=dx,κ=κ)
    #Define chem pot 
    μ = similar(ϕ)
    μ[1] = dfdphi(ϕ[1]) - (2*κ/(dx^2))*(ϕ[2]-ϕ[1])
    μ[end] = dfdphi(ϕ[end]) - (2*κ/(dx^2))*(ϕ[end-1]-ϕ[end])
    μ[2:end-1] = dfdphi.(ϕ[2:end-1]) - (κ/(dx^2)).*(ϕ[3:end] -2 .*ϕ[2:end-1] .+ ϕ[1:end-2])

    #Define LHS
    f = similar(ϕ)

    f[1] = (2/(dx^2))*(M_func_half(ϕ[1],ϕ[2])*(μ[2] - μ[1]))
    f[end] = (2/(dx^2))*(M_func_half(ϕ[end],ϕ[end-1]))*(μ[end-1]-μ[end])
    f[2:end-1] = (1/(dx^2)) .*(M_func_half.(ϕ[2:end-1],ϕ[3:end]).*(μ[3:end].-μ[2:end-1]) .-
                                M_func_half.(ϕ[2:end-1],ϕ[1:end-2]).*(μ[2:end-1] .- μ[1:end-2])
                                ) 
    
    return f
end

#Wrap for ODE solver
function ode_system!(du,u,p,t)
    du .= CH(u,dx,κ)
end

# Define initial condition
ϕ₀ = 0.5 .+ 0.05 .* randn(nx)

# Set up the problem
prob = ODEProblem(ode_system!, ϕ₀, (0.0, tf))

# Solve the problem
sol = solve(prob, ROS3P(), reltol=1e-6, abstol=1e-8)

# Plot the results
anim = @animate for t in sol.t
    plot(xvals, sol(t), ylim=(-0.2, 1.2), xlabel="x", ylabel="u", title="Time = $(round(t, digits=2))")
end

gif(anim, "nonlinear_heat_equation.gif", fps=15)