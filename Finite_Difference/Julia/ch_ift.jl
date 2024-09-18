using Pkg
Pkg.activate("CahnHilliard")
using DifferentialEquations
using Plots
using Sundials
using LaTeXStrings
using ForwardDiff
# using Dierckx
using BSplineKit
using NonlinearSolve

function spline_generator(χ,N1,N2,knots)

    #Def log terms 
    log_terms(ϕ) =  (ϕ./N1).*log.(ϕ) .+ ((1 .-ϕ)./N2).*log.(1 .-ϕ)

    function tanh_sinh_spacing(n, β)
        points = 0.5 * (1 .+ tanh.(β .* (2 * collect(0:n-1) / (n-1) .- 1)))
        return points
    end
    
    phi_vals_ = collect(tanh_sinh_spacing(knots-2,14))
    f_vals_ = log_terms(phi_vals_)

    #Append boundary values
    phi_vals = pushfirst!(phi_vals_,0)
    f_vals = pushfirst!(f_vals_,0)
    push!(phi_vals,1)
    push!(f_vals,0)

    
    # spline = Spline1D(phi_vals,f_vals)
    # d_spline(phi) = Dierckx.derivative(spline,phi)
    spline = BSplineKit.interpolate(phi_vals, f_vals,BSplineOrder(4))
    d_spline = Derivative(1)*spline

    df_spline(phi) = d_spline.(phi) .+ χ.*(1 .- 2*phi)

    return df_spline
end

# Function to compute interfacial tension
function compute_interfacial_tension(ϕ, xvals, κ)
    # Compute the derivative dϕ/dx using central differences
    dϕdx = similar(ϕ)
    dϕdx[1] = (ϕ[2] - ϕ[1]) / (xvals[2] - xvals[1])
    dϕdx[end] = (ϕ[end] - ϕ[end-1]) / (xvals[end] - xvals[end-1])
    # dϕdx[2:end-1] = (ϕ[3:end] - ϕ[1:end-2]) / (2 * dx)
    for i in 2:(length(xvals)-1)
        h1 = xvals[i] - xvals[i-1]
        h2 = xvals[i+1] - xvals[i]
        dϕdx[i] = (ϕ[i+1] - ϕ[i-1]) / (h1 + h2)
    end

    # Compute the integrand κ * (dϕ/dx)^2
    integrand = κ .* (dϕdx.^2)

    # Perform trapezoidal integration
    σ = trapezoidal_integration(xvals, integrand)

    return ForwardDiff.value(σ)
end

# Function for trapezoidal integration
function trapezoidal_integration(x, y)
    n = length(x)
    integral = 0.0
    for i in 1:(n-1)
        integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
    end
    return integral
end

# Define the CH equations with spatial discretization
function CH(ϕ, dx, params)
    χ, κ, N₁, N₂, energy_method = params
    spline = spline_generator(χ,N₁,N₂,100)

    dfdphi = ϕ -> begin 
        if energy_method == "analytical"
            -2 .* χ .* ϕ .+ χ - (1/N₂).*log.(1-ϕ) .+ (1/N₁).*log.(ϕ)
        else
            spline.(ϕ)
        end
    end
    mobility(ϕ) = ϕ .* (1 .- ϕ)

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

    # Define chemical potential
    μ = similar(ϕ)
    μ[1] = dfdphi(ϕ[1]) - (2 * κ / (dx^2)) * (ϕ[2] - ϕ[1])
    μ[end] = dfdphi(ϕ[end]) - (2 * κ / (dx^2)) * (ϕ[end-1] - ϕ[end])
    μ[2:end-1] = dfdphi.(ϕ[2:end-1]) - (κ / (dx^2)) .* (ϕ[3:end] - 2 .* ϕ[2:end-1] .+ ϕ[1:end-2])

    # Define LHS (time derivative of ϕ)
    f = similar(ϕ)
    f[1] = (2 / (dx^2)) * (M_func_half(ϕ[1], ϕ[2]) * (μ[2] - μ[1]))
    f[end] = (2 / (dx^2)) * (M_func_half(ϕ[end], ϕ[end-1]) * (μ[end-1] - μ[end]))
    f[2:end-1] = (1 / (dx^2)) .* (M_func_half.(ϕ[2:end-1], ϕ[3:end]) .* (μ[3:end] .- μ[2:end-1]) .-
                                   M_func_half.(ϕ[2:end-1], ϕ[1:end-2]) .* (μ[2:end-1] .- μ[1:end-2]))

    return f
end

# Function to solve the system and return values at t = 20
function solve_CH_at_time_20(L, nx, params, t_end)
    dx = L / (nx - 1)
    xvals = LinRange(0, L, nx)

    # Define the sigmoid function for initial conditions
    function sigmoid(x, a, b)
        return 1 ./ (1 .+ exp.(-a .* (x - b)))
    end

    # Parameters for the sigmoid
    a = 5        # Controls the slope of the transition
    b = L / 2    # Midpoint of the transition

    # Define the initial condition across the domain
    ϕ₀ = 0.2 .+ 0.6 .* sigmoid.(xvals, a, b)

    # Define ODE system
    function ode_system!(du, u, p, t)
        du .= CH(u, dx, params)
    end

    # Set up the problem
    prob = ODEProblem(ode_system!, ϕ₀, (0.0, t_end))

    # Solve the problem
    # sol = solve(prob, Rosenbrock23(), reltol=1e-6, abstol=1e-8,nlsolve_kwargs=RobustMultiNewton())
    sol = solve(prob, Rosenbrock23(), reltol=1e-6, abstol=1e-8)
    # Get the solution at t = t_end
    ϕ_at_t_end = sol(t_end)

    # Compute the interfacial tension at t = t_end
    σ_at_t_end = compute_interfacial_tension(ϕ_at_t_end, xvals, κ)

    return ϕ_at_t_end, σ_at_t_end
end

# Parameters
L = 5.0    #Length of the domain
nx = 401   #Number of spatial points
χ = 0.3  # Parameter χ
κ = (2/3) * χ  # Parameter κ
N₁ = 100
N₂ = 50
t_end = 1000
energy_method="spline"

params = (χ, κ, N₁, N₂, energy_method)

# Call the function to solve and get values at t = 20
ϕ_at_20, σ_at_20 = solve_CH_at_time_20(L, nx, params,t_end)
println(ϕ_at_20)


# Print the results
println("Interfacial tension at t = 20: ", σ_at_20)
# # println("Values of ϕ at t = 20: ", ϕ_at_20)

#Generate IFT values and plot interfacial profile
# chi_vals_analytical1 = range(3,25,10)
# ift_vals_analytical1 = zeros(size(chi_vals_analytical1))

# p1 = plot()

# for i in 1:length(chi_vals_analytical1)
#     chi = chi_vals_analytical1[i]
#     ϕ_vals, σ = solve_CH_at_time_20(5.0,401,(chi,(2/3)*chi, 1,1,"spline"))
#     println(σ)
#     ift_vals_analytical1[i] = σ
#     plot!(p1,range(0,10,401),ϕ_vals, label=L"\chi = %$chi")
# end
# plot(p1,size=(500,500), grid=false, xlabel=L"x", ylabel=L"\phi",
#     tickfont=Plots.font("Computer Modern", 10),
#     legendfont=Plots.font("Computer Modern",8),dpi=300)

# println(ift_vals_analytical1)

