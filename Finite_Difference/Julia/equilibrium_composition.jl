using LinearAlgebra
using Plots
using NonlinearSolve
using NaNMath
using LaTeXStrings

# Define the function for the binodal calculations
function binodal!(F, phi, params)
    chi, N1, N2 = params

    function fh_deriv(phi, chi, N1, N2)
        df = (1/N1)*NaNMath.log(phi) + (1/N1) - (1/N2)*NaNMath.log(1-phi) - (1/N2) - 2*chi*phi + chi
        return df
    end

    function osmotic(phi, chi, N1, N2)
        osmo = phi*((1/N1)-(1/N2)) - (1/N2)*NaNMath.log(1-phi) - chi*phi^2
        return osmo
    end

    phiA = phi[1]
    phiB = phi[2]
    dF = fh_deriv(phiA, chi, N1, N2) - fh_deriv(phiB, chi, N1, N2)
    dO = osmotic(phiA, chi, N1, N2) - osmotic(phiB, chi, N1, N2)
    F[1] = dF
    F[2] = dO
end

# Function to solve for the binodal points
function binodal_solver(N1, N2, chi_scale_vals)
    phic = sqrt(N2)/(sqrt(N1) + sqrt(N2))
    chic = (sqrt(N1) + sqrt(N2))^2 / (2*N1*N2)

    phi1bn = zeros(length(chi_scale_vals))
    phi2bn = zeros(length(chi_scale_vals))

    x0 = [1e-3, 1 - 1e-3]

    for i in 1:length(chi_scale_vals)
        chi_scaled = chi_scale_vals[i] * chic
        params = (chi_scaled, N1, N2)
        prob = NonlinearProblem(binodal!, x0, params)
        sol = solve(prob, RobustMultiNewton())
        phi1bn[i] = sol.u[1]
        phi2bn[i] = sol.u[2]
        println("beep")
        x0 = copy(sol.u)
    end

    return [phic, chic], phi1bn, phi2bn
end

# Test
chi_vals = range(10.0, 1.0, length=1000)
critvals, phi1_, phi2_ = binodal_solver(100, 50, chi_vals)

# Plot the results
p = plot(phi1_, chi_vals .^ -1, color="blue", label="Binodal", lw=1)
plot!(p, phi2_, chi_vals .^ -1, color="blue", lw=1, label="")
scatter!(p, [critvals[1]], [1.0], color="black", ms=3, label="Critical Point")
xlims!(p, 0, 1)
ylims!(p, 0.0, 1.1)
xlabel!(p, L"\phi")
ylabel!(p, L"\frac{\chi_c}{\chi}")
plot!(p, legend=:topright, size=(600, 600), 
    tickfont=Plots.font("Computer Modern", 12), grid=false,
    legendfont=Plots.font("Computer Modern",8),dpi=300)

# Display the plot
# display(p)
