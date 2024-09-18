using Plots
using NonlinearSolve
using NaNMath
using LaTeXStrings
# using Dierckx
using BSplineKit


#Compute the binodal
function solve_binodal(chi, N1, N2)
    # Define the binodal function
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

    # Set up the problem and initial guess
    initial_guess = [1e-6, 1-1e-6]
    params = (chi, N1, N2)
    problem = NonlinearProblem(binodal!, initial_guess, params)

    # Solve using NonlinearSolve
    solution = solve(problem, RobustMultiNewton())
    eq_vals = sort!(solution.u)

    # Return the values of phi
    return eq_vals
end

function spline_generator(chi,N1,N2,knots)

    #Def log terms 
    log_terms(phi) =  (phi./N1).*log.(phi) .+ ((1 .-phi)./N2).*log.(1 .-phi)

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
    spline = BSplineKit.interpolate(phi_vals, f_vals,BSplineOrder(4))
    f_spline(phi) = spline(phi) .+ chi.*phi.*(1 .- phi)

    return f_spline
end

function fh(chi,N1,N2, u)
    f = ((1/N1).*u.*log.(u) + (1/N2).*(1 .- u).*log.(1 .- u) .+ chi.*u.*(1 .- u))
    return f
end

function trapezoidal_integration(x, y)
    n = length(x)
    integral = 0.0
    for i in 1:(n-1)
        # integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
        # Replace NaN values in y with 0
        y_i = isnan(y[i]) ? 0.0 : y[i]
        y_ip1 = isnan(y[i+1]) ? 0.0 : y[i+1]

        # Trapezoidal rule integration step
        integral += 0.5 * (x[i+1] - x[i]) * (y_ip1 + y_i)
    end
    return integral
end

function interfacial_tension(chi,N1,N2,energy_method)
    #Solve equilibrium composition
    eq_vals = solve_binodal(chi, N1, N2)

    #Define kappa
    asym_factor = N2/N1
    if asym_factor < 0.1
        kappa = 1/3*chi
    else 
        kappa = 2/3*chi
    end

    spline = spline_generator(chi,N1,N2,100)

    f = u -> begin
        if energy_method == "analytical"
            fh(chi,N1,N2,u)
        else
            spline.(u)
        end
    end
    df(u) = f(u) -(u.*(f(eq_vals[2]) - f(eq_vals[1])) .+ (eq_vals[2]*f(eq_vals[1]) - eq_vals[1]*f(eq_vals[2])))./(eq_vals[2]-eq_vals[1])
    int_fun(u) = NaNMath.sqrt.(2 .* kappa .*df.(u))

    #Compute integral
    phi_vals = range(eq_vals[1],eq_vals[2],1000)
    int_vals = int_fun(phi_vals)
    val = trapezoidal_integration(phi_vals,int_vals)

    return val
end


#Generate IFT Plot
chi_vals1 = range(2.5,50,length=100)
chi_vals2 = range(0.02,0.3,length=100)
ift_vals_analytical1 = zeros(length(chi_vals1))
ift_vals_spline1 = zeros(length(chi_vals1))

ift_vals_analytical2 = zeros(length(chi_vals2))
ift_vals_spline2 = zeros(length(chi_vals2))

for i in 1:length(chi_vals1)
    ift_vals_analytical1[i]= interfacial_tension(chi_vals1[i],1,1,"analytical")
    ift_vals_spline1[i]= interfacial_tension(chi_vals1[i],1,1,"spline")

    ift_vals_analytical2[i]= interfacial_tension(chi_vals2[i],100,50,"analytical")
    ift_vals_spline2[i]= interfacial_tension(chi_vals2[i],100,50,"spline")
end

#IFT Values
chi_vals_ch_analytical1 = [3,4,5,6,7,10,12,14,16,18,20]
ift_vals_ch_analytical1 = [0.40375096501515584,0.9118986512530513,1.4087228114012793,1.8933843980393206,
                        2.3690717601911735,3.7667183016974364,4.6866237321860105,5.60224180916003,6.515380617645344,
                        7.426946903193865,8.337509477042541
                        ]

chi_vals_ch_spline1 = [3,4,5,6,7,10,12,14,16,18,20,22,25,28,30,35,40,45,48,50]
ift_vals_ch_spline1 = [0.40360817639169844,0.9118919924767018,1.4087042450771559,
                        1.8933632836311558,2.3690538195461386,3.7667009552830026,
                        4.686608274271352, 5.6022210693245675, 6.515361668188908,
                        7.426936692958273, 8.337496043885011, 9.247320498811423,
                        10.61104162382206, 11.974068017663338, 12.882033175332348,
                        15.151998962305925, 17.42196978709527, 19.68964890786439,
                        21.044889821948775,21.953235773514287
                        ]


chi_vals_ch_analytical2 = [0.035, 0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16, 0.18, 0.2, 0.22]
ift_vals_ch_analytical2 = [0.0005574063142541528,0.008874457714711155,0.013865145690294121,
                        0.018787530532206723,0.023641994813200288,0.028437263790641334,0.033184579881446734,
                        0.042573320782148275,0.05186579560539525,0.061096754971589395, 0.07028714573044764,
                        0.07944981431785422, 0.08859332662219838
                        ]

chi_vals_ch_spline2 = [0.035, 0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3]
ift_vals_ch_spline2 = [0.0005518920998109167,0.008873582373262376,0.013865008488235622,
                        0.018787368937804652,0.02364176755391312,0.02843700154134357,
                        0.033184297141725536,0.042573036741737715,0.05186556692718115,
                        0.061096527321797725,0.07028691606961184,0.07944960389662528,
                        0.08859308503011626, 0.09772207776069804,0.10684121057425087,
                        0.1159520318482597, 0.1250551365685994
                            ]

chicrit(N1,N2) = sqrt((N1^2) + (N2^2))/(2*N1*N2)

p = plot()
plot!(p,chi_vals1./2,ift_vals_analytical1,label=L"\textrm{Analytical (Full)}, N_{1} = N_{2} = 1",color="black")
plot!(p,chi_vals1./2,ift_vals_spline1,label=L"\textrm{Analytical (Spline)}, N_{1} = N_{2} = 1",linestyle=:dash,color="red")
plot!(p,chi_vals_ch_analytical1./2,ift_vals_ch_analytical1,label=L"\textrm{Simulation (Full)}, N_{1} = N_{2} = 1",seriestype=:scatter,markershape=:diamond, color="black",markeralpha=0.5)
plot!(p,chi_vals_ch_spline1./2,ift_vals_ch_spline1,label=L"\textrm{Simulation (Spline)}, N_{1} = N_{2} = 1",seriestype=:scatter,markershape=:+,color="black",markersize=7)

plot!(p,chi_vals2./(chicrit(100,50)),ift_vals_analytical2,label=L"\textrm{Analytical}, N_{1} = 100, N_{2} = 50",color="blue")
plot!(p,chi_vals2./(chicrit(100,50)),ift_vals_spline2,label=L"\textrm{Spline}, N_{1} = 100, N_{2} = 50",linestyle=:dashdot,color="red")
plot!(p,chi_vals_ch_analytical2./chicrit(100,50),ift_vals_ch_analytical2,label=L"\textrm{Simulation (Full)}, N_{1} = 100, N_{2} = 50",seriestype=:scatter,markershape=:diamond, color="blue",markeralpha=0.5)
plot!(p,chi_vals_ch_spline2./chicrit(100,50),ift_vals_ch_spline2,label=L"\textrm{Simulation (Spline)}, N_{1} = 100, N_{2} = 50",seriestype=:scatter,markershape=:+,color="blue",markersize=7)
plot(p,legend=:topleft, size=(500, 500),
    grid=false,xlabel=L"\chi_{12}/\chi_{C}", ylabel=L"\textrm{Scaled \ Interfacial \ Tension} \ \tilde{\sigma}",
    tickfont=Plots.font("Computer Modern", 10),
    legendfont=Plots.font("Computer Modern",8),dpi=300
    )