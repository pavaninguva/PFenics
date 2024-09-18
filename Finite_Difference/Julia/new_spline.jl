using Plots
using LaTeXStrings
using Statistics
# using Dierckx
using BSplineKit


@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 100)
    n = length(y)
    sx, sy = x[1:step:n], y[1:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := []
        y := []
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end


#Define spline approximation

function spline_generator(chi,N1,N2,knots)

    #Def log terms 
    log_terms(phi) =  (phi./N1).*log.(phi) .+ ((1 .-phi)./N2).*log.(1 .-phi)

    function tanh_sinh_spacing(n, β)
        points = 0.5 * (1 .+ tanh.(β .* (2 * collect(0:n-1) / (n-1) .- 1)))
        return points
    end
    
    phi_vals_ = collect(tanh_sinh_spacing(knots-2,14))
    # phi_vals_ = collect(range(1e-16,1-1e-16, knots-2))

    f_vals_ = log_terms(phi_vals_)

    #Append boundary values
    phi_vals = pushfirst!(phi_vals_,0)
    f_vals = pushfirst!(f_vals_,0)

    push!(phi_vals,1)
    push!(f_vals,0)

    #This code uses Dierckx.jl

    # spline = Spline1D(phi_vals,f_vals)
    # d_spline(phi) = Dierckx.derivative(spline,phi)

    # f_spline(phi) = spline(phi) .+ chi.*phi.*(1 .- phi)
    # df_spline(phi) = d_spline.(phi) .+ chi.*(1 .- 2*phi)

    #This code uses BSplinekit.jl
    spline = BSplineKit.interpolate(phi_vals, f_vals,BSplineOrder(4))
    
    # Define derivative function
    d_spline = Derivative(1)*spline

    # Define the final spline and its derivative functions
    f_spline(phi) = spline(phi) .+ chi .* phi .* (1 .- phi)
    df_spline(phi) = d_spline(phi) .+ chi .* (1 .- 2 * phi)


    return f_spline, df_spline
end

#Define fh and fh_deriv

function fh(phi,chi,N1,N2)
    f = (1/N1).*phi.*log.(phi) .+ (1/N2).*(1 .- phi).*log(1 .- phi) .+ chi.*phi.*(1 .- phi)
    return f
end

function fh_deriv(phi,chi,N1,N2)
    df = (1/N1).*log.(phi) .+ (1/N1) .- (1/N2).*log.(1 .- phi) .- (1/N2) .- 2*chi.*phi .+ chi
    return df
end

phi_vals = collect(range(1e-16,1-1e-16,1000))

#Generate splines
f1, df1 = spline_generator(50,1,1,100)
f2, df2 = spline_generator(1,100,50,100)

p1 = plot(legend=:bottom, 
grid=false,xlabel=L"\phi", ylabel=L"f")
#Plot Analytical 
plot!(p1,phi_vals,fh.(phi_vals,50,1,1),label="",lw=2,linecolor="black",alpha=1.0,seriestype=:samplemarkers,m = (5, :white, stroke(1, :blue)))
plot!(p1,phi_vals,fh.(phi_vals,1,100,50),label="",lw=2,linecolor="black",alpha=1.0,seriestype=:samplemarkers,)

plot!(p1,phi_vals,f1.(phi_vals),label="",linecolor="red",linestyle=:dot,linewidth=2,alpha=0.5)
plot!(p1,phi_vals,f2.(phi_vals),label="",linecolor="red",linestyle=:dot,linewidth=2,alpha=0.5)

#Plot df
p2 = plot(legend=:bottomleft, size=(500, 500),
grid=false,xlabel=L"\phi", ylabel=L"\frac{\partial f}{\partial \phi}")
plot!(p2,phi_vals,fh_deriv.(phi_vals,50,1,1),label="Analytical, "*L"\chi=50, N_{1}=N_{2} = 1",lw=2,linecolor="black",alpha=1.0,seriestype=:samplemarkers,m = (5, :white, stroke(1, :blue)))
plot!(p2,phi_vals,fh_deriv.(phi_vals,1,100,50),label="Analytical, "*L"\chi=10, N_{1}= 100, N_{2} = 1",lw=2,linecolor="black",alpha=1.0,seriestype=:samplemarkers,)

plot!(p2,phi_vals,df1.(phi_vals),label="Spline, "*L"n=100",linecolor="red",linestyle=:dot,linewidth=2,alpha=0.5)
plot!(p2,phi_vals,df2.(phi_vals),label="",linecolor="red",linestyle=:dot,linewidth=2,alpha=0.5)


plot(p1,p2, layout= (1,2), size=(800, 400),
    tickfont=Plots.font("Computer Modern", 10),
    legendfont=Plots.font("Computer Modern",7),dpi=300,
    bottom_margin = 3Plots.mm, left_margin = 3Plots.mm, right_margin=3Plots.mm
    )
# # savefig("./spline_plot_bsplinekit.png")


### Plot RMSE
# knots_values = [10,20,25,40,50,75,100,150,200]
# rmse_f = zeros(length(knots_values))
# mae_f = zeros(length(knots_values))

# rmse_df = zeros(length(knots_values))
# mae_df = zeros(length(knots_values))

# rmse_df_int = zeros(length(knots_values))
# mae_df_int = zeros(length(knots_values))


# for j in 1:length(knots_values)
#     knots = knots_values[j]

#     phi_test_ = range(1e-16,1-1e-16,100)
#     phi_test_int = range(1e-2,1-1e-2,100)
    
#     f_test, df_test = spline_generator(50,1,1,knots)

#     f_ana_vals = fh.(phi_test_,50,1,1)
#     df_ana_vals = fh_deriv.(phi_test_,50,1,1)

#     df_ana_int_vals = fh_deriv.(phi_test_int,50,1,1)

#     #Compute errors
#     rmse_f[j] = sqrt(mean((f_test.(phi_test_).- f_ana_vals).^2))
#     mae_f[j] = maximum(abs.(f_test.(phi_test_).- f_ana_vals))

#     rmse_df[j] = sqrt(mean((df_test.(phi_test_).- df_ana_vals).^2))
#     mae_df[j] = maximum(abs.(df_test.(phi_test_).- df_ana_vals))

#     rmse_df_int[j] = sqrt(mean((df_test.(phi_test_int).- df_ana_int_vals).^2))
#     mae_df_int[j] = maximum(abs.(df_test.(phi_test_int).- df_ana_int_vals))

# end


# p3 = plot()
# plot!(p3,knots_values,rmse_f,label=L"\textrm{RMSE}, f",color="black")
# plot!(p3,knots_values,mae_f,label=L"\textrm{MAE}, f",linestyle=:dash,color="black")
# plot!(p3,knots_values,rmse_df,label=L"\textrm{RMSE}, \frac{\partial f}{\partial \phi}",color="red")
# plot!(p3,knots_values,mae_df,label=L"\textrm{MAE}, \frac{\partial f}{\partial \phi}",linestyle=:dash,color="red")
# plot!(p3,knots_values,rmse_df_int,label=L"\textrm{RMSE, Interior}, \frac{\partial f}{\partial \phi}",color="blue")
# plot!(p3,knots_values,mae_df_int,label=L"\textrm{MAE, Interior}, \frac{\partial f}{\partial \phi}",linestyle=:dash,color="blue")


# plot(p3,legend=:bottomleft, size=(500, 500),
#     grid=false,xlabel=L"N", ylabel=L"\textrm{Error}",
#     tickfont=Plots.font("Computer Modern", 10),
#     legendfont=Plots.font("Computer Modern",7),dpi=300,
#     xaxis=:log, yaxis=:log
#     )



# savefig("./spline_error_bspline.png")

