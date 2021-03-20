using FractalTools 
using Makie 

# Generate interpolation data 
fs = [FractalTools.sinusoid, FractalTools.wen]
mses = map(fs) do f 
    xi, dx, xf = 0., 0.01, 1. 
    x = collect(xi : dx : xf) 
    y = f.(x )
    pts = collect.(zip(x, y))
    xt = collect(xi : 0.1dx : xf) 
    fvals = f.(xt) 

    # Contruct fractal interpolation f
    freevars = collect(0.01 : 0.01 : 0.9) 
    mse = map(freevars) do freevar 
        interp = interpolate(pts, Interp1D(freevar))
        ivals = interp.(xt) 
        evals = abs.(fvals - ivals)
        sum(evals.^2) / length(evals)
    end 
end 

# Plots 
fig = Figure() 
ax = [Axis(fig[i, 1]) for i in 1 : 2]
for i in 1 : 2 
    axx = ax[i]
    stem!(axx, freevars, mses[i], color=:black) 
    axx.xlabel = "Free Variable"
    axx.ylabel = "MSE"
    axx.title = "$(fs[i]) Interp1D MSE"
end 

save(joinpath(@__DIR__, "$(f)_interp1d_mse.png"), fig)
display(fig)
