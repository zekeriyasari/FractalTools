using FractalTools 
import Interpolations
using Makie 

# Generate interpolation data 
f = FractalTools.wen
xi, dx, xf = 0., 0.01, 1. 
x = collect(xi : dx : xf) 
y = f.(x )
pts = collect.(zip(x, y))

# Contruct fractal interpolation 
freevar = 0.1
interpf = interpolate(pts, Interp1D(freevar))

# Construct spline interpolation 
interps = Interpolations.CubicSplineInterpolation(range(x[1], x[end], step=x[2] - x[1]), y)

# Test domain 
xt = collect(xi : 0.1dx : xf) 
fvals = f.(xt) 
ivals = interpf.(xt) 
svals = interps.(xt)
evalsf = abs.(fvals - ivals)
evalss = abs.(fvals - svals)

msef = sum(evalsf.^2) / length(evalsf)
mses = sum((fvals - svals).^2) / length(svals)

@show msef 
@show mses

# Plots 
fig = Figure(resolution=(1920, 1080)) 
ax1 = fig[1, 1] = Axis(fig) 
ax2 = fig[2, 1] = Axis(fig, 
    title="freevar = $freevar, MSE Fractal Interp=$(round(msef, digits=3)), " *  
    "MSE Spline Interp=$(round(mses, digits=3))") 
lines!(ax1, xt, fvals, color=:red, label="True Values")
lines!(ax1, xt, ivals, color=:blue, label="Fractal Interp")
lines!(ax1, xt, svals, color=:green, label="Spline Interp")
scatter!(ax1, x, y, markersize=5, color=:black, label="Interp data")
lines!(ax2, xt, evalsf, color=:blue, label="Fractal Interp Error")
lines!(ax2, xt, evalss, color=:red, label="Spline Interp Error")
lb1 = fig[1, 2] = Legend(fig, ax1, framevisible = false)
lb2 = fig[2, 2] = Legend(fig, ax2, framevisible = false)
save(joinpath(@__DIR__, "interp1d_comparison.png"), fig)
display(fig)
