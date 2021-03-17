# This file investigates mean square error (MSE) with respect to number of interpolation points in 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x)      
xi = 0. 
xf = 1. 
freevar = 0.001                
npts = 10 : 2 : 50
ntpts = 10 * npts[end]

# Compute errors 
mse = map(npts) do npt 
    # Construct interpolation data 
    x = collect(range(xi, xf, length=npt))
    y = f.(x)
    pts = collect.(zip(x, y))

    # Construct interpolation 
    interp = interpolate(pts, Interp1D(freevar))

    # Construct test data 
    xt = collect(range(xi, xf, length=ntpts))

    # Compute error
    sum((f.(xt) - interp.(xt) ).^2) / length(xt)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Number of Points", ylabel="MSE", title="1D Interpolation MSE") 
stem!(ax, npts, mse, color=:black)
save(joinpath(@__DIR__, "interp1d_mse.png"), fig)
display(fig)

