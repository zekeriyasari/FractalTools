# This file investigates mean square error (MSE) with respect to number of interpolation points in 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x)      
g(x) = cos(2π * x)
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
    z = g.(x)
    pts = collect.(zip(x, y, z))

    # Construct interpolation 
    interp = interpolate(pts, HInterp1D(fill(freevar, 2, 2)))

    # Construct test data 
    xt = collect(range(xi, xf, length=ntpts))

    # Compute error
    fvals = getindex.(f.(xt), 1)
    ivals = getindex.(interp.(xt), 1)
    sum((fvals - ivals).^2) / length(xt)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Number of Points", ylabel="MSE", title="1D Hidden Interpolation MSE") 
stem!(ax, npts, mse, color=:black)
save(joinpath(@__DIR__, "hinterp1d_mse.png"), fig)
display(fig)
