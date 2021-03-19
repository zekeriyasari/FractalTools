# This file investigates integration error with respect to number of points in 1D hidden integration 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x) + 1
g(x) = cos(2π * x) + 1 
xi = 0. 
xf = 1. 
freevar = 0.001                
npts = 10 : 2 : 50
ntpts = 10 * npts[end]

# Compute theoretical value 
fval = (cos(2π * xi) - cos(2π * xf)) / (2π) + (xf - xi)

# Compute errors 
mse = map(npts) do npt 
    # Construct interpolation data 
    x = collect(range(xi, xf, length=npt))
    y = f.(x)
    z = g.(x) 
    pts = collect.(zip(x, y, z))

    # Compute integration
    ival = integrate(pts, HInterp1D(fill(freevar, 2, 2)))

    # Compute error
    abs(fval - ival)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Number of Points", ylabel="Error", title="1D Hidden Integration Error") 
stem!(ax, npts, mse, color=:black)
save(joinpath(@__DIR__, "hinteg1d_error.png"), fig)
display(fig)
