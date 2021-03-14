# This file includes an example file for HInterp1D interpolation. 

using FractalTools 
using Makie 

# Generate data 
xi, dx, xf = 0., 0.1, 10.
x = collect(xi : dx : xf) 
y = sin.(x) 
z = cos.(x)
pts = collect.(zip(x, y, z))
freevars = [0.01 * rand() * ones(2,2) for i in 1 : length(x) - 1]
interp = interpolate(pts, HInterp1D(freevars))

# Calculate interpolant 
xt = collect(xi : 0.1dx : xf)
_yt = interp.(xt)
yt = getindex.(_yt, 1)

# Plot interpolation 
fig, ax, plt = lines(xt, yt, color=yt, linewidth=5)
scatter!(ax, x, y, markersize=5, color=:black)
fig 
