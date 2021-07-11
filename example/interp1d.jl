# This file includes an example file for 1D interpolation. 

using FractalTools 
using GLMakie 

# Generate data 
xi, dx, xf = 0., 0.1, 10.
x = collect(xi : dx : xf) 
y = sin.(x) 
pts = collect.(zip(x, y))
interp = interpolate(pts, Interp1D(0.1))

# Calculate interpolant 
xt = collect(xi : 0.1dx : xf)
yt = interp.(xt)

# Plot interpolation 
fig, ax, plt = lines(xt, yt, color=yt, linewidth=5)
scatter!(ax, x, y, markersize=5, color=:black)
