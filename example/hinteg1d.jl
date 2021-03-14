# This file includes an example file for 1D integrate. 

using FractalTools 
using Makie 

# Generate data 
xi, dx, xf = 0., 0.01, 10.
x = collect(xi : dx : xf) 
y = sin.(x) 
z = cos.(x) 
pts = collect.(zip(x, y, z))
res = integrate(pts, HInterp1D(0.001 * ones(2, 2)))
