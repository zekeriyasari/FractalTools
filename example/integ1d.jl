# This file includes an example file for 1D integrate. 

using FractalTools 
using GLMakie 

# Generate data 
xi, dx, xf = 0., 0.01, 10.
x = collect(xi : dx : xf) 
y = sin.(x) 
pts = collect.(zip(x, y))
res = integrate(pts, Interp1D(0.0001))
