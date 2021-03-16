# This file investigates the effect of free variables on the interpolation performance. 

using FractalTools 
using Makie 

# Settings 
f(x) = sin(2π * x)
xi, dx, xf = 0., 0.1, 1.
freevar = 0.01

# Construt interpolation 
x = collect(xi : dx : xf) 
y = f.(x) 
interp = interpolate(collect.(zip(x, y)), Interp1D(freevar))

# Define error function 
function erf(x) 
    fval = f(x) 
    ival = interp(x) 
    abs(fval - ival) / abs(fval) * 100
end 

# Get test data 
xt = collect(xi : 0.1dx : xf) 
fval = f.(xt) 
ival = interp.(xt) 
