# This file includes one dimensional interpolation. 

using FractalTools 
using Plots 
default(:size, (600,500))
theme(:default)

# Generate interpolation data
f(x) = sin(2π*x) 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = f.(x)
d = 0.1 * ones(length(x) - 1)

# Construct interpolant 
interp = interpolate(x, y, d, niter=10)

# Plot data 
xd = collect(xi : 0.1dx : xf)
yd = interp.(xd)
plt = plot(xd, yd)
plot!(xd, f.(xd))
scatter!(x, y, marker=(:circle, 2))
plt 
