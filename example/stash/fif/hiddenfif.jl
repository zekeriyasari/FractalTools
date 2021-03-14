# This file includes one dimensional interpolation. 

using FractalTools 
using Plots 
plotly()

# Generate interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x .+ π/6)
z = cos.(2π * x)
u = 0.2
d = u * ones(length(x) - 1)
h = u * ones(length(x) - 1)
l = u * ones(length(x) - 1)
m = u * ones(length(x) - 1)

# Construct interpolant 
# ifs = gethiddenifs(x, y, z, d, h, l, m)
interp = hiddenfif(x, y, z, d, h, l, m, niter=40)

# Plot data 
xd = collect(xi : 0.01dx : xf)
yd = interp.(xd)

plt = plot(layout=(2,1))
plot!(xd, getindex.(yd, 1), label="original", subplot=1)
scatter!(x, y, marker=(:circle, 3), subplot=1)
plot!(xd, getindex.(yd, 2), label="hidden", subplot=2)
scatter!(x, z, marker=(:circle, 3), subplot=2)

