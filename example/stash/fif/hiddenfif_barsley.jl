# This file includes one dimensional interpolation. 

using FractalTools 
using Plots 

# Generate interpolation data 
x = [0, 30,60, 100]
y = [0, 50, 40, 10]
z = [0, 30, 60, 100]
d = [0.3, 0.3, 0.3] * 0.01 # d must be small enaugh.
h = [0.2, 0.2, 0.1]
l = [-0.1, -0.1, -0.1]
m = [0.3, 0., -0.1]

# Construct interpolant 
# ifs = gethiddenifs(x, y, z, d, h, l, m)
interp = hiddenfif(x, y, z, d, h, l, m, niter=10)

# Plot data 
xd = sort(unique([x; collect(range(x[1], x[end], length=length(x) * 500))]))
yd = interp.(xd)
plot(xd, getindex.(yd, 1))
scatter!(x, y, marker=(:circle, 2))

