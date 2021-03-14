# Test file 

using FractalTools
using Plots
pyplot()
# default(:size, (1900, 1000))

# Construct interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
z = cos.(2π * x)
u = 0.9
d = u * ones(length(x) - 1)
h = u * ones(length(x) - 1)
l = u * ones(length(x) - 1)
m = u * ones(length(x) - 1)

plt = plot(layout=(5,1))
for niter in 1 : 5
    interp = interpolate(x, y, d, niter=1)
    # interp = hiddenfif(x, y, z, d, h, l, m)
    xd = collect(xi : 0.1dx : xf)
    plot!(xd, getindex.(interp.(xd), 1), subplot=niter)
    scatter!(x, y, subplot=niter)
end
display(plt)

