# This file includes an example to complare FIF and HiddenFIF

using FractalTools
using Plots 

# Generate interpolation data 
xi, dx, xf = 0, 0.01, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
z = cos.(2π * x)
z2 = tan.(2π * x)
u = 0.1
d = u * ones(length(x) - 1)
h = 0.1 * ones(length(x) - 1)
l = 0.1 * ones(length(x) - 1)
m = 0.1 * ones(length(x) - 1)

# Interplate data 
interp = interpolate(x, y, d)
hiddeninterp = hiddenfif(x, y, z, d, h, l, m)
hiddeninterp2 = hiddenfif(x, y, z2, d, h, l, m)

# Calculate interpolation data 
xd = collect(xi : 0.01dx : xf)
yd = interp.(xd)
hidden_yd = hiddeninterp.(xd)
hidden_yd2 = hiddeninterp2.(xd)

# Plots 
plt = plot(layout=(2,1))
plot!(xd, getindex.(yd, 1), label="interpolate", subplot=1, lw=5)
scatter!(x, y, marker=(:circle, 5), subplot=1)
plot!(xd, getindex.(hidden_yd, 1), label="hidden", subplot=1)
scatter!(x, y, marker=(:circle, 5), subplot=1)
plot!(xd, getindex.(hidden_yd2, 1), label="hidden2", subplot=1)
scatter!(x, y, marker=(:circle, 5), subplot=1)
plot!(xd, getindex.(hidden_yd, 2), label="hidden_extra", subplot=2)
scatter!(x, z, marker=(:circle, 5), subplot=2)
plot!(xd, getindex.(hidden_yd2, 2), label="hidden_extra", subplot=2)
scatter!(x, z2, marker=(:circle, 5), subplot=2)


