# This file includes mapping of interpolation interval to different subintervals. 

using FractalTools
using Plots 

# Construct interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
z = cos.(2π * x)
d = 0.9 * ones(length(x) - 1)
h = 0.9 * ones(length(x) - 1)
l = 0.9 * ones(length(x) - 1)
m = 0.9 * ones(length(x) - 1)

# Interpolate the data 
@time hiddeninterp = hiddenfif(x, y, z, d, h, l, m) 

# Get transformations 
ws = hiddeninterp.ifs.ws

# Get an example transformation 
pts = hcat.(eachcol([x'; y'; z']))

# Plots 
plt = plot(layout=(4,1))
plot!(x, sin.(2π*x), subplot=1, label="")
vline!(x, subplot=1, label="", ls=:dot)
for (i, w) in enumerate(ws)
    v = getindex.(w.(pts), 1)
    plot!(v, ones(length(v)) * i, label="", subplot=2)
end
for (i, w) in enumerate(ws)
    v = getindex.(w.(pts), 2)
    plot!(ones(length(v)) * i, v, label="", subplot=3, lw=3)
end
for (i, w) in enumerate(ws)
    v = getindex.(w.(pts), 3)
    plot!(ones(length(v)) * i, v, label="", subplot=4)
end
display(plt)

# # Checks 
# ifs = hiddeninterp.ifs

# # Analytical 
# ei = x[1 : end - 1]
# ai = x[2 : end] - ei 
# fi = y[1 : end - 1] - h 
# ci = y[2 : end] - h - fi
# gi = z[1 : end - 1] - h 
# ki = z[2 : end] - h - gi 

# # Computed 
# e = [w.b[1] for w in ws] 
# a = [w.A[1, 1] for w in ws] 
# f = [w.b[2] for w in ws] 
# c = [w.A[2, 1] for w in ws] 
# g = [w.b[3] for w in ws] 
# k = [w.A[3, 1] for w in ws] 

# # Compare 
# e .== ei 
# a .== ai
# f .== fi
# c .== ci
# g .== gi
# k .≈ ki
