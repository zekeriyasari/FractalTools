# This file includes mapping of interpolation interval to different subintervals. 

using FractalTools
using Plots 

# Construct interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
d = 0.01 * ones(length(x) - 1)

# Interpolate the data 
interp = interpolate(x, y, d) 

# Get transformations 
ws = interp.ifs.ws

# Get an example transformation 
pts = hcat.(eachcol([x'; y']))

# Plots 
plt = plot(layout=(2,1))
for (i, w) in enumerate(ws)
    v = getindex.(w.(pts), 1)
    plot!(v, ones(length(v)) * i, label="", subplot=1)
end
for (i, w) in enumerate(ws)
    v = getindex.(w.(pts), 2)
    plot!(ones(length(v)) * (i - 1), v, label="", subplot=2)
end
display(plt)
