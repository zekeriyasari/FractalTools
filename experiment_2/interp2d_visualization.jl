# This file plots two dimensional interpolation with irregular data 

using FractalTools 
using GeometryBasics
using Makie 

# Construct interpolation data 
f = FractalTools.ackley
α = 10
ngon = [BigFloat.([-1., -1.]), BigFloat.([1., -1]), BigFloat.([0., 1.])] * α
npts = 100
pts = getdata(f, ngon, npts)

# Construct interpolant 
freevar = 0.001
interp = interpolate(pts, Interp2D(freevar))

# Define error function 
erf(x, y) = abs(f(x, y) - interp(x, y))

# Construct test data 
tpts = getdata(ngon, 7 * npts)

# Plots 
fig = Figure()
ls1 = fig[1, 1] = LScene(fig)
ls2 = fig[1, 2] = LScene(fig)
trisurf!(ls1, tpts, f, meshcolor3=:red)
trisurf!(ls1, tpts, interp, meshcolor3=:blue)
scatter!(ls1, tpts)
trisurf!(ls2, tpts, erf, meshcolor3=:green)
