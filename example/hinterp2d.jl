# This file includes and example file for HInterp2D interpolation 

using FractalTools 
using GeometryBasics
using Makie 

# Generate data 
f(x, y) = [
    x^2 + y^2 + 1, 
    x^2 - y^2
    ]
ngon = Triangle(
    Point(BigFloat(0.), BigFloat(0.)), 
    Point(BigFloat(1.), BigFloat(0.)), 
    Point(BigFloat(0.5), BigFloat(1.)))
npts = 100
pts = getdata(f, ngon, npts)
interp = interpolate(pts, HInterp2D(0.01 * ones(2,2)))

function err(x, y)
    fval = f(x, y) 
    ival = interp(x, y) 
    map(item -> abs.(item), fval - ival) ./ map(item -> abs.(item), fval) * 100
end

tpts = getdata(ngon, npts)

fig, ax, plt = trisurf(tpts, (x, y) -> f(x, y)[1], meshcolor3=:red)
trisurf!(ax, tpts, (x, y) -> interp(x, y)[1], meshcolor3=:blue)
fig, ax, plt = trisurf(tpts, (x, y) -> err(x, y)[1], meshcolor3=:red)
