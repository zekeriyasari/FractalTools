using Interpolations 
using Makie 
using GeometryBasics 
import GeometryBasics: Ngon 
using StaticArrays
using FractalTools 

# Construct interpolation data 
f(x, y) = x^2 + y^2
rect = Rect2D(BigFloat(-1.), BigFloat(-1.), BigFloat(2.), BigFloat(2.)) 
msh2 = GeometryBasics.mesh(Tesselation(rect, (10, 10)))
pnts3d = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in msh2.position]

# Construct interpolant 
fitp = FractalTools.interpolate(pnts3d)
x = getindex.(pnts3d, 1) |> unique
y = getindex.(pnts3d, 2) |> unique
z = getindex.(pnts3d, 3) |> val -> reshape(val, 10, 10) |> collect
iitp = Interpolations.interpolate((x, y), z, Gridded(Linear()))
