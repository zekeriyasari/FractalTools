# This file includes an example file to plot the interpolation errors.

using FractalTools 
using Makie 
using StaticArrays 
using GeometryBasics
import GeometryBasics: Ngon 
using Random 
Random.seed!(0)

# Construct ngon
ngon = Ngon(SVector([Point(BigFloat(cos(θ)), BigFloat(sin(θ))) for θ in 0 : π / 2 : 2π - π / 2]...))
f(x, y) = x^2 + y^2
pnts3d = getdata(f, ngon, 100)

# Construct interpolant 
itp = interpolate(pnts3d)

# Calculate interpolation error 
pnts2d1 = project(pnts3d)
fvals1 = map(pnt -> Point(pnt..., f(pnt...)), pnts2d1)
ivals1 = map(pnt -> Point(pnt..., itp(pnt...)), pnts2d1)
evals1 = map(pnt -> Point(pnt..., relerf(pnt...)), pnts2d1)

@show maximum(getindex.(evals1, 3))

pnts3d2 = getdata(f, ngon, 400)
pnts2d2 = project(pnts3d2) 
fvals2 = map(pnt -> Point(pnt..., f(pnt...)), pnts2d2)
ivals2 = map(pnt -> Point(pnt..., itp(pnt...)), pnts2d2)
evals2 = map(pnt -> Point(pnt..., relerf(pnt...)), pnts2d2)

@show maximum(getindex.(evals2, 3))

# Plots 
fig = Figure() 
ax11 = LScene(fig) 
ax12 = LScene(fig)
ax13 = LScene(fig) 
ax21 = LScene(fig) 
ax22 = LScene(fig) 
ax23 = LScene(fig) 
fig[1, 1] = vgrid!(hgrid!(ax11, ax12, ax13),hgrid!(ax21, ax22, ax23)) 


msh21 = project(triangulate(pnts3d)[2])
scatter!(ax11, pnts3d, color=:black) 
wireframe!(ax11, msh21)
trisurf!(ax11, fvals1, meshcolor3=last.(fvals1)) 
trisurf!(ax12, ivals1, meshcolor3=last.(ivals1)) 
wireframe!(ax12, msh21)
trisurf!(ax13, evals1, meshcolor3=last.(evals1)) 
wireframe!(ax13, msh21)

msh22 = project(triangulate(pnts3d2)[2])
scatter!(ax21, pnts3d2, color=:black) 
trisurf!(ax21, fvals2, meshcolor3=last.(fvals2)) 
wireframe!(ax21, msh22)
trisurf!(ax22, ivals2, meshcolor3=last.(ivals2)) 
wireframe!(ax22, msh22)
trisurf!(ax23, evals2, meshcolor3=last.(evals2)) 
wireframe!(ax23, msh22)
display(fig) 

