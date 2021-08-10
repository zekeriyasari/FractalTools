using FractalTools 
using GLMakie 
using Triangulate

# Dataset 
f(x, y) = x^2 + y^2 
Ω = uniformdomain(3, BigFloat) 
ipts = interiorpoints(f, Ω, 100)
bpts = boundarypoints(f, Ω, 10) 
pts = [bpts; ipts]
dataset = Dataset(pts, bpts)
# Construct interpolant 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluations
tipts = interiorpoints(Ω, 100)
tbpts = boundarypoints(Ω, 10)
ivals = map(pnt ->  interp(pnt...), tipts)
ivals = map(pnt ->  interp(pnt...), tbpts)

# Test location 
tpts = project(dataset.points)
tpt = tpts[1] 
tess = tessellate(dataset, method)
idx = locate(tess, tpt)

msh = tomesh(dataset)
msh2 = project(msh)
fig, ax, plt = wireframe(msh2)
scatter!(msh2.position, color=:black)
scatter!(collect(msh2[idx].points))
scatter!([tpt[1]], [tpt[2]])


