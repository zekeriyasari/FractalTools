# Import Meshes.atol to package 
# Meshes.atol(::Type{BigFloat}) = 100 * eps()

using FractalTools 
using GLMakie 
using Triangulate

# Dataset 
f(x, y) = x^2 + y^2 
Ω = [
    [BigFloat("-1.0"), BigFloat("-1.0")],
    [BigFloat("1.0"), BigFloat("-1.0")],
    [BigFloat("1.0"), BigFloat("1.0")],
    [BigFloat("0.0"), BigFloat("0.5")],
    [BigFloat("-1.0"), BigFloat("1.0")],
]
ipts = interiorpoints(f, Ω, 100)
bpts = boundarypoints(f, Ω, 10) 
pts = [bpts; ipts]
dataset = Dataset(pts, bpts)

# Construct interpolant 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluations
tipts = interiorpoints(Ω, 100)
ivals = map(pnt ->  interp(pnt...), tipts)


# tbpts = boundarypoints(Ω, 10)
# ivals = map(pnt ->  interp(pnt...), tbpts)

# # Test location 
# tpts = project(dataset.points)
# tpt = tpts[1] 
# tess = tessellate(dataset, method)
# idx = locate(tess, tpt)

msh = tomesh(dataset)
msh2 = project(msh)

fig, ax, plt = mesh(msh, color=last.(msh.position))
wireframe!(msh, linewidth=3)
mesh!(msh2, color=last.(msh2.position))
wireframe!(msh2, color=:black)
scatter!(msh2.position, color=:black)
scatter!(collect(msh2[idx].points))
scatter!([tpt[1]], [tpt[2]])


fig, ax, plt = wireframe(msh2) 

pts = project(dataset.points)
idx = 10 
tpnt = pts[idx]
tval = interp(tpnt[1], tpnt[2])
scatter!(Point(tpnt...))

indexes = map(enumerate(pts)) do (i, pnt) 
    try
        tval = interp(pnt...)
        return true 
    catch 
        return false
    end
end 
n = length(pts)
nfound = length(pts[indexes])
nnotfound = length(pts[.!indexes])
