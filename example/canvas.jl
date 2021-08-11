using FractalTools
using GLMakie
using Triangulate 

# Data generation 
f(x, y) = x^2 + y^2
pts = [
    [-2., 0.], 
    [-1., 1.], 
    [0., 0.], 
    [1., 1.], 
    [2., 0.], 
    [2., -1.], 
    [1., 0.], 
    [0., -1.], 
    [-1., 0.], 
    [-2., -1.], 
]
bpts = boundarypoints(f, pts, 0.1)
tpts = interiorpoints(f, pts, 0.1)
dataset = Dataset(tpts, bpts)


# Tessellation 
points = tpts 
domain = bpts 
triin = Triangulate.TriangulateIO()
triin.pointlist = hcat(points...)
list = 1 : length(domain)
triin.segmentmarkerlist = collect(list) 
triin.segmentlist = hcat(collect.(zip(list, circshift(list, -1)))...)
args = verbose ? "p" : "pQ"
triout, vorout = Triangulate.triangulate(args, triin)
tess = tri.Triangulation(triout.pointlist[1, :], triout.pointlist[2, :], triout.trianglelist' .- 1)

# Plots 
push!(bpts, bpts[1]) 
fig, ax, plt = mesh(msh, color=last.(msh.position))
scatterlines(getindex.(pts, 1), getindex.(pts, 2), markersize=10)
catterlines!(getindex.(bpts, 1), getindex.(bpts, 2), markersize=2)
scatter!(getindex.(tpts, 1), getindex.(tpts, 2), markersize=2)
fig 

