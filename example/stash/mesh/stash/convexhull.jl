using PyCall 
using Makie 

# Construct data 
spt = pyimport_conda("scipy.spatial", "scipy") 

# Generate random points 
pnts = [Point(rand(2)...) for i in 1 : 100] 

# Find the convex hull 
hull = spt.ConvexHull(collect(hcat(collect.(pnts)...)') )
vtx =  [Point(val[1], val[2]) for val in eachrow(hull.points[hull.vertices .+ 1 , :])]

fig, ax, plt = scatter(pnts, color=:black, markersize=5)
lines!(vtx, linewidth=3)
lines!(vtx[[1, end]], linewidth=3)
fig 

