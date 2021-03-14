
using FractalTools 
using Makie 
using AbstractPlotting
import AbstractPlotting.GeometryBasics: TriangleFace 

# Construct delaunay 
p1, p2, p3 = [-1., √3.], [-1., -√3], [2., 0.]
tridln = TriDelaunay(p1, p2, p3, addboundarypoints=true)

while npoints(tridln) ≤ 1000 
    addpoint!(tridln) 
end 

finegrain!(tridln, 100)

# Define function 
f(x, y) = x^2 + y^2

# Plot the result 
fig, ax, plt = tridelaunayplotf(tridln, f, vmarkersize=10)
wireframe!(tomesh(tridln.delaunay), linewidth=3)
display(fig)
