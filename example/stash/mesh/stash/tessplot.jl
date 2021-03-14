using FractalTools 
using Makie 

# Construct a tesselation 
pts = [0., 0.], [1., 0], [0.5, 1.] 
tridln = TriDelaunay(pts...) 

# Plot the delaunay. 
tridelaunayplot(tridln)
