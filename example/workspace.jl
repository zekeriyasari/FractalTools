using FractalTools 
using GLMakie 

# Dataset 
f(x, y) = x^2 + y^2 
Ω = uniformdomain(3) 
ipts = interiorpoints(f, Ω, 100)
bpts = boundarypoints(f, Ω, 10) 
pts = [bpts; ipts]
dataset = Dataset(pts, bpts)

# msh = tomesh(dataset)

# _msh = project(msh) 
# mesh(msh, color=last.(msh.position)) 
# wireframe!(msh, color=:black) 
# scatter!(msh.position)

# Construct interpolant 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluations
tipts = interiorpoints(Ω, 100)
tbpts = boundarypoints(Ω, 10)
ivals = map(pnt ->  interp(pnt...), tipts)
ivals = map(pnt ->  interp(pnt...), tbpts)

