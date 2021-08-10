# Import Meshes.atol to package 
# Meshes.atol(::Type{BigFloat}) = 100 * eps()

using FractalTools 
using GLMakie 

# Dataset 
f(x, y) = x^2 + y^2 + 1
Ω = [
    [-1.0, -1.0],
    [1.0, -1.0],
    [1.0, 1.0],
    [0.0, 0.5],
    [-1.0, 1.0],
]
Ω = @. map(item -> BigFloat(string(item)), Ω)  # Convert to string and bigfloat.
ipts = interiorpoints(f, Ω, 100)
bpts = boundarypoints(f, Ω, 10) 
pts = [bpts; ipts]
dataset = Dataset(pts, bpts)

# Construct interpolant 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluations
tipts = interiorpoints(Ω, 1000)
tbpts = boundarypoints(Ω, 10)
tpts = [tbpts; tipts]
ivals = map(pnt ->  interp(pnt...), tpts)
fvals = map(pnt ->  f(pnt...), tpts)
abserr = abs.(fvals - ivals) 
relerr = abserr ./ abs.(fvals) * 100

# Plots 
fig = Figure() 
ls11 = Axis3(fig[1, 1])
ls12 = Axis3(fig[1, 2])
ls21 = Axis3(fig[2, 1])
ls22 = Axis3(fig[2, 2])
plt11 = trisurf!(ls11, tpts, fvals, tbpts)
plt12 = trisurf!(ls12, tpts, ivals, tbpts)
plt21 = trisurf!(ls21, tpts, abserr, tbpts)
plt21 = trisurf!(ls22, tpts, relerr, tbpts)
fig 
