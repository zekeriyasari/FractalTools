# This file includes and example file for HInterp2D interpolation 

using FractalTools 
using GeometryBasics
using GLMakie 

# Generate data 
f(x, y) = [
    x^2 + y^2 + 1, 
    x^2 - y^2
    ]
Ω = uniformdomain(4, BigFloat)
N = 100
pts = getdata(f, Ω, N)

# Construct intepolant 
f̃ = interpolate(pts, HInterp2D(0.01 * ones(2,2)))

# Construct test data 
tpts = getdata(Ω, N)

# Evaluations 
fvals = map(p -> f(p...) |> first, tpts)
ivals = map(p -> f̃(p...) |> first, tpts)
evals = abs.(fvals - ivals) ./ abs.(fvals) * 100

# Plots 
fig = Figure() 
ls1 = Axis3(fig[1, 1])
scatter!(ls1, project(pts))
wireframe!(ls1, tomesh(project(pts)))
scatter!(ls1, project(pts, 2))
wireframe!(ls1, tomesh(project(pts, 2)))

ls2 = Axis3(fig[1, 2])
trisurf!(ls2, tpts, fvals)
wireframe!(ls2, tomesh(project(pts, 2)))

ls3 = Axis3(fig[2, 1])
trisurf!(ls3, tpts, ivals)
wireframe!(ls3, tomesh(project(pts, 2)))

ls4 = Axis3(fig[2, 2])
trisurf!(ls4, tpts, evals)
wireframe!(ls4, tomesh(project(pts, 2)))

fig 
