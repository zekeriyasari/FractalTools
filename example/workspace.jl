using FractalTools
using GLMakie 

using ConcaveHull

# Data generation 
α = 10.
f(x, y) = FractalTools.ackley(x, y) 
n = 6
Ω = [[cos(θ), sin(θ)] for θ in range(0, 2π, step=2π/n)] * α
lc = 0.2 * α
pts = getdata(f, Ω, lc)

Ω = [
    [0., 0.],
    [1., 0.],
    [1., 1.],
    [0.5, 0.5],
    [0., 1.],
]
lc = 0.1
pts = getdata(Ω, lc)
hull = concave_hull(pts, 100)

# Interpolation 
method = Interp2D(1e-3)
interp = interpolate(pts, method, niter=10)

# Evaluations 
tpts = getdata(Ω, lc/4)
ivals = map(pnt -> interp(pnt...), tpts)
fvals = map(pnt -> f(pnt...), tpts)
abserr = abs.(fvals - ivals)
relerr = abserr ./ abs.(fvals) * 100

# Plots 
msh = tomesh(tpts)
fig = Figure() 
ls = [Axis3(fig[i, j]) for i in 1 : 2, j in 1 : 2]
trisurf!(ls[1, 1], tpts, fvals)
trisurf!(ls[1, 2], tpts, ivals)
trisurf!(ls[2, 1], tpts, abserr)
trisurf!(ls[2, 2], tpts, relerr)
for li in ls 
    wireframe!(li,msh) 
end 
fig 

# # Plot 
# msh = tomesh(pts) |> project
# fig, ax, plt = scatterlines(getindex.(Ω, 1), getindex.(Ω, 2))
scatter(getindex.(pts, 1), getindex.(pts, 2))
scatterlines!(getindex.(hull.vertices, 1), getindex.(hull.vertices, 2))
# wireframe!(msh)
# scatter!(getindex.([pm], 1), getindex.([pm], 2), color=:red)
# isempty(ppts) ||  scatter!(getindex.(ppts, 1), getindex.(ppts, 2), color=:red)
# @info length(ppts)
# fig 
