using ConcaveHull
using Meshes 
using FractalTools

pts = [[th, -sin(th)] .+ (0.4*rand(2) .- 0.2) for th in range(0, stop=2pi, length=5000)]

# Concave hull 
hpts = concave_hull(pts, 100).vertices .|> collect
ctr = Ngon(Point.(hpts)) |> centroid |> coordinates |> collect
_hpts = map(item -> item - ctr, hpts)
cpts = angle.(getindex.(_hpts, 1) + getindex.(_hpts, 2) * 1im ) 
ds = [angle => i for (i, angle) in enumerate(cpts)]
sds = sort(ds, by = x -> x[1])
Ω = hpts[last.(sds)]
# Ω = @. map(item -> BigFloat(string(item)), Ω) 

# Dataset 
f(x, y) = x^2 + y^2 + 1
ipts = interiorpoints(f, Ω, 100)
bpts = boundarypoints(f, Ω, 2) 
pts = [bpts; ipts]
dataset = Dataset(pts, bpts)

# Construct interpolant 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluations
tipts = interiorpoints(Ω, 100)
tbpts = Ω
tpts = [tbpts; tipts]
# tpts = tipts
ivals = map(pnt ->  interp(pnt...), tipts)
fvals = map(pnt ->  f(pnt...), tpts)
abserr = abs.(fvals - ivals) 
relerr = abserr ./ abs.(fvals) * 100

# Plots 

msh = tomesh(combine(tpts, fvals), tbpts) |> project 

fig = Figure() 
ls11 = Axis3(fig[1, 1])
ls12 = Axis3(fig[1, 2])
ls21 = Axis3(fig[2, 1])
ls22 = Axis3(fig[2, 2])
plt11 = trisurf!(ls11, tpts, fvals, tbpts)
plt12 = trisurf!(ls12, tpts, ivals, tbpts)
plt21 = trisurf!(ls21, tpts, abserr, tbpts)
plt21 = trisurf!(ls22, tpts, relerr, tbpts)
for ls in [ls11, ls12, ls21, ls22]
    wireframe!(ls, msh, linewidth=2)
end 
fig 

