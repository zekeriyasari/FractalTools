using FractalTools
using GLMakie 


# Data generation 
α = 1.
f(x, y) = FractalTools.paraboloid(x, y) 
n = 6
Ω = [[cos(θ), sin(θ)] for θ in range(0, 2π, step=2π/n)] * α
lc = 0.1 * α
pts = getdata(f, Ω, lc)

# Interpolation 
method = Interp2D(1e-3)
interp = interpolate(pts, method)

# Evaluations 
tpts = getdata(Ω, lc/2)
ivals = map(pnt -> interp(pnt...; d = 1e-9), tpts)
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

