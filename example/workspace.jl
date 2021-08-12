# This file includes example surface interpolation performance

using FractalTools
using GLMakie 

# Dataset 
f(x, y) = x^2 + y^2 + 1
Ω = [
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
lc = 0.1
dataset = Dataset(f, Ω, lc)

# Perform interpolation 
method = Interp2D(0.001)
interp = interpolate(dataset, method)

# Evaluate interpolation 
tpts, tbpts = getdata(Ω, lc/2)
ivals = map(pnt -> interp(pnt...), tpts)
fvals = map(pnt -> f(pnt...), tpts)
abserr = abs.(fvals - ivals) 
relerr = abserr ./ abs.(fvals) * 100

# Plots 
tess = Tessellation(tpts, tbpts)
msh = tomesh(tpts, tbpts, tess)
fig = Figure()
ls = [LScene(fig[i,j]) for i in 1 : 2, j in 1 : 2]
trisurf!(ls[1, 1], tpts, fvals, tbpts)
trisurf!(ls[1, 2], tpts, ivals, tbpts)
trisurf!(ls[2, 1], tpts, abserr, tbpts)
trisurf!(ls[2, 2], tpts, relerr, tbpts)
for lsi in ls 
    wireframe!(lsi, msh)
end 
fig
