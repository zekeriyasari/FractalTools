# This file investigates 2D hiddden interpolation mean square error(MSE) with respect to inteprolation point number.

using FractalTools 
using Makie 

# Construct interpolation data 
f(x, y) = [
    x^2 + y^2 + 1,
    x^2 - y^2,
]
vtx = [
    BigFloat.([0.0, 0.0]), 
    BigFloat.([1.0, 0.0]), 
    BigFloat.([0.5, 1.0])
    ]
freevar = 0.001
npts    = 50 : 5 : 150
ntpts   = 2 * npts[end]

# Construct test data 
tpts = getdata(vtx, ntpts)

# Compute errors 
mse = map(npts) do npt
    @info npt 
    # Construct interpolation data 
    pts = getdata(f, vtx, npt)

    # Construct interpolant 
    interp = interpolate(pts, HInterp2D(fill(freevar, 2, 2)))

    # Compute error 
    fvals = getindex.(map(pnt -> f(pnt...), tpts), 1)
    ivals = getindex.(map(pnt -> interp(pnt...), tpts), 1)
    sum((fvals - ivals).^2) / length(tpts)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Number of Points", ylabel="MSE", title="2D Hidden Interpolation MSE") 
stem!(ax, npts, mse, color=:black)
save(joinpath(@__DIR__, "hinterp2d_mse.png"), fig)
display(fig)
