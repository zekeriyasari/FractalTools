# This file investigates the effect of number of points on the 2D hidden integration error. 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x, y) = [
    x^2 + y^2 + 1, 
    x^y - y^2 
]
vtx = [
    BigFloat.([0.0, 0.0]), 
    BigFloat.([1.0, 0.0]), 
    BigFloat.([0.5, 1.0])
    ]
freevar = 0.001
npts    = 50 : 5 : 200
ntpts   = 2 * npts[end]

# Compute errors 
fval = 35 / 48
mse = map(npts) do npt
    @info npt 
    # Construct interpolation data 
    pts = getdata(f, vtx, npt)

    # Construct interpolant 
    ival = integrate(pts, HInterp2D(fill(freevar, 2, 2)))

    # Compute error 
    abs(fval - ival)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Number of Points", ylabel="Error", title="2D Hidden Integration Error") 
stem!(ax, npts, mse, color=:black)
save(joinpath(@__DIR__, "hinteg2d_mse.png"), fig)
display(fig)
