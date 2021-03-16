# This file investigates integration error with respect to free variables 

using FractalTools 
using GeometryBasics
using Makie 

# Construct interpolation data 
f(x, y) = x^2 + y^2 + 1
ngon = Triangle(
    Point(BigFloat(0.), BigFloat(0.)), 
    Point(BigFloat(1.), BigFloat(0.)), 
    Point(BigFloat(0.5), BigFloat(1.)))
npts = 100
pts = getdata(f, ngon, npts)

# Compute errors 
fvals = 35 / 48
freevars = 0.001 : 0.001 : 0.025 
mse = map(freevars) do freevar 
    ival = integrate(pts, Interp2D(freevar))
    abs(fval - ival)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="Integration Error", title="2D Integration Error") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "integ2d_error.png"), fig)
display(fig)
