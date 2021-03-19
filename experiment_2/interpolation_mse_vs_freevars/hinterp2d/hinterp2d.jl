# This file investigates MSE with respect to free variables 

using FractalTools 
using GeometryBasics
using Makie 

# Construct interpolation data 
_f(x, y) = -20*exp(-0.2 * √(0.5*(x^2+y^2))) -
            exp(0.5*(cos(2π*x)+cos(2π*y))) + 
            MathConstants.e + 20 

# _f(x, y) = x^2 + y^2 + 1 

# _f(x, y) = sin(x) * cos(y)

# function f(x, y) 
#     val = _f(x, y)
#     [val, val]
# end


f(x, y) = [
    _f(x, y) , 
    x^2 - y^2 - 1
    ]
ngon = Triangle(
    Point(BigFloat(-5.), BigFloat(-5.)), 
    Point(BigFloat(5), BigFloat(-5.)), 
    Point(BigFloat(0), BigFloat(5.)))

# ngon = Triangle(
#     Point(BigFloat(0.), BigFloat(0.)), 
#     Point(BigFloat(1), BigFloat(0.)), 
#     Point(BigFloat(0.5), BigFloat(1.)))

# ngon = [BigFloat.([-5,5]),
#         BigFloat.([-5,-5]),
#         BigFloat.([5,-5]),
#         BigFloat.([5,5])
# ]
npts = 100
pts = getdata(f, ngon, npts)

# Construct test data 
tpts = getdata(ngon, npts)
ntpts = length(tpts)

# Compute errors 
fvals = getindex.(map(pt -> f(pt...), tpts), 1)
freevars = 0.001 : 0.001 : 0.05
mse = map(freevars) do freevar
    interp = interpolate(pts, HInterp2D(fill(freevar, 2, 2)))
    ivals = getindex.(map(pt -> interp(pt...), tpts), 1)
    sum((fvals - ivals).^2) / ntpts
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="2D Hidden Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "hinterp2d_mse.png"), fig)
display(fig)
