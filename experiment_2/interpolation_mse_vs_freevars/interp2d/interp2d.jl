# This file investigates MSE with respect to free variables 

using FractalTools 
using GeometryBasics
using Makie 

# Construct interpolation data 
f = FractalTools.ackley
α = 10
ngon = [BigFloat.([-1., -1.]), BigFloat.([1., -1]), BigFloat.([0., 1.])] * α
npts = 100
pts = getdata(f, ngon, npts)

# Construct test data 
tpts = getdata(ngon, npts)
ntpts = length(tpts)

# Compute errors 
fvals = map(pt -> f(pt...), tpts) 
freevars = 0.001 : 0.001 : 0.05
mse = map(freevars) do freevar 
    interp = interpolate(pts, Interp2D(freevar))
    ivals = map(pt -> interp(pt...), tpts)
    sum((fvals - ivals).^2) / ntpts
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="2D Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "interp2d_mse.png"), fig)
display(fig)

# using FractalTools 
# using GeometryBasics
# using Makie 

# # Construct interpolation data 
# f(x, y) = -20*exp(-0.2 * √(0.5*(x^2+y^2))) -
#             exp(0.5*(cos(2π*x)+cos(2π*y))) + 
#             MathConstants.e + 20 

# # f(x, y) = x^2 + y^2 + 1

# # f(x, y) = sin(x) * cos(y)
# ngon = Triangle(
#     Point(BigFloat(-5.), BigFloat(-5.)), 
#     Point(BigFloat(5), BigFloat(-5.)), 
#     Point(BigFloat(0), BigFloat(5.)))

# # ngon = [BigFloat.([-5,5]),
# #         BigFloat.([-5,-5]),
# #         BigFloat.([5,-5]),
# #         BigFloat.([5,5])
# # ]
# npts = 100
# pts = getdata(f, ngon, npts)

# # Construct test data 
# tpts = getdata(ngon, npts)
# ntpts = length(tpts)

# # Compute errors 
# fvals = map(pt -> f(pt...), tpts) 
# freevars = 0.001 : 0.001 : 0.05
# mse = map(freevars) do freevar 
#     interp = interpolate(pts, Interp2D(freevar))
#     ivals = map(pt -> interp(pt...), tpts)
#     sum((fvals - ivals).^2) / ntpts
# end 

# # Plot mse 
# fig = Figure() 
# ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="2D Interpolation MSE") 
# stem!(ax, freevars, mse, color=:black)
# save(joinpath(@__DIR__, "interp2d_mse.png"), fig)
# display(fig)
