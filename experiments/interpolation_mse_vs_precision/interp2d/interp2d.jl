# This file investigates MSE with respect to free variables 

using FractalTools 
using GeometryBasics
using Makie 

# Construct interpolation data 
f(x, y) = x^2 + y^2 + 1
ngon = [
    [0.0, 0.0], 
    [1.0, 0.0], 
    [0.5, 1.0]
]
npts = 100
pts = getdata(f, ngon, npts)

# Construct test data 
tpts = getdata(ngon, npts)
ntpts = length(tpts)

# Compute errors 
fvals = map(pt -> f(pt...), tpts) 
freevars = 0.001 : 0.001 : 0.025 
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
