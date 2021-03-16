# This file investigates mean square error (MSE) with respect to free variables in hidden 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x)                      
g(x) = cos(2π * x)
dx   = 0.1                              
xi   = 0.                               
xf   = 1.                               
x    = collect(xi : dx : xf)            
y    = f.(x)                            
z    = g.(x)
pts  = collect.(zip(x, y, z))     

# Construct test data 
xt   = collect(xi : 0.1dx : xf) 
npts = length(xt) 

# Compute errors 
fvals = f.(xt) 
freevars = 0.001 : 0.01 : 0.999
mse = map(freevars) do freevar 
    interp = interpolate(pts, HInterp1D(fill(freevar, 2, 2)))
    ivals = getindex.(interp.(xt), 1)
    sum((fvals - ivals).^2) / npts 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="1D Hidden Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "hinterp1d_mse.png"), fig)
display(fig)

