# This file investigates mean square error (MSE) with respect to free variables in 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x)                      
dx   = 0.1                              
xi   = 0.                               
xf   = 1.                               
x    = collect(xi : dx : xf)            
y    = f.(x)                            
pts  = collect.(zip(x, y))     

# Construct test data 
xt   = collect(xi : 0.1dx : xf) 
npts = length(xt) 

# Compute errors 
fvals = f.(xt) 
freevars = 0.001 : 0.01 : 0.999
mse = map(freevars) do freevar 
    interp = interpolate(pts, Interp1D(freevar))
    ivals = interp.(xt) 
    sum((fvals - ivals).^2) / npts 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="1D Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "interp1d_mse.png"), fig)
display(fig)

