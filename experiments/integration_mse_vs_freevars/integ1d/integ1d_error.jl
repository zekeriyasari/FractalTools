# This file investigates integration error with respect to free variables in 1D integration 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x) = sin(2π * x) + 1            
dx   = 0.1                              
xi   = 0.                               
xf   = 1.                               
x    = collect(xi : dx : xf)            
y    = f.(x)                            
pts  = collect.(zip(x, y))     

# Compute errors 
fval = (cos(2π * xi) - cos(2π * xf)) / (2π) + (xf - xi)
freevars = 0.001 : 0.01 : 0.999
mse = map(freevars) do freevar 
    ival = integrate(pts, Interp1D(freevar))
    abs(fval - ival)
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="Integration Error", title="1D Integration Error") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "integ1d_error.png"), fig)
display(fig)

