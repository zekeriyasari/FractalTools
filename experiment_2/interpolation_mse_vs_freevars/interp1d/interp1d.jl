# This file investigates mean square error (MSE) with respect to free variables in 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
f = FractalTools.wen            # For irregular data 
# f = FractalTools.sinusoid     # For regular data 
line = [[0.], [1.]]
pts = getdata(f, line, 11) 

# Construct test data 
tpts = getdata(line, 101)
npts = length(tpts) 

# Compute errors 
fvals = map(pnt -> f(pnt[1]), tpts)
freevars = 0.001 : 0.01 : 0.999
mse = map(freevars) do freevar 
    interp = interpolate(pts, Interp1D(freevar))
    ivals = map(pnt -> interp(pnt[1]), tpts)
    sum((fvals - ivals).^2) / npts 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="1D Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "interp1d_mse.png"), fig)
display(fig)

