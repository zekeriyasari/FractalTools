# This file investigates mean square error (MSE) with respect to free variables in hidden 1D interpolation 

using FractalTools 
using Makie 

# Construct interpolation data 
# f(x) = [FractalTools.sinusoid(x), FractalTools.sinusoid(x)]  # For regular data 
f(x) = [FractalTools.wen(x), FractalTools.weierstrass(x)]  # For irregular data 
line = [[0.], [1.]]
pts = getdata(f, line, 11)   

# Construct test data 
tpts = getdata(line, 101)
npts = length(tpts) 

# Compute errors 
fvals = getindex.(map(pnt -> f(pnt[1]), tpts), 1)
freevars = 0.001 : 0.01 : 0.5
mse = map(freevars) do freevar 
    interp = interpolate(pts, HInterp1D(fill(freevar, 2, 2)))
    ivals = getindex.(map(pnt -> interp(pnt[1]), tpts), 1)
    sum((fvals - ivals).^2) / npts 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="1D Hidden Interpolation MSE") 
stem!(ax, freevars, mse, color=:black)
save(joinpath(@__DIR__, "hinterp1d_mse.png"), fig)
display(fig)

