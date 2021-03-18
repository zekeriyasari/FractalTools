# This file investigates MSE with respect to free variables 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x, y) = x^2 + y^2 + 1
precisions = 2 .^ (7 : 11)

npts = 100
mse = map(precisions) do prec 
    # Set precision 
    setprecision(prec) do 
        freevar = BigFloat(0.001)

        # Construct interpolation data 
        ngon = [
            BigFloat.([0.0, 0.0]), 
            BigFloat.([1.0, 0.0]), 
            BigFloat.([0.5, 1.0])
        ]
        pts = getdata(f, ngon, npts)

        # Construct interpolant 
        interp = interpolate(pts, Interp2D(freevar))        

        # Construct test data 
        tpts = getdata(ngon, npts)
        ntpts = length(tpts)
        @info tpts[1][1].prec        
        @info interp.ifs.ws[1].A[1, 1].prec

        # Compute MSE 
        fvals = map(pt -> f(pt...), tpts) 
        ivals = map(pt -> interp(pt...), tpts)
        sum((fvals - ivals).^2) / ntpts
    end 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Free Variable", ylabel="MSE", title="2D Interpolation MSE") 
stem!(ax, precisions, mse, color=:black)
save(joinpath(@__DIR__, "interp2d_mse.png"), fig)
display(fig)
