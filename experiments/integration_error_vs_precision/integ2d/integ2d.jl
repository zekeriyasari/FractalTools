# This file investigates integration error with respect to floating point precison. 

using FractalTools 
using Makie 

# Construct interpolation data 
f(x, y) = x^2 + y^2 + 1
precisions = 2 .^ (7 : 11)


# Specify theoretical value 
fval = 35 / 48

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

        # Compute integration.
        ival = integrate(pts, Interp2D(freevar))        

        # Compute errro 
        abs(fval - ival)
    end 
end 

# Plot mse 
fig = Figure() 
ax = fig[1, 1] = Axis(fig, xlabel="Precision", ylabel="Integration Error", title="2D Integration Error") 
stem!(ax, precisions, mse, color=:black)
save(joinpath(@__DIR__, "integ2d_error.png"), fig)
display(fig)
