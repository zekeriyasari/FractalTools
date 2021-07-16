using FractalTools
using Cubature
using Plots 

# Define function 
f(x, y) = x^2 + y^2

# Define ranges 
cubaval = hcubature(x -> f(x...), [0, 0], [1., 1]) |> first 

# An example of one dimensional ifs
Ω = Attractor(Square(), chunksize=10)

pts = [take!(Ω) for i in 1 : 100]
plt = plot() 
for pt in pts 
    scatter!(plt, getindex.(pt,1), getindex.(pt, 2), markersize=2)
end 
plt 

# Evaluate integral 
elton_integral(f, indicator_measure, Ω, (), ([0.5, 0.5], 0.5), 1e-6, 1e8)
