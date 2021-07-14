# This file compares Elton integral with quadrature integrral 

using FractalTools 
using Cubature 
using DataFrames 

# Define function and integration domain. 
f(x, y) = x^2 + y^2 
Ω = Attractor(Square()) # Fractal domain is unit square [0. 0] -> [1. 1]

# Compute integral using cubature 
cubaval = hcubature(p -> f(p...), [0., 0.], [1., 1.]) |> first 

# Compute integral using elton 
eltonval = elton(f, Ω, 100000)

# Absolute error 
relerr = abs(cubaval - eltonval) / abs(cubaval) * 100

# Display result 
df = DataFrame(
    name = ["Cubature Value", "Elton Value", "Relative Error"], 
    value = [cubaval, eltonval, relerr]
)
