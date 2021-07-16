using FractalTools

# Define function 
f(x,k) = k
f2(x,w=1) = sin(w * x)

# An example of one dimensional ifs
Ω = Attractor(IFS([
    Transformation(reshape([1/2], 1, 1), [0.]), 
    Transformation(reshape([1/2], 1, 1), [1 / 2]), 
]), chunksize=10)

# Evaluate integral 
elton_integral(f, indicator_measure, Ω, (1,), ([0.5], 0.25), 1e-10, 1e7)
