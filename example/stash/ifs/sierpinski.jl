# This file is for the attractor of Sierpinski 

using FractalTools
using Plots
using BenchmarkTools
using Plots
# Construct ifs 
ifs = Sierpinski() 
initset = [rand(2)]
atr = attractor(ifs, initset, alg=RandAlg(), numiter=Int(1e4), allocated=true) 

# Plot the attractor 
scatter(getindex.(atr.set,1), getindex.(atr.set,2), markersize=2)

# ch = Channel(0)
# task = @async FractalTools._randalg_sequential(ch, only(initset), ifs.ws, 1000, ifs.probs)

