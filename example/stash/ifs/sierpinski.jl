# This file is for the attractor of Sierpinski 

using FractalTools
using Plots
using BenchmarkTools

# Construct ifs 
ifs = Sierpinski() 
initset = [rand(2)]
atr = attractor(ifs, initset, alg=RandAlg(), numiter=Int(1e10), allocated=false) 

function mysum(atr)
    total = 0 
    for point in atr.set 
        total += sum(point)
    end
    total 
end 

mysum(atr)

map(sum, atr.set)

# ch = Channel(0)
# task = @async FractalTools._randalg_sequential(ch, only(initset), ifs.ws, 1000, ifs.probs)

# # Plot the attractor 
# scatter(getindex.(atr.set,1), getindex.(atr.set,2), markersize=3)
