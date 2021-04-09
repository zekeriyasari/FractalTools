# This file is for the attractor of Sierpinski 

using FractalTools
using Plots

# Construct ifs 
ifs = Sierpinski() 
initset = [rand(2)]
atr = attractor(ifs, initset, alg=RandAlg(), numiter=1000, allocated=false) 

# # Plot the attractor 
# scatter(getindex.(atr.set,1), getindex.(atr.set,2), markersize=3)
