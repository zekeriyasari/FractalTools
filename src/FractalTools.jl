module FractalTools

using DocStringExtensions
using Distributed
using PyCall
using LinearAlgebra
using GeometryBasics
using Clustering
using StaticArrays

import AbstractPlotting
import AbstractPlotting: @recipe
import GeometryBasics: Ngon
import Base: show, display
import StatsBase: sample, Weights

function __init__()
    global spt = pyimport_conda("scipy.spatial", "scipy")
end

include("datagenerators.jl")
include("recipes.jl")
include("ifs.jl")
include("interpolation.jl")
include("integration.jl")

end # module
