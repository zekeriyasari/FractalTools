module FractalTools

using DocStringExtensions
using Distributed
using PyCall
using LinearAlgebra
using GeometryBasics
using Clustering
using StaticArrays

using Makie
import Makie: plot!, convert_arguments
import GeometryBasics: Ngon
import Base: show, display
import StatsBase: sample, Weights

const MAXPREC = 1024    # Maximum precision for BigFloat arithmetic 

function __init__()
    global spt = pyimport_conda("scipy.spatial", "scipy")
end

include("testfunctions.jl")
include("datagenerators.jl")
include("recipes.jl")
include("ifs.jl")
include("elton.jl")
include("interpolation.jl")
include("integration.jl")


end # module
