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

include("datagenerators.jl")
include("recipes.jl")
include("ifs.jl")
include("interpolation.jl")
include("integration.jl")

# include("ifs/transformations.jl")
# include("ifs/ifs.jl")
# include("ifs/prototypes.jl")
# include("interpolation/fif.jl")
# include("interpolation/hiddenfif.jl")
# include("interpolation/fis.jl")
# include("interpolation/hiddenfis.jl")
# include("integration/integration1d.jl")
# include("integration/integration2d.jl")

function __init__()
    global spt = pyimport_conda("scipy.spatial", "scipy")
end

end # module
