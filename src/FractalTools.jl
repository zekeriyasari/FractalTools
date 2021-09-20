module FractalTools

using DocStringExtensions
using Distributed
using PyCall
using LinearAlgebra
using GeometryBasics
using StaticArrays

import Gmsh: gmsh

using Makie
import Makie: plot!, convert_arguments
import GeometryBasics: Ngon
import Base: show, display
import StatsBase: sample, Weights

const MAXPREC = 1024                # Maximum precision for BigFloat arithmetic 
const MAX_LOCATION_COUNT = 1000      # Maximum number of iteration for point location. See `locate` function
const POINT_LOCATION_PERTUBATION = 1e-6     # Amoun of perturbatin for point loation. See `locate` function

function __init__()
    global spt = pyimport_conda("scipy.spatial", "scipy")
end

include("testfunctions.jl")
include("utils.jl")
include("dataset.jl")
include("recipes.jl")
include("ifs.jl")
include("elton.jl")
include("interpolation.jl")
include("integration.jl")
include("integration_helpers.jl")
include("elton_integral_1D.jl")
include("elton_integral.jl")





end # module
