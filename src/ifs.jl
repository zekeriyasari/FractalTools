# This file includes IFS tools.
using LinearAlgebra
using StatsBase

export IFS, attractor, DetAlg, RandAlg, dimension, contfactor, Sierpinski, Fern, Square, Tree

export w, Transformation, randalg_sequential_for_generator, randalg_sequential_generator

""" 
    $(TYPEDEF) 
Affine transformation 

# Fields 
    $(TYPEDFIELDS)
"""
struct Transformation{T1<:AbstractMatrix{<:Real}, T2<:AbstractVector{<:Real}}
    "Transformation matrix"
    A::T1 
    "Transformation vector"
    b::T2 
end 

(w::Transformation)(x) = w.A * x + w.b

"""
    $SIGNATURES 

Returns dimension of `w`.
"""
dimension(w::Transformation) = size(w.A,1)

"""
    $SIGNATURES 

Returns contraction factor of `w`. Contraction factor is computed as the norm of `w.A`.
"""
contfactor(w::Transformation) = norm(w.A)


""" 
    $(TYPEDEF) 

Iterated fucntion sytem (IFS) 

# Fields 

    $(TYPEDFIELDS)
"""
struct IFS{T1<:AbstractVector{<:Transformation}, T2<:AbstractVector{<:Real},T3}
    "Vector of transformations of IFS"
    ws::T1 
    "Vector of probabilities of IFS"
    probs::T2
    "Attractor of the ifs"
    attractor::T3
end

# TODO: Document args and kwargs... 
function IFS(ws::AbstractVector{<:Transformation}, 
             probs::AbstractVector{<:Real}, 
             initset::AbstractVector{<:AbstractVector}=[rand(dimension(ws[1]))];
             kwargs...)  
    # Note: For the floating point numbers, aproximation(≈), instead of exact equal (==), should be considered
    sum(probs) ≈ 1 || throw(ArgumentError("Sum of probabilities must be 1."))
    attractor = Attractor(ws, probs, initset; kwargs...)
    IFS(ws, probs, attractor)
end
IFS(ws, args...) = (n = length(ws); IFS(ws, 1  / n * ones(n), args...))

""" 
    $SIGNATURES

Conctructs an IFS for Sierpinski triangle.
"""
Sierpinski() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [0; 0]),
    Transformation([0.5 0.0; 0. 0.5], [0; 1/2.]),
    Transformation([0.5 0.0; 0. 0.5], [1/2.; 1/2.])
    ], [1/3., 1/3., 1/3.])

"""
    $SIGNATURES

Constructs and IFS for a sqaure.
"""
Square() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [1.; 1.]),
    Transformation([0.5 0.0; 0. 0.5], [50.; 1.]),
    Transformation([0.5 0.0; 0. 0.5], [1.; 50.]),
    Transformation([0.5 0.0; 0. 0.5], [50.; 50.])
    ], [0.25, 0.25, 0.25, 0.25])

"""
    $SIGNATURES

Constructs and IFS for a fern.
"""
Fern() = IFS([
    Transformation([0 0; 0 0.16], [0.; 0.]),
    Transformation([0.85 0.04; -0.04 0.85],[0.; 1.6]),
    Transformation([0.2 -0.26; 0.23 0.22], [0.; 1.6]),
    Transformation([-0.15 0.28; 0.26 0.24], [0.; 0.44])
    ], [0.01, 0.85, 0.07, 0.07])

"""
    $SIGNATURES

Constructs and IFS for a fractal tree.
"""
Tree() = IFS([
    Transformation([0 0; 0 0.5], [0.; 0.]),
    Transformation([0.42 -0.42; 0.42 0.42], [0.; 0.2]),
    Transformation([0.42 0.42; -0.42 0.42], [0.; 0.2]),
    Transformation([0.1 0; 0 0.1], [0.; 0.2])
    ], [0.05, 0.40, 0.40, 0.15])

"""
    $SIGNATURES

Returns dimension of `ifs`.
"""
dimension(ifs::IFS) = dimension(ifs.ws[1])

"""
    $SIGNATURES

Returns the contraction factor of `IFS`.
"""
contfactor(ifs::IFS) = maximum(contfactor.(ifs.ws))

"""
    $TYPEDEF

Abstract algorithm See [`DetAlg`](@ref), [`RandAlg`](@ref)  
"""
abstract type AbstractAlgorithm end

"""
    $TYPEDEF

A type signifying that deterministic algorithm is used when calculating the attractor of and IFS.
"""
struct DetAlg <: AbstractAlgorithm end

"""
    $TYPEDEF

A type signifying that random algorithm is used when calculating the attractor of and IFS.
"""
struct RandAlg <: AbstractAlgorithm end

"""
    $TYPEDEF

Attractor of `IFS` type 

# Fields

    $TYPEDFIELDS
"""
mutable struct Attractor{S  <: AbstractAlgorithm, 
                         R1 <: AbstractVector{<:AbstractVector}, 
                         R2 <: Union{<:AbstractVector{<:AbstractVector}, AbstractChannel}}
    "Type of algorithm to be used to compute attractor(Options are DetAlg and RandAlg"
    alg::S
    "Initial set of attractor"
    initset::R1
    "Set of the attractor"
    set::R2
    "Number of iterations"
    numiter::Int
    "Sequential or parallel"
    parallel::Bool
    "Number of chunk size"
    chunksize::Int
end

# TODO: Document kwargs...
"""
    $SIGNATURES

Computes the attractor of `ifs`. If `alg` is of type `DetAlg`, the deterministic algorithm is used. If `alg` is of type
`RandAlg`, random algorithm is used. `kwargs` may include

* `numiter::Int` : Number of iterations to used to calcuate the attractor (defaults to 10)

* `numtransient::Int` : Number of transient iterations to used to calcuate a transient set. When the transient set is
  constructed, the computation of attractor is continued with distributed computation if `alg` is `RandAlg` and `parallel` is
  `true`. (defaults to 10)

* `parallel::Bool`: If  `true`, the attractor is computed using distrbuted computation. (defaults to false)

* `placedependent::Bool` : If `true`, place dependent attractor is computed if α and β is given accordingly. (default to
  false)

* `α::AbstractVector` : Place-dependent probility coefficient(defaults to nothing)

* `β::AbstractVector` : Place-dependent probility coefficient. (default to nothing)
"""
function Attractor(ws::AbstractVector{<:Transformation}, 
                   probs::AbstractVector{<:Real}, 
                   initset::AbstractVector{<:AbstractVector};
                   alg=DetAlg(), 
                   kwargs...) 
    if typeof(alg) == DetAlg 
        detalg(ws, probs, initset; kwargs...)
    else
        randalg(ws, probs, initset; kwargs...)
    end 
end 

function attractor(ws::AbstractVector{<:Tranformation}, 
                   probs::AbstractVector{<:Real}, 
                   initset::AbstractVector{<:AbstractVector}; 
                   alg=DetAlg(), 
                   kwargs...)
    msg = "`attractor(ws, probs, initset; alg=DetAlg(), kwargs...)` has been deprecated". 
    msg *= "Use `Attractor(ws, probs, initset; alg=DetAlg(), kwargs...)` instead"
    @warn msg 
    Attractor(ws, probs, initset; alg=alg, kwargs...) 
end

"""
    $TYPEDSIGNATURES

Computes the attractor of `ifs` with deterministic algorithm.`numiter` is number of iterations. (Defaults to 10). If
`parallel` is true, attractor is computed via parallel computation.
"""
function detalg(ws::AbstractVector{<:Transformation}, 
                probs::AbstractVector{<:Real}, 
                initset::AbsractVector{<:AbstractVector}; 
                numiter::Int=10, 
                parallel::Bool=false, 
                chunksize::Int=10)
    copiedset = copy(initset)
    set = if parallel
        detalg_parallel(ws, copiedset, numiter)
    else 
        detalg_sequential(ws, copiedset, numiter)
    end 
    Attractor(DetAlg(), initset, set, numiter, parallel, chunksize)
end

# Computes the attractor of an ifs via deterministic algorithm sequentially. 
function detalg_sequential(ws::AbstractVector{<:Transformation}, 
                           set::AbstractVector{<:AbstractVector}, 
                           numiter::Int)
    for i in 1 : numiter
        set = vcat(map(w -> w.(set), ws)...)
    end
    set
end

# Computes the attractor of an ifs via deterministic algorithm in parallel. 
function detalg_parallel(ws::AbstractVector{<:Transformation}, 
                         set::AbstractVector{<:AbstractVector}, 
                         numiter::Int)
    loadprocs()
    for i in 1 : numiter
        set = vcat(map(w -> pmap(w, set), ws)...)
    end
    set
end

"""
    $TYPEDSIGNATURES

Computes the attractor of `ifs` with random algorithm.`numiter` is number of iterations. (Defaults to 100). `numtransient` is
the number of transient iterations. If `parallel` is true, attractor is computed via parallel computation. If
`placedependent` is true, the probabilties of the ifs are dependent on the coordinates `x`. This dependency `p(x)` is given
via the parameters `α` and `β` where p(x) = α x + β.
"""
function randalg(ws::AbstractVector{<:Tranformation},           # IFS transformations 
                 probs::AbstractVector{<:Real},                 # IFS tranformation probabilities 
                 initset::AbstractVector{<:AbstractVector};     # Initial set to compute attractor 
                 numiter::Int=100,                              # Maximum iterations number 
                 posterrorbound::Real=NaN,                      # Posterior error bound
                 allocated::Bool=false,                         # If true, attractor points are allocated   
                 chunksize::Int=10,                             # Chunk size when attractor is processed (eg. elton integ.). 
                 numtransient::Int=10,                          # Number of transient steps before constructing attracttor
                 parallel::Bool=false,                          # If true, attractor is constructed using parallel computing 
                 placedependent::Bool=false,                    # If true, place dependent attractor is constructed. 
                 α=nothing,                                     # Place dependent transformation is in form p(x) = αx + β 
                 β=nothing)                                     # Place dependent transformation is in form p(x) = αx + β 
    if parallel
        # FIXME: The order of the argument positions must be checked.
        if placedependent
            transient = randalg_sequential_pd(ws, probs, copy(initset), numtransient, α, β, allocated)
            set = randalg_parallel_pd(ws, probs, transient, numiter, α, β, allocated)
        else
            transient = randalg_sequential(ws, probs, copy(initset), numtransient, allocated)
            set = randalg_parallel(ws, probs, transient, numiter, allocated)
        end
    else
        if placedependent
            set = randalg_sequential_pd(ws, probs, copy(initset), numiter, α, β, allocated)
        else
            if posterrorbound === NaN
                set = randalg_sequential(ws, probs, copy(initset), numiter, chunksize, allocated)
            else 
                set = randalg_sequential(ws, probs, copy(initset), numiter, posterrorbound, chunksize, allocated)   
            end  
        end
    end
    Attractor(RandAlg(), initset, set, numiter, parallel, chunksize)
end

function randalg_sequential(ws::AbstractVector{<:Tranformation}, 
                            probs::AbsractVector{<:Real}, 
                            set::AbstractVector{<:AbstractVector}, 
                            numiter::Int, 
                            allocated::Bool=false)
    if allocated
        # Compute whole attractor and allocate attractor points. 
        _randalg_sequential(ws, probs, set, numiter)
    else
        # Compute whole attractor and do not allocate attractor points. 
        channel = Channel(0)
        task = @async _randalg_sequential(ws, probs, channel, only(set), numiter)
        bind(channel, task)
        channel
    end 
end

# Compute whole attractor and without allocating attractor points via Channel-Task interface 
function randalg_sequential(ws::AbstractVector{<:Transformation}, 
                            probs::AbstractVector{<:Real},
                            set::AbstractVector{<:AbstractVector},
                            numiter::Int, 
                            posterrorbound::Real,
                            chunksize::Int)
    channel = Channel(0)
    task = @async _randalg_sequential(ws, probs, channel, only(set), numiter, posterrorbound, chunksize)
    bind(channel, task)
    channel
end

# Compute whole attractor and allocate attractor points. 
function _randalg_sequential(ws::AbstractVector{<:Transformation},
                             probs::AbstractVector{<:Real},
                             set::AbstractVector{<:AbstractVector}, 
                             numiter::Int)
    weights = Weights(probs)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(ws, weights)
        xi = trfmi(xi)
        push!(set, xi)      # Allocation
    end
    set
end


# Compute whole attractor and do not allocate attractor points. 
function _randalg_sequential(ws::AbstractVector{<:Transformation},
                             probs::AbstractVector{<:Real},
                             channel::AbstractChannel, 
                             xinit::AbstractVector{<:Real},
                             numiter::Int)
    weights = Weights(probs)
    for i = 1 : numiter
        trfmi = sample(ws, weights)
        xnew = trfmi(xinit)
        put!(channel, xnew)
        xinit = xnew
    end
    channel
end

function _randalg_sequential(ws::AbstractVector{<:Transformation},  
                             probs::AbstractVector{<:Real},
                             channel::AbstractChannel, 
                             xinit::AbstractVector{<:Real}, 
                             numiter::Int=nothing, 
                             posterrorbound::Real=1e-8, 
                             chunksize::Int=10)
    n = length(xinit)

    # # Compute number of iterations
    # σ, index = findmax(contfactor.(ws))
    # if numiter === nothing
    #     # Compute num_iter with respect to posterrorbound
    #     x1 = ws[index](xinit)
    #     _k = (log(posterrorbound) - log(norm(x1 - xinit))) / log(σ) + 1
    #     k = Int(ceil(_k))
    # else
    #     # Assign num_iter directly
    #     k = numiter
    # end

    # # Compute transients
    # weights = Weights(probs)
    # xnew = xinit
    # for i = 1 : k
    #     trfmi = sample(ws, weights)
    #     xnew = trfmi(xnew)
    # end

    # Compute attractor
    chunk = zeros(n, chunksize) 
    while true
        # chunksize = take!(channel)
        # chunksize === NaN && break
        for i = 1 : chunksize
            trfmi = sample(ws, weights)
            xnew = trfmi(xnew)
            chunk[:,i] = xnew
        end
        put!(channel, chunk)
    end
end

function transientsteps!(ws::AbstractVector{<:Transformation},
                         probs::AbstractVector{<:Real}, 
                         xinit::AbstractVector{<:Real},
                         numiter::Int, 
                         posterrorbound::Real) 
    # Compute number of iterations
    σ, index = findmax(contfactor.(ws))
    if numiter === nothing
        # Compute num_iter with respect to posterrorbound
        x1 = ws[index](xinit)
        # _k = (log(posterrorbound) - log(norm(x1 - xinit))) / log(σ) + 1
        _k = (log(posterrorbound) + log(1-σ) - log(norm(x1 - xinit))) / log(σ) + 1
        k = Int(ceil(_k))
    else
        # Assign num_iter directly
        k = numiter
    end

    # Compute transients
    weights = Weights(probs)
    xnew = xinit
    for i = 1 : k
        trfmi = sample(ws, weights)
        xnew = trfmi(xnew)
    end
end 

# Computes the attractor of an ifs via random algorithm sequentially with placedependent probabilties.
function randalg_sequential_pd(ws, set, numiter, probs, α, β, allocated::Bool=false)
    # TODO: Impletement no allocation method 
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(ws, Weights(probs))
        xi = trfmi(xi)
        probs = α * xi + β
    end
    set
end

# Computes the attractor of an ifs via random algorithm in parallel. 
function randalg_parallel(ws, set, numiter, probs, allocated::Bool=false)
    # TODO: Impletement no allocation method 
    weights = Weights(probs)
    loadprocs()
    vcat(pmap(process_chunk, [(ws, set, floor(Int, numiter / nworkers()), weights) for i =  1 : nworkers()])...)
end

# Computes the attractor of an ifs via random algorithm in parallel with placedependent probabilties.
function randalg_parallel_pd(ws, set, numiter, probs, α, β)
    loadprocs()
    vcat(pmap(process_chunk_pd, [(ws, set, floor(Int, numiter / nworkers()), probs, α, β) for i =  1 : nworkers()])...)
end

# `process_chunk` is the worker function that is used in all processes(both in master and worker process). when calculating
# the attractor if alg is `RandAlg` and `parallel` is true.
function process_chunk(ws_set_niter_weights)
    ws, set, niter, weights = ws_set_niter_weights
    xi = set[end]
    for i = 1 : niter
        wi = sample(ws, weights)
        xi = wi(xi)
        push!(set, xi)
    end
    set
end

# `process_chunk` is the worker function that is used in all processes(both in master and worker process). when calculating
# the attractor if alg is `RandAlg` and `parallel` is true with placedependent probabilties.
function process_chunk_pd(ws_set_niter_probs_alpha_beta)
    ws, set, niter, probs, α, β = ws_set_niter_probs_alpha_beta
    xi = set[end]
    for i = 1 : niter
        wi = sample(ws, Weights(probs))
        xi = wi(xi)
        probs = α * xi + β
        push!(set, xi)
    end
    set
end

# Load worker processes and load FractalTools to those worker processes.
function loadprocs(numprocs=Base.Sys.CPU_THREADS - 1 - nprocs())
    addprocs(numprocs)
    @everywhere @eval using FractalTools
end

