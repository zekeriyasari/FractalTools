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

function IFS(ws::AbstractVector{<:Transformation}, probs::AbstractVector{<:Real}, args...)  
    # Note: For the floating point numbers, aproximation(≈), instead of exact equal (==), should be considered
    sum(probs) ≈ 1 || throw(ArgumentError("Sum of probabilities must be 1."))
    initset = [rand(dimension(ws[1]))]
    attractor = Attractor(ws, probs, initset)
    IFS(ws, probs, attractor)
end
IFS(ws,args...) = (n = length(ws); IFS(ws, 1  / n * ones(n), args...))

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

A type signifying that deterministic algorithm is used when calculating the attractor of and IFS.
"""
struct DetAlg end

"""
    $TYPEDEF

A type signifying that random algorithm is used when calculating the attractor of and IFS.
"""
struct RandAlg end

"""
    $TYPEDEF

Attractor of `IFS` type 

# Fields

    $TYPEDFIELDS
"""
struct Attractor{S, R1, R2}
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
    chunk_size::Int
end

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
Attractor(ws, probs, initset; alg=DetAlg(), kwargs...) = typeof(alg) == DetAlg ? 
                                                   detalg(ws, probs, initset; kwargs...) : 
                                                   randalg(ws, probs,initset; kwargs...)

function attractor(ws, probs, initset; alg=DetAlg(), kwargs...)
    @warn "`attractor(ws, probs, initset; alg=DetAlg(), kwargs...)` has been deprecated. Use `Attractor(ws, probs, initset; alg=DetAlg(), kwargs...)` instead"
    Attractor(ws, probs, initset; alg=alg, kwargs...) 
end

"""
    $SIGNATURES

Computes the attractor of `ifs` with deterministic algorithm.`numiter` is number of iterations. (Defaults to 10). If
`parallel` is true, attractor is computed via parallel computation.
"""
function detalg(ws, probs, initset; numiter=10, parallel=false, chunk_size=10)
    copiedset = copy(initset)
    set = parallel ? 
          detalg_parallel(ws, copiedset, numiter) : 
          detalg_sequential(ws, copiedset, numiter)
    Attractor(DetAlg(), initset, set, numiter, parallel, chunk_size)
end

# Computes the attractor of an ifs via deterministic algorithm sequentially. 
function detalg_sequential(ws, set, numiter)
    for i in 1 : numiter
        set = vcat(map(w -> w.(set), ws)...)
    end
    set
end

# Computes the attractor of an ifs via deterministic algorithm in parallel. 
function detalg_parallel(ws, set, numiter)
    loadprocs()
    for i in 1 : numiter
        set = vcat(map(w -> pmap(w, set), ws)...)
    end
    set
end

"""
    $SIGNATURES

Computes the attractor of `ifs` with random algorithm.`numiter` is number of iterations. (Defaults to 100). `numtransient` is
the number of transient iterations. If `parallel` is true, attractor is computed via parallel computation. If
`placedependent` is true, the probabilties of the ifs are dependent on the coordinates `x`. This dependency `p(x)` is given
via the parameters `α` and `β` where p(x) = α x + β.
"""
function randalg(ws, probs, initset; numiter=100, numtransient=10, parallel=false, placedependent=false, α=nothing, β=nothing, allocated::Bool=false, chunk_size = 10, ϵ::Real=NaN)
    if parallel
        if placedependent
            transient = randalg_sequential_pd(ws, copy(initset), numtransient, probs, α, β, allocated)
            set = randalg_parallel_pd(ws, transient, numiter, probs, α, β, allocated)
        else
            transient = randalg_sequential(ws, copy(initset), numtransient, probs, allocated)
            set = randalg_parallel(ws, transient, numiter, probs, allocated)
        end
    else
        if placedependent
            set = randalg_sequential_pd(ws, copy(initset), numiter, probs, α, β, allocated)
        else
            if ϵ === NaN
                set = randalg_sequential(ws, copy(initset), numiter, probs, allocated)
            else 
                set = randalg_sequential(ws, copy(initset), numiter, probs, allocated, ϵ)   
            end  
        end
    end
    Attractor(RandAlg(), initset, set, numiter, parallel, chunk_size)
end

function _randalg_sequential(set::AbstractVector, ws, numiter::Int)

    weights = Weights(probs)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(ws, weights)
        xi = trfmi(xi)
        push!(set, xi)
    end
    return set
end

# Computes the attractor of an ifs via random algorithm sequentially. 
function _randalg_sequential(ch::AbstractChannel, xinit, ws, probs, numiter::Int)
    weights = Weights(probs)
    #TODO: Cancel the numiter, determine the post-error, K(contraction)
    for i = 1 : numiter
        trfmi = sample(ws, weights)
        xnew = trfmi(xinit)
        put!(ch, xnew)
        xinit = xnew
    end
    ch
end

function _randalg_sequential(ch::AbstractChannel, xinit, ws, probs, numiter=nothing, ϵ::Real=1e-8, chunk_size=10)
    # Compute initial set with a single point.
    if xinit === nothing
        n = size(ws[1].b)[1]
        xinit = rand(n)       
    else
        n = size(xinit)[1]
    end
    n = size(xinit)[1]

    # Compute number of iterations
    σ, index = findmax(contfactor.(ws))
    if numiter === nothing
        # Compute num_iter with respect to ϵ
        x1 = ws[index](xinit)
        _k = (log(ϵ) - log(norm(x1 - xinit))) / log(σ) + 1
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

    # Compute attractor
    chunk = zeros(n, chunk_size) 
    while true
        # chunk_size = take!(ch)
        # chunk_size === NaN && break
        for i = 1 : chunk_size
            trfmi = sample(ws, weights)
            xnew = trfmi(xnew)
            chunk[:,i] = xnew
        end
        put!(ch, chunk)
    end
end

function randalg_sequential(ws, set, numiter, probs, allocated::Bool=false)
    if allocated
        _randalg_sequential(set, ws, numiter, probs)
    else
        # NOTE: The set is assumed to have just a single initial point. If the set consists of more element, then we need a fix.
        ch = Channel(0)
        task = @async _randalg_sequential(ch, only(set), ws, numiter, probs)
        # task = @async _randalg_sequential(ch, ws, probs)
        bind(ch, task)
        ch
    end 
end

function randalg_sequential(ws, set, numiter, probs, ϵ::Real)
    ch = Channel(0)
    # ch = Channel{Vector{Float64}}()
    task = @async _randalg_sequential(ch, ws, probs, args...)
    bind(ch, task)
    ch
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

