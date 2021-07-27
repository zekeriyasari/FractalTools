# This file includes IFS tools 

export  Transformation, IFS, Attractor, Tree, Line, Square, Fern, Sierpinski, DetAlg, RandAlg, dimension, contfactor 

#--- ------------------------------------- Transformation ------------------------------------------------ # 

"""
    $TYPEDEF

Affine transformation of the form ``w(x) = A x + b `` 

# Fields 
    $TYPEDFIELDS
"""
struct Transformation{T1<:AbstractMatrix, T2<:AbstractVector} 
    "Transformation matrix"
    A::T1 
    "Transformation vector"
    b::T2 
end 
(w::Transformation)(x) = w.A * x + w.b 

"""
    $TYPEDSIGNATURES

Returns dimension of of `w`.
"""
dimension(w::Transformation) = size(w.A, 1)

"""
    $TYPEDSIGNATURES

Returns contraction factor of `w` 
"""
contfactor(w::Transformation, p::Int=2) = norm(w.A, p)

#--- -------------------------------------------- IFS ---------------------------------------------------- # 

"""
    $TYPEDEF

Iterated function system.

# Fields 
    $TYPEDFIELDS
"""
struct IFS{T1<:AbstractVector{<:Transformation}, T2<:AbstractVector{<:Real}} 
    "Transformations"
    ws::T1 
    "Probilities of transformations"
    probs::T2 
    function IFS(ws::T1, probs::T2) where {T1, T2}
        sum(probs) ≈ 1 || error("Sum of probs must be 1")
        new{T1, T2}(ws, probs)
    end
end 

IFS(ws) = (n = length(ws); IFS(ws, ones(n) / n))

"""
    $TYPEDSIGNATURES

Returns dimension of `ifs` 
"""
dimension(ifs::IFS) = dimension(ifs.ws[1])

"""
    $TYPEDSIGNATURES

Returns contraction factor of `ifs` 
"""
contfactor(ifs::IFS) = maximum(contfactor.(ifs.ws))

"""
    $SIGNATURES

Returns the IFS of a line between [0, 1].
"""
Line() = IFS([
    Transformation(fill(1 / 2, 1, 1), fill(0, 1)), 
    Transformation(fill(1 / 2, 1, 1), fill(1 / 2, 1))
])

"""
    $SIGNATURES

Returns the IFS of a Sierpinski triangle
"""
Sierpinski() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [0; 0]),
    Transformation([0.5 0.0; 0. 0.5], [0; 1/2.]),
    Transformation([0.5 0.0; 0. 0.5], [1/2.; 1/2.])
    ], [1/3., 1/3., 1/3.])

"""
    $SIGNATURES

Returns the IFS of a unit square [0, 1] - [0, 1]
"""
Square() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [0., 0.]),
    Transformation([0.5 0.0; 0. 0.5], [0.5, 0.]),
    Transformation([0.5 0.0; 0. 0.5], [0., 0.5]),
    Transformation([0.5 0.0; 0. 0.5], [0.5, 0.5])
    ], [0.25, 0.25, 0.25, 0.25])

"""
    $SIGNATURES

Returns the IFS of a fern
"""
Fern() = IFS([
    Transformation([0 0; 0 0.16], [0.; 0.]),
    Transformation([0.85 0.04; -0.04 0.85],[0.; 1.6]),
    Transformation([0.2 -0.26; 0.23 0.22], [0.; 1.6]),
    Transformation([-0.15 0.28; 0.26 0.24], [0.; 0.44])
    ], [0.01, 0.85, 0.07, 0.07])

"""
    $SIGNATURES

Returns the IFS of a tree. 
"""
Tree() = IFS([
    Transformation([0 0; 0 0.5], [0.; 0.]),
    Transformation([0.42 -0.42; 0.42 0.42], [0.; 0.2]),
    Transformation([0.42 0.42; -0.42 0.42], [0.; 0.2]),
    Transformation([0.1 0; 0 0.1], [0.; 0.2])
    ], [0.05, 0.40, 0.40, 0.15])

#--- ------------------------------------- Attracttor ------------------------------------------------ # 

abstract type Algorithm end

"""
    $TYPEDEF

Deterministic algoirthm to compute attractor of an ifs. 
"""
struct DetAlg <: Algorithm end

"""
    $TYPEDEF

Random iteration algoirthm to compute attractor of an ifs.
"""
struct RandAlg <: Algorithm end

"""
    $TYPEDEF

Attractor of an ifs. 

# Fields 
    $TYPEDFIELDS
"""
mutable struct Attractor{T1<:IFS, T2<:Algorithm, T3} 
    "IFS of the attractor"
    ifs::T1 
    "Algorithm to compute attractor. May be `DetAlg` of `RandAlg`"
    alg::T2 
    "Generator of the attactor that generates points from the attractor"
    generator::T3
    "Size of the point chunks drawn from the generator"
    chunksize::Int 
    "State of the attractor: May be `:open` or `:closed`"
    state::Symbol 
end 

Attractor(ifs, alg, generator, chunksize=1) = Attractor(ifs, alg, generator, chunksize, :open)

getinitset(ifs) = [rand(dimension(ifs))]

"""
    $TYPEDSIGNATURES

Returns a chunk of points drawn from the attractor
"""
function Base.take!(atr::Attractor)
    if atr.state == :open 
        take!(atr.generator) 
    else 
        @warn "Attractor is closed"
    end 
end

# While an instance of `Attractor` is constructed, transient steps are taken so that when new points are drawn from the
# attractor, they are guarenteed to be drawn from the attractor. 
function transients(ifs, initset, numiter, posterrorbound)
    # Compute number of iterations
    ws = ifs.ws
    probs = ifs.probs 
    σ, index = findmax(contfactor.(ws))
    wm = ws[index]
    if numiter === nothing 
        x0 = initset[end]
        x1 = wm(x0)
        k = ceil(Int, (log(posterrorbound) + log(1 - σ) - log(norm(x1 - x0))) / log(σ))
    else 
        k = numiter 
    end 

    # Take transient steps
    weights = Weights(probs)
    xnew = only(initset)
    for i in 1 : k 
        xnew = xnew |> sample(ws, weights)
    end 
    [xnew] 
end 

# Worker function for random iteration algorithm 
function worker(alg::RandAlg, ifs, initset, channel, chunksize) 
    ws = ifs.ws 
    probs = ifs.probs 
    weights = Weights(probs)
    xi = initset[end] 
    while true  
        chunk = map(1 : chunksize) do i 
            xi = xi |> sample(ws, weights)
        end 
        put!(channel, chunk)
    end 
end 

# Worker function for deterministic algorithm
function worker(alg::DetAlg, ifs, set, channel) 
    ws = ifs.ws
    while true 
        set = vcat(map(w -> w.(set), ws)...)
        put!(channel, set)
    end
end 

function Base.setproperty!(atr::Attractor, name::Symbol, val) 
    name == :chunksize && updatetask!(atr, val)
    setfield!(atr, name, val) 
end 

Base.iterate(atr::Attractor, state...) = (take!(atr), nothing)

function updatetask!(atr, chunksize)
    initset = take!(atr.generator)
    close(atr.generator) 
    atr.generator = Channel{eltype(atr.generator)}(0)
    task = @async worker(atr.alg, atr.ifs, initset, atr.generator, chunksize)
    bind(atr.generator, task)
    atr
end 

Attractor(ifs; kwargs...) = Attractor(ifs, RandAlg(); kwargs...)

function Attractor(ifs, alg::RandAlg; initset=getinitset(ifs), numiter=nothing, posterrorbound=1e-10, chunksize=1, parallel=false)
    if parallel
        # TODO: Compute attractor in parallel 
    else 
        initset = transients(ifs, initset, numiter, posterrorbound)
        generator = Channel{typeof(initset)}(0) 
        task = @async worker(alg, ifs, initset, generator, chunksize)
        bind(generator, task)
        Attractor(ifs, alg, generator, chunksize)
    end 
end

function Attractor(ifs, alg::DetAlg; initset=getinitset(ifs), numiter=nothing, posterrorbound=1e-10, parallel=false)
    if parallel
        # TODO: Compute attractor in parallel 
    else 
        initset = transients(ifs, initset, numiter, posterrorbound)
        generator = Channel{typeof(initset)}(0) 
        task = @async worker(alg, ifs, initset, generator)
        bind(generator, task)
        Attractor(ifs, alg, generator, 0)
    end 
end

