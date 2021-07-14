# This file includes IFS tools 

export  Transformation, IFS, Attractor, Tree, Square, Fern, Sierpinski, DetAlg, RandAlg, dimension, contfactor 

#--- ------------------------------------- Transformation ------------------------------------------------ # 

struct Transformation{T1, T2} 
    A::T1 
    b::T2 
end 
(w::Transformation)(x) = w.A * x + w.b 

dimension(w::Transformation) = size(w.A, 1)
contfactor(w::Transformation, p::Int=2) = norm(w.A, p)

#--- -------------------------------------------- IFS ---------------------------------------------------- # 

struct IFS{T1, T2} 
    ws::T1 
    probs::T2 
    function IFS(ws::T1, probs::T2) where {T1, T2}
        sum(probs) ≈ 1 || error("Sum of probs must be 1")
        new{T1, T2}(ws, probs)
    end
end 

IFS(ws) = (n = length(ws); IFS(ws, ones(n) / n))

dimension(ifs::IFS) = dimension(ifs.ws[1])
contfactor(ifs::IFS) = maximum(contfactor.(ifs.ws))

Sierpinski() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [0; 0]),
    Transformation([0.5 0.0; 0. 0.5], [0; 1/2.]),
    Transformation([0.5 0.0; 0. 0.5], [1/2.; 1/2.])
    ], [1/3., 1/3., 1/3.])

Square() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [0., 0.]),
    Transformation([0.5 0.0; 0. 0.5], [0.5, 0.]),
    Transformation([0.5 0.0; 0. 0.5], [0., 0.5]),
    Transformation([0.5 0.0; 0. 0.5], [0.5, 0.5])
    ], [0.25, 0.25, 0.25, 0.25])

Fern() = IFS([
    Transformation([0 0; 0 0.16], [0.; 0.]),
    Transformation([0.85 0.04; -0.04 0.85],[0.; 1.6]),
    Transformation([0.2 -0.26; 0.23 0.22], [0.; 1.6]),
    Transformation([-0.15 0.28; 0.26 0.24], [0.; 0.44])
    ], [0.01, 0.85, 0.07, 0.07])

Tree() = IFS([
    Transformation([0 0; 0 0.5], [0.; 0.]),
    Transformation([0.42 -0.42; 0.42 0.42], [0.; 0.2]),
    Transformation([0.42 0.42; -0.42 0.42], [0.; 0.2]),
    Transformation([0.1 0; 0 0.1], [0.; 0.2])
    ], [0.05, 0.40, 0.40, 0.15])

#--- ------------------------------------- Attracttor ------------------------------------------------ # 

abstract type Algorithm end
struct DetAlg <: Algorithm end
struct RandAlg <: Algorithm end

mutable struct Attractor{T1, T2, T3} 
    ifs::T1 
    alg::T2 
    generator::T3 
    state::Symbol # :open, :closed
end 

Attractor(ifs, alg, generator) = Attractor(ifs, alg, generator, :open)

getinitset(ifs) = [rand(dimension(ifs))]

function Base.take!(atr::Attractor)
    if atr.state == :open 
        take!(atr.generator) 
    else 
        @warn "Attractor is closed"
    end 
end


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

function worker(alg::DetAlg, ifs, set, channel) 
    ws = ifs.ws
    while true 
        set = vcat(map(w -> w.(set), ws)...)
        put!(channel, set)
    end
end 

function Base.setproperty!(atr::Attractor, name::Symbol, val) 
    if name == :chunksize 
        updatetask!(atr, val)
    else 
        setfield!(atr, name, val) 
    end 
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
        Attractor(ifs, alg, generator)
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
        Attractor(ifs, alg, generator)
    end 
end



