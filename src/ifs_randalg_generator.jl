# This file includes IFS tools.

export IFS, attractor, RandAlg, dimension, contfactor


struct Transformation{T1<:AbstractMatrix{<:Real}, T2<:AbstractVector{<:Real}}
    "Transformation matrix"
    A::T1 
    "Transformation vector"
    b::T2 
end 

(w::Transformation)(x) = w.A * x + w.b


struct IFS{T1<:AbstractVector{<:Transformation}, T2<:AbstractVector{<:Real}}
    "Vector of transformations of IFS"
    ws::T1 
    "Vector of probabilities of IFS"
    probs::T2
end

function IFS(ws::T1, probs::T2) where {T1, T2} 
    # Note: For the floating point numbers, aproximation(≈),instead of exact equal (==), should be considered
    sum(probs) ≈ 1 || throw(ArgumentError("Sum of probabilitiesmust be 1."))
    new{T1, T2}(ws, probs) 
end
IFS(ws) = (n = length(ws); IFS(ws, 1  / n * ones(n)))


dimension(w::Transformation) = size(w.A,1)
contfactor(w::Transformation) = norm(w.A)
dimension(ifs::IFS) = dimension(ifs.ws[1])
contfactor(ifs::IFS) = maximum(contfactor.(ifs.ws))

# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [0]
A2 = reshape([1/2],1,1)
b2 = [1/2]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)


# function estimate_contraction_factor(ws)
#     _length = length(ws)
#     σ = zeros(_length)
#     for i = 1 : _length
#         σ[i] = norm(ws.A, inf)
#     end
# end

function _randalg_sequential(ch::AbstractChannel, ws, probs, xinit=nothing; chunk_size=1, num_iter=nothing, ϵ = 1e-10)
    if xinit=nothing
        n = size(ws[1].b)[1]
        xinit =rand(n)
    else
        n = size(xinit)[1]
    end
    σ, index= findmax(contfactor.(ws))
    if num_iter = nothing
        x1 = ws[index](xinit)
        _k = (log(ϵ) - log(norm(x1-xinit))) / log(σ) + 1
        k = Int(floor(_k))
    else
        k = num_iter
    end
    weights = Weights(probs)
    xnew = xinit
    for i = 1 : k
        trfmi = sample(ws, weights)
        xnew = trfmi(xnew)
    end
    chunk = zeros(n, chunk_size) 
    while true
        for i = 1 : chunk_size
            trfmi = sample(w, weights)
            xnew = trfmi(xnew)
            chunk[:,i] = xnew
        end
        put!(ch, chunk)
    end
end

function randalg_sequential_generator(ws, probs, args...)
    ch = Channel(0)
    task = @async _randalg_sequential(ch, ws, probs, args...)
    bind(ch, task)
    ch
end

struct DetAlg end
struct RandAlg end
struct Attractor{T, S, R1, R2}
    "IFS of Attractor"
    ifs::T
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
end

attractor(ifs, initset; alg=RandAlg(), kwargs...) = typeof(alg) == DetAlg ? 
                                                   randalg(ifs,initset; kwargs...) : 
                                                   detalg(ifs,initset; kwargs...)

function randalg(ifs, initset; numiter=100, numtransient=10, parallel=false, placedependent=false, α=nothing, β=nothing, allocated::Bool=false, args...)
    ws = ifs.ws
    probs = ifs.probs
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
            # set = randalg_sequential(ws, copy(initset), numiter, probs, allocated)
            # *************** check ************* #
            set = randalg_sequential_generator(ws, probs, args...)
        end
    end
    Attractor(ifs, RandAlg(), initset, set, numiter, parallel)
end








