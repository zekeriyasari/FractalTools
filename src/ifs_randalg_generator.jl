using LinearAlgebra
using StatsBase

export IFS, attractor, RandAlg, dimension, contfactor

#TODO : Affine Transformation
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

#TODO : \sigma  of the IFS's parameter

function _randalg_sequential(ch::AbstractChannel, ws, probs, xinit=nothing; num_iter=nothing, chunk_size=10, ϵ = 1e-8)
    # Compute initial set with a single point.
    if xinit === nothing
        n = size(ws[1].b)[1]
        xinit = rand(n)       
    else
        n = size(xinit)[1]
    end

    # Compute number of iterations
    σ, index = findmax(contfactor.(ws))
    if num_iter === nothing
        # Compute num_iter with respect to ϵ
        x1 = ws[index](xinit)
        _k = (log(ϵ) - log(norm(x1 - xinit))) / log(σ) + 1
        k = Int(ceil(_k))
    else
        # Assign num_iter directly
        k = num_iter
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

# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [0]
A2 = reshape([1/2],1,1)
b2 = [1/2]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)

generator = randalg_sequential_generator(ifs.ws, ifs.probs)

take!(generator)



