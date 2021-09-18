# This file includes utility functions 

export combine, tovector, topoint, bigfloat, meshgrid, getline 

"""
    $SIGNATURES 

Combines `vals`.

# Example 
```julia 

julia> foo = [rand(2) for i in 1 : 3]
3-element Vector{Vector{Float64}}:
 [0.07304969216674784, 0.13922697837562814]
 [0.17391711119297937, 0.01724058232513692]
 [0.6244409832793123, 0.9377684091200582]

julia> bar = [rand(1) for i in 1 : 3]
3-element Vector{Vector{Float64}}:
 [0.694956323617361]
 [0.5868716964238885]
 [0.5625412213942036]

julia> combine(foo, bar) 
3-element Vector{Vector{Float64}}:
 [0.07304969216674784, 0.13922697837562814, 0.694956323617361]
 [0.17391711119297937, 0.01724058232513692, 0.5868716964238885]
 [0.6244409832793123, 0.9377684091200582, 0.5625412213942036]
```
"""
combine(vals::AbstractVector...) = [vcat(val...) for val in zip(vals...)]

tovector(pnt::AbstractPoint) = [pnt...]

topoint(vector::AbstractVector) = Point(vector...)

randomstring(n::Int=10) = String(rand('A' : 'z', n))

randommeshpath(n::Int=10) = joinpath(tempdir(), randomstring(n) * ".msh")

bigfloat(Ω::AbstractArray) = @. map(item -> BigFloat(string(item)), Ω)

function orientation(p1, p2, p3) 
    A = p2 - p1 
    B = p3 - p1 
    sgn = sign(A[1] * B[2] - A[1] - B[1])
    if sgn > 0 
        return :leftturn
    elseif sgn < 0 
        return :righturn 
    else 
        return :collinear
    end 
end 

function projection(p1, p2, p3) 
    v = p2 - p1 
    q = p3 - p1 
    (q' * v) / (v' * v) * v + p1
end 

function orthogonalunitvector(v)
    w = [1, -v[1] / v[2]]
    w / norm(w) 
end 

distance(p1, p2, p3) = norm(p3 - projection(p1, p2, p3))

function movepoint(p1, p2, p3, d=distance(p1, p2, p3) + 100eps())
    o0 = orientation(p1, p2, p3) 
    w = orthogonalunitvector(p2 - p1) 
    o1 = orientation(p1, p2, w) 
    (o0 == o1) && (w = -w)
    p3 + d * w 
end 

meshgrid(x, y) = ones(length(y)) * x', y * ones(length(x))'

getline(xi, yi, xf, yf) = x -> (yf - yi) / (xf - xi) * (x - xi) + yi
