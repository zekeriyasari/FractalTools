# This file inludes fractal integration tools. 

export integrate 

"""
    $SIGNATURES

Calculates numerical integration of `pts` using the fractal integration method `method` that may be 

* [`Interp1D`](@ref) 
* [`HInterp1D`](@ref) 
* [`Interp2D`](@ref) 
* [`HInterp2D`](@ref) 

"""
integrate(pts::AbstractVector{<:AbstractVector{<:Real}}, method::AbstractInterp) = 
    integrate(map(pnt -> Point(pnt...), pts), method)

function integrate(pts::PointVector, method::AbstractInterp) 
    tess = tessellate(pts, method)
    transforms = gettransforms(pts, method)
    val = evaluate(transforms, method, pts) 
end 

# Interp1D
function evaluate(transforms::AbstractVector{<:Transformation}, method::Interp1D, pts::PointVector{2})
    x  = getindex.(pts, 1) 
    (a11, a21, _, a22), (b1, b2) = extract(transforms)

    J   = abs.(a11)
    K1  = (x[end]^2 - x[1]^2) / 2 
    K2  = x[end] - x[1]
    
    num = sum(J .* (a21 * K1 + b2 * K2))
    denum = 1 - sum(J .* a22)
    num / denum
end 

# HInterp1D
function evaluate(transforms::AbstractVector{<:Transformation}, method::HInterp1D, pts::PointVector{3})
    x = getindex.(pts, 1)
    (a11, a21, a31, _, a22, a32, _, a23, a33), (b1, b2, b3) = extract(transforms)

    J   = abs.(a11)
    K1  = (x[end]^2 - x[1]^2) / 2 
    K2  = (x[end] - x[1])

    W11 = sum(a22 .* J);    W12 = sum(a23 .* J);    Λ1  = sum((a21 * K1 + b2 * K2) .* J);
    W21 = sum(a32 .* J);    W22 = sum(a33 .* J);    Λ2  = sum((a31 * K1 + b3 * K2) .* J);

    A = [1 - W11    -W12; 
         -W21       1 - W22]
    b = [Λ1, Λ2]
    I = A \ b 
    I[1]    # Hidden integration value I[2] is not returned. 
end 

# Interp2D
function evaluate(transforms::AbstractVector{<:Transformation}, method::Interp2D, pts::PointVector{3})
    (x1, y1), (x2, y2), (x3, y3) = getboundary(pts, method) 

    k11 = x2 - x1;      k12 = x3 - x1;      l1 = x1 
    k21 = y2 - y1;      k22 = y3 - y1;      l2 = y1
    
    JT = abs(k11 * k22 - k21 * k12)
    Δ1 = k11 + k12 + 3 * l1
    Δ2 = k21 + k22 + 3 * l2

    (a11, a21, a31, a12, a22, a32, _, _, a33), (b1, b2, b3) = extract(transforms)
    JL = abs.(a11 .* a22 - a21 .* a12)

    num = JT / 6 * sum(JL .* (a31 * Δ1 + a32 * Δ2 + 3 * b3))
    denum = 1 - sum(JL .* a33)
    num / denum
end 

# HInterp2D
function evaluate(transforms::AbstractVector{<:Transformation}, method::HInterp2D, pts::PointVector{4})
    (x1, y1), (x2, y2), (x3, y3) = getboundary(pts, method) 
    k11 = x2 - x1;      k12 = x3 - x1;      l1 = x1 
    k21 = y2 - y1;      k22 = y3 - y1;      l2 = y1
    
    JT = abs(k11 * k22 - k21 * k12)
    Δ1 = k11 + k12 + 3 * l1
    Δ2 = k21 + k22 + 3 * l2

    (a11, a21, a31, a41, a12, a22, a32, a42, _, _, a33, a43, _, _, a34, a44), (b1, b2, b3, b4) = extract(transforms)
    JL = abs.(a11 .* a22 - a21 .* a12)

    W11 = sum(a33 .* JL);    W12 = sum(a34 .* JL);    Λ1  = JT / 6 * sum((a31 * Δ1 + a32 * Δ2 + 3 * b3) .* JL);
    W21 = sum(a43 .* JL);    W22 = sum(a44 .* JL);    Λ2  = JT / 6 * sum((a41 * Δ1 + a42 * Δ2 + 3 * b4) .* JL);

    A = [1 - W11    -W12; 
         -W12       1 - W22]
    b = [Λ1, Λ2]
    I = A \ b 
    I[1]    # Hidden integration value I[2] is not returned. 
end 

function extract(transforms::AbstractVector{<:Transformation})
    As = getfield.(transforms, :A)
    bs = getfield.(transforms, :b)
    A = map(idx -> getindex.(As, idx...), [(i, j) for i in 1 : size(As[1], 1), j in 1 : size(As[1], 2)])
    b = map(idx -> getindex.(bs, idx), 1 : length(bs[1]))
    (A, b)
end 
