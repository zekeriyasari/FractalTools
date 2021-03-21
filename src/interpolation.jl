# This includes interpolation methods 

export Interp1D, Interp2D, HInterp1D, HInterp2D, interpolate, project

const PointVector{Dim} = AbstractVector{<:AbstractPoint{Dim, T}} where {T}
const Tessellation = Union{<:LineString, <:PyObject}

abstract type AbstractInterp end
abstract type AbstractCurveInterp   <: AbstractInterp end       # One dimensional interpolations 
abstract type AbstractSurfaceInterp <: AbstractInterp end       # Two dimensional interpolations

"""
    $TYPEDEF

One dimensional fractal interpolation. IFS is defined as
```math 
    w_n(x, y) = \\begin{bmatrix} 
        a_{11, n} & 0 \\
        a_{21, n} & a_{22, n} 
    \\end{bmatrix} 
    \\begin{bmatrix} 
        b_{1, n} \\
        b_{2, n}
    \\end{bmatrix} \\quad n = 1, 2, \\ldots, N 
```
Here ``a_{22, n}, n = 1, 2, \\ldots, N`` are the free variables. While constructing an `Interp1D`, the free varaibles can be 
specfied as either a single real value or a vector of real values. If a single real value is specified, it assumed that 
``a_{22} = a_{22, n}, n = 1, 2, \\ldots, N``.

# Fields 

    $TYPEDFIELDS 
"""
struct Interp1D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: AbstractCurveInterp
    """Free variables"""
    freevars::T 
end 

"""
    $TYPEDEF

Two dimensional fractal interpolation. IFS is defined as
```math 
    w_n(x, y) = \\begin{bmatrix} 
        a_{11, n}   & a_{12, n}     & 0         \\
        a_{21, n}   & a_{22, n}     & 0         \\ 
        a_{31, n}   & a_{32, n}     & a_{33, n} \\  
    \\end{bmatrix} 
    \\begin{bmatrix} 
        b_{1, n} \\
        b_{2, n} \\ 
        b_{3, n} \\ 
    \\end{bmatrix} \\quad n = 1, 2, \\ldots, N 
```
Here ``a_{33,n}, n = 1, 2, \\ldots, N` are the free variables. While constructing an `InterpdD`, the free variables can be 
specfied as either a single real value or a vector of real values. If a single real value is specified, it assumed that 
``a_{33, n}, n = 1, 2, \\ldots, N``.

# Fields 

    $TYPEDFIELDS 
"""
struct Interp2D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: AbstractSurfaceInterp
    """Free variables"""
    freevars::T
end 


"""
    $TYPEDEF

One dimensional hidden fractal interpolation. IFS is defined as
```math 
    w_n(x, y) = \\begin{bmatrix} 
        a_{11, n}   & 0             & 0         \\
        a_{21, n}   & a_{22, n}     & a_{23, n} \\
        a_{31, n}   & a_{32, n}     & a_{33, n} \\ 
    \\end{bmatrix} 
    \\begin{bmatrix} 
        b_{1, n} \\ 
        b_{2, n} \\ 
        b_{3, n} \\ 
    \\end{bmatrix} \\quad n = 1, 2, \\ldots, N 
```
Here 
```math 
    \\begin{bmatrix}
        a_{22, n} & a_{23, n} \\ 
        a_{32, n} & a_{33, n}  
    \\end{bmatrix} 
```
are the free variables. While constructing an `HInterp1D`, the free variables can be specfied as either a single real matrix 
or a vector of real matrices. If a single real matrix is specified, it assumed that 
```math 
    \\begin{bmatrix}
        a_{22} & a_{23} \\ 
        a_{32} & a_{33}  
    \\end{bmatrix} 
    \\begin{bmatrix}
        a_{22, n} & a_{23, n} \\
        a_{32, n} & a_{33, n}  
    \\end{bmatrix} 
    \\quad \\forall n = 1, 2, \\ldots, N
```

# Fields 

    $TYPEDFIELDS 
"""
struct HInterp1D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: AbstractCurveInterp
    """Free variables"""
    freevars::T 
end 


"""
    $TYPEDEF

Two dimensional hidden fractal interpolation. IFS is defined as
```math 
    w_n(x, y) = \\begin{bmatrix} 
        a_{11, n}   & 0             & 0             & 0         \\
        a_{21, n}   & a_{22, n}     & 0             & 0         \\ 
        a_{31, n}   & a_{32, n}     & a_{33, n}     & a_{34, n} \\ 
        a_{41, n}   & a_{42, n}     & a_{43, n}     & a_{44, n} \\ 
    \\end{bmatrix} 
    \\begin{bmatrix} 
        b_{1, n} \\ 
        b_{2, n} \\
        b_{3, n} \\ 
        b_{4, n} \\ 
    \\end{bmatrix} \\quad n = 1, 2, \\ldots, N 
```
Here 
```math 
    \\begin{bmatrix}
        a_{33, n} & a_{34, n} \\ 
        a_{43, n} & a_{44, n}  
    \\end{bmatrix} 
```
are the free variables. While constructing an `HInterp2D`, the free variables can be specfied as either a single real matrix 
or a vector of real matrices. If a single real matrix is specified, it assumed that 
```math 
    \\begin{bmatrix}
        a_{33} & a_{34} \\ 
        a_{43} & a_{44}  
    \\end{bmatrix} 
    \\begin{bmatrix}
        a_{33, n} & a_{34, n} \\ 
        a_{43, n} & a_{44, n}  
    \\end{bmatrix} 
    \\quad \\forall n = 1, 2, \\ldots, N
```

# Fields 

    $TYPEDFIELDS 
"""
struct HInterp2D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: AbstractSurfaceInterp
    """Free variables"""
    freevars::T 
end 

# Note: For plotting recipes, Interpolant should a subtype of Function
"""
    $TYPEDEF

Interpolation object. 

# Fields 

    $TYPEDFIELDS
"""
struct Interpolant{T1<:IFS, T2, T3<:AbstractInterp} <: Function 
    """IFS """
    ifs::T1 
    """Interpolation function"""
    itp::T2 
    """Interpolation method"""
    method::T3 
end 

(interp::Interpolant)(x...) = interp.itp(x...)

"""
    $SIGNATURES

Returns an [`Interpolant`](@ref) by interpolating `pts` using `method`. `f0` is the initial function. `niter` is the number 
of iterations during the construction of the interpolant. `method` can be 

* [`Interp1D`](@ref) 
* [`HInterp1D`](@ref) 
* [`Interp2D`](@ref) 
* [`HInterp2D`](@ref) 

# Example 
```julia 
julia> xi, dx, xf = 0., 0.01, 10.;

julia> x = collect(xi : dx : xf);

julia> y = sin.(x);

julia> pts = collect.(zip(x, y)); 

julia> interp = integrate(pts, Interp1D(0.001));  # One dimensional interpolation. 
```
"""
interpolate(pts::AbstractVector{<:AbstractVector{<:Real}}, method::AbstractInterp; f0 = getinitf(method), niter::Int = 10) = 
    interpolate(map(pnt -> Point(pnt...), pts), method, f0=f0, niter=niter)

function interpolate(pts::PointVector, method::AbstractInterp; f0 = getinitf(method), niter::Int = 10) 
    tess = tessellate(pts, method)
    transforms = gettransforms(pts, method)
    mappings = getmappings(transforms, method)
    itp = wrap(f0, tess, mappings, niter)[1]
    Interpolant(IFS(transforms), itp, method)
end 


getinitf(::Interp1D)   = x -> 0. 
getinitf(::HInterp1D)  = x -> [0., 0.]
getinitf(::Interp2D)   = (x, y) -> 0. 
getinitf(::HInterp2D)  = (x, y) -> [0., 0.] 

"""
    $SIGNATURES

Projects `pts` on a lower dimensioanl space by dropping the `drop` number of indexes. 
"""
project(pts::PointVector, drop::Int=1)     = [Point(pnt[1 : end - drop]...) for pnt in pts]
project(pts::PointVector{2}, ::Interp1D)   = project(pts, 1)
project(pts::PointVector{3}, ::HInterp1D)  = project(pts, 2)
project(pts::PointVector{3}, ::Interp2D)   = project(pts, 1)
project(pts::PointVector{4}, ::HInterp2D)  = project(pts, 2)

# We use tessellation to locate the points. The fractal interpolant that is returned by `interpolate` function is a piecewise
# function. So to calculate the interpolant `interp` at a point `pnt`, we need to locate `pnt` to find the correct subdomain.
# We use tessellations to locate the points. See also `locate` function 

tessellate(pts::PointVector{2}, method::Interp1D)  = LineString(project(pts, method))
tessellate(pts::PointVector{3}, method::HInterp1D) = LineString(project(pts, method)) 
tessellate(pts::PointVector{3}, method::Interp2D)  = spt.Delaunay(project(pts, method))
tessellate(pts::PointVector{4}, method::HInterp2D) = spt.Delaunay(project(pts, method))

_locate(pnt::AbstractPoint{1, T}, tess::LineString) where {T} = findfirst(((p1, p2),) -> p1[1] ≤ pnt[1] ≤ p2[1], tess)
_locate(pnt::AbstractPoint{2, T}, tess::PyObject)   where {T} = tess.find_simplex(pnt)[1] + 1  
function locate(pnt::AbstractPoint, tess::Tessellation)
    n = _locate(pnt, tess)
    n == 0 || return n 
    pnt = doubleprecision(pnt)
    locate(pnt, tess) 
end 

doubleprecision(pnt::AbstractPoint{N, <:Real}) where {N} = Point(BigFloat.(pnt)...)
function doubleprecision(pnt::AbstractPoint{N, <:BigFloat}) where {N}
    prec = 2 * precision(BigFloat)
    prec ≤ MAXPREC ? setprecision(prec) : error("Exceeded maximum allowed precision $MAXPREC") 
    Point(BigFloat.(pnt)...)
end 

# IFS coefficients of the interpolant is found by using a linear algrebraic equation system using the boundary conditions.
# Each subtransformation in the transformations of a IFS maps a larger domain (Line in case of curve interpolation and
# Triangle in case of surface interpolation) to a smaller domain. These mappings maps the boundary points of the larger
# domains to boundary points of the smaller domains. `partition` function returns these smaller domains. 

partition(pts::PointVector, method::AbstractCurveInterp) = LineString(pts) 
function partition(pts::PointVector, method::AbstractSurfaceInterp, tess::Tessellation=tessellate(pts, method)) 
    trifaces = [TriangleFace(val...) for val in eachrow(tess.simplices .+ 1)]
    GeometryBasics.Mesh(pts, trifaces)
end 

# In case of curve interpolations, the boundary is the line that connects the endpoints of the interpolation domain. In case
# of surface interpolations, the interpolation domain is triangle, than that triangle is the boundary. If the interpolation
# domain is ngon (such as, tetragon, pentagon, hexagon, etc.) the triangle that can be drawn inside the convex hull of the
# points and that has the largest area is returned.

getboundary(pts::PointVector, ::AbstractCurveInterp) = Line(pts[1], pts[end])
function getboundary(pts::PointVector, method::AbstractSurfaceInterp) 
    hull = spt.ConvexHull(collect(hcat(collect.(project(pts, method))...)') )
    if length(hull.vertices) == 3   # Triangle  
        Triangle(Point.(pts[hull.vertices .+ 1])...)
    else    # Ngon (tetragon, pentagon, hexagon)
        polygon = Ngon(
            SVector{length(hull.vertices)}(Point.(pts[hull.vertices .+ 1]))
        )
        msh = GeometryBasics.mesh(polygon)
        idx = argmax([area(trig.points) for trig in msh])
        msh[idx]
    end 
end 

# `gettransforms` returns tranforms that maps outer domains to smaller domains in the interpolation domain. 
function gettransforms(pts::PointVector, method::AbstractInterp) 
    parts = partition(pts, method)
    n = length(parts)
    freevars = typeof(method.freevars) <: AbstractVector ? method.freevars : fill(method.freevars, n)
    boundary = getboundary(pts, method)
    map(((domain, freevar),) -> _gettransform(boundary, domain, freevar), zip(parts, freevars))
end 


# Interp1D
function _gettransform(outline::Line, inline::Line, freevar::Real) 
    outmat = collect(hcat(coordinates(outline)...)')
    inmat = collect(hcat(coordinates(inline)...)')
    inmat[:, end] -= outmat[:, end] * freevar
    outmat[:, end] .= 1
    a11, b1, a21, b2 = outmat \ inmat
    A = [a11    0; 
         a21    freevar]
    b = [b1, b2]
	Transformation(A, b)
end 

# HInterp1D
function _gettransform(outline::Line, inline::Line, freevar::AbstractMatrix) 
    outmat = collect(hcat(coordinates(outline)...)')
    inmat = collect(hcat(coordinates(inline)...)')
    inmat[:, 2 : 3] -= outmat[:, 2 : 3] * freevar'
    outmat = [outmat[:, 1] ones(2)]
    a11, b1, a21, b2, a31, b3 = outmat \ inmat
    A = [a11    0               0; 
         a21    freevar[1,1]    freevar[1,2];
         a31    freevar[2,1]    freevar[2,2]
         ]
    b = [b1, b2, b3]
	Transformation(A, b)
end 

# Interp2D
function _gettransform(outtrig::Triangle, intrig::Triangle, freevar::Real) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, end] -= outmat[:, end] * freevar
    outmat[:, end] .= 1
    a11, a12, b1, a21, a22, b2, a31, a32, b3 = outmat \ inmat
    A = [a11    a12    0; 
         a21    a22    0;
         a31    a32    freevar
         ]
    b = [b1, b2, b3]
	Transformation(A, b)
end 

# HInterp2D
function _gettransform(outtrig::Triangle, intrig::Triangle, freevar::AbstractMatrix) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, 3:4] -= outmat[:, 3:4] * freevar'
    outmat = [outmat[:, 1:2] ones(3)]
    a11, a12, b1, a21, a22, b2, a31, a32, b3, a41, a42, b4 = outmat \ inmat
    A = [a11    a12    0                0; 
         a21    a22    0                0;
         a31    a32    freevar[1, 1]    freevar[1, 2];
         a41    a42    freevar[2, 1]    freevar[2, 2];
         ]
    b = [b1, b2, b3, b4]
	Transformation(A, b)
end 

getmappings(transforms, method) = map(transform -> _getmapping(transform, method), transforms)

function _getmapping(transform, method::Interp1D)
    (a11, a21, _, a22), (b1, b2) = transform.A, transform.b
    linv = x -> (x - b1) / a11 
    F = (x, y) -> a21 * x + a22 * y + b2
    (linv, F)
end 

function _getmapping(transform, method::HInterp1D)
    (a11, a21, a31, _, a22, a32, _, a23, a33), (b1, b2, b3) = transform.A, transform.b
    linv = x -> (x - b1) / a11 
    F = (x, y, z) ->  [a21 a22 a23; a31 a32 a33] * [x, y, z] + [b2, b3]
    (linv, F)
end 

function _getmapping(transform, method::Interp2D)
    (a11, a21, a31, a12, a22, a32, _, _, a33), (b1, b2, b3) = transform.A, transform.b
    linv = (x, y) -> [a11 a12; a21 a22] \ ([x, y] - [b1, b2])
    F = (x, y, z) ->  a31 * x + a32 * y + a33 * z + b3
    (linv, F)
end 

function _getmapping(transform, method::HInterp2D)
    (a11, a21, a31, a41, a12, a22, a32, a42, _, _, a33, a34, _, _, a43, a44), (b1, b2, b3, b4) = transform.A, transform.b
    linv = (x, y) -> [a11 a12; a21 a22] \ ([x, y] - [b1, b2])
    F = (x, y, z, t) ->  [a31 a32 a33 a34; a41 a42 a43 a44] * [x, y, z, t] + [b3, b4]
    (linv, F)
end 

# Main iteration of the initial function `f0` for `niter` iterations. 
wrap(f0, tess::Tessellation, mappings::AbstractVector{<:Tuple{T, S}}, niter::Int) where {T, S} = 
    ((f0, tess, mappings)) |> ∘((wrapper for i in 1 : niter)...) 

function wrapper((f, tess, mappings))
    function fnext(x...) 
        pnt = Point(x...) 
        n = locate(pnt, tess)
        # n == 0 && error("Point $pnt cannot be found.")
        linv, F = mappings[n]
        val = linv(x...) 
        F(val..., f(val...)...)
    end, tess, mappings
end
