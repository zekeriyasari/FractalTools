# This includes interpolation methods 

export Interp1D, Interp2D, HInterp1D, HInterp2D, interpolate, tessellate, project

const PointVector{Dim} = AbstractVector{<:AbstractPoint{Dim, T}} where {T}

# ------------------------------------ Interpolation types ------------------------------------ # 

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
        a_{11, n}   & a_{11, n}     & 0             & 0         \\
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


# ------------------------------------ Interpolation methods------------------------------------ # 

interpolate(points::AbstractVector, domain::AbstractVector, method::AbstractInterp; kwargs...) = 
    interpolate(Dataset(points, domain), method; kwargs...)

function interpolate(dataset::Dataset, method::AbstractInterp; f0 = getinitf(method), niter::Int = 10) 
    tess = tessellate(dataset, method)
    transforms = gettransforms(dataset, method, tess)
    mappings = getmappings(transforms, method)
    itp = wrap(f0, tess, mappings, niter)[1]
    Interpolant(IFS(transforms), itp, method)
end 

getinitf(::Interp1D)   = x -> 0. 
getinitf(::HInterp1D)  = x -> [0., 0.]
getinitf(::Interp2D)   = (x, y) -> 0. 
getinitf(::HInterp2D)  = (x, y) -> [0., 0.] 

project(points::AbstractVector, ::Interp1D)  = project(points, 1)
project(points::AbstractVector, ::HInterp1D) = project(points, 2)
project(points::AbstractVector, ::Interp2D)  = project(points, 1)
project(points::AbstractVector, ::HInterp2D) = project(points, 2)

# We use tessellation to locate the points. The fractal interpolant that is returned by `interpolate` function is a piecewise
# function. So to calculate the interpolant `interp` at a point `pnt`, we need to locate `pnt` to find the correct subdomain.
# We use tessellations to locate the points. See also `locate` function 

tessellate(dataset::Dataset, method::AbstractInterp) = Tessellation(project(dataset.points, method), dataset.domain)

# In case of curve interpolations, the boundary is the line that connects the endpoints of the interpolation domain. In case
# of surface interpolations, the interpolation domain is triangle, than that triangle is the boundary. If the interpolation
# domain is ngon (such as, tetragon, pentagon, hexagon, etc.) the triangle that can be drawn inside the convex hull of the
# points and that has the largest area is returned.

getboundary(dataset::Dataset, ::AbstractCurveInterp, ::Tessellation) = 
    GeometryBasics.Line(Point(dataset.points[1]...), Point(dataset.points[end]...))

function getboundary(dataset::Dataset, ::AbstractSurfaceInterp, tess::Tessellation) 
    domain = dataset.domain
    if length(domain) == 3  # domain is already triangle 
        Triangle([Point(point...) for point in domain]...)
    else
        #= Note
        For the case, hidden surface fractal interpolation, the points `points` are in four-dimensional space.  However, 
        GeometryBasics.mesh requires at most three-dimensional points. Thus, we cannot use the method `GeometryBasics.mesh` 
        to constuct the mesh for the convex hull point. Instead, we first find the hull points, then construct d-dimensional 
        mesh (where d = 3 for surface interpolation and d = 4 for hidden surface interpolation. Then, we find index of the 
        triangle with the largest area  in the interpolation domain (which is to 2-dimensional domain) and we return the 
        triangle in the d-dimensioanl mesh  
        =#
        msh = tomesh(domain)
        n = length(dataset.domain[1])
        idx = argmax([area(trig.points) for trig in project(msh, n - 2)]) 
        msh[idx]
    end 
end 

# `gettransforms` returns tranforms that maps outer domains to smaller domains in the interpolation domain. 
function gettransforms(dataset::Dataset, method::AbstractInterp, tess::Tessellation) 
    # IFS coefficients of the interpolant is found by using a linear algrebraic equation system using the boundary 
    # conditions.Each subtransformation in the transformations of a IFS maps a larger domain (Line in case of curve 
    # interpolation and Triangle in case of surface interpolation) to a smaller domain. These mappings maps the boundary 
    # points of the larger domains to boundary points of the smaller domains. `partition` function returns these smaller 
    # domains. 
    parts = tomesh(dataset, tess)
    n = length(parts)
    freevars = typeof(method.freevars) <: AbstractVector ? method.freevars : fill(method.freevars, n)
    boundary = getboundary(dataset, method, tess)
    map(((domain, freevar),) -> _gettransform(boundary, domain, freevar), zip(parts, freevars))
end 


# Interp1D
function _gettransform(outline::GeometryBasics.Line, inline::GeometryBasics.Line, freevar::Real) 
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
function _gettransform(outline::GeometryBasics.Line, inline::GeometryBasics.Line, freevar::AbstractMatrix) 
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

# Here transforms and mappings are different. The type of a transform is FractalTools.Transformation. A
# FractalTools.Transformation object basically consists of the matrix A and vector b of an affine transformation w(x) = A * x
# + b. However, a mapping is a subtransformation (defined as Ln: Ω ↦ Ωₙ and Fn: Ω × Rᵐ ↦ Rᵖ, where Ω is the interpolation
# domain and Ωₙ is the subdomain obtained by partitioning the interpolation domian) of an affine transformation in an IFS.
# So, we compute the mappings (i.e. subtransformations) from transforms and use directly the mappings to compute the
# interpolant. 

getmappings(transforms, method) = map(transform -> _getmapping(transform, method), transforms)

function _getmapping(transform, ::Interp1D)
    (a11, a21, _, a22), (b1, b2) = transform.A, transform.b
    linv = x -> (x - b1) / a11 
    F = (x, y) -> a21 * x + a22 * y + b2
    (linv, F)
end 

function _getmapping(transform, ::HInterp1D)
    (a11, a21, a31, _, a22, a32, _, a23, a33), (b1, b2, b3) = transform.A, transform.b
    linv = x -> (x - b1) / a11 
    F = (x, y, z) ->  [a21 a22 a23; a31 a32 a33] * [x, y, z] + [b2, b3]
    (linv, F)
end 

function _getmapping(transform, ::Interp2D)
    (a11, a21, a31, a12, a22, a32, _, _, a33), (b1, b2, b3) = transform.A, transform.b
    linv = (x, y) -> [a11 a12; a21 a22] \ ([x, y] - [b1, b2])
    F = (x, y, z) ->  a31 * x + a32 * y + a33 * z + b3
    (linv, F)
end 

function _getmapping(transform, ::HInterp2D)
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
        point = [x...] 
        # Because of finite precison arithmetic, the location of a valid point may fail. To overcome this, we can switch to
        # arbitrary # precision arithmetic (by using BigFloat) in the expense of increasing the computationonal complexity.
        # However, not to # increase the computationonal complexity directly, the strategy employed here is this:  before
        # switching completely to # arbitrary precision arithmetic, we first try to locate the point with finite precision
        # (Float64 precision). If the point # cannot be located, we double the precision and try to locate the point again.
        # We do this until maximum allowed # precison(specified to be 1024 bits as MAXPREC). If MAXPREC is reached and the
        # point still cannot be found, the point # location fails with an error message.   
        n = locate(tess, point)
        # n == 0 && error("Point $point cannot be found.")
        linv, F = mappings[n]
        val = linv(x...) 
        F(val..., f(val...)...)
    end, tess, mappings
end

