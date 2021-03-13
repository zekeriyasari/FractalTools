# This includes interpolation methods 

export Interp1D, Interp2D, HInterp1D, HInterp2D, interpolate, project

const PointVector{Dim} = AbstractVector{<:AbstractPoint{Dim, T}} where {T}
const Tessellation = Union{<:LineString, <:PyObject}

abstract type AbstractInterp end
abstract type AbstractCurveInterp   <: AbstractInterp end
abstract type AbstractSurfaceInterp <: AbstractInterp end

struct Interp1D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: AbstractCurveInterp
    freevars::T 
end 

struct Interp2D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: AbstractSurfaceInterp
    freevars::T
end 

struct HInterp1D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: AbstractCurveInterp
    freevars::T 
end 

struct HInterp2D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: AbstractSurfaceInterp
    freevars::T 
end 

# Note: For plotting recipes, Interpolant should a subtype of Function
struct Interpolant{T1<:IFS, T2, T3<:AbstractInterp} <: Function 
    ifs::T1 
    itp::T2 
    method::T3 
end 

(interp::Interpolant)(x...) = interp.itp(x...)

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


project(pts::PointVector, drop::Int=1) = [Point(pnt[1 : end - drop]...) for pnt in pts]
project(pts::PointVector{2}, ::Interp1D)   = project(pts, 1)
project(pts::PointVector{3}, ::HInterp1D)  = project(pts, 2)
project(pts::PointVector{3}, ::Interp2D)   = project(pts, 1)
project(pts::PointVector{4}, ::HInterp2D)  = project(pts, 2)


tessellate(pts::PointVector{2}, method::Interp1D)  = LineString(project(pts, method))
tessellate(pts::PointVector{3}, method::HInterp1D) = LineString(project(pts, method)) 
tessellate(pts::PointVector{3}, method::Interp2D)  = spt.Delaunay(project(pts, method))
tessellate(pts::PointVector{4}, method::HInterp2D) = spt.Delaunay(project(pts, method))


partition(pts::PointVector, method::AbstractCurveInterp) = LineString(pts) 
function partition(pts::PointVector, method::AbstractSurfaceInterp, tess::Tessellation=tessellate(pts, method)) 
    trifaces = [TriangleFace(val...) for val in eachrow(tess.simplices .+ 1)]
    GeometryBasics.Mesh(pts, trifaces)
end 


getboundary(pts::PointVector, ::AbstractCurveInterp) = Line(pts[1], pts[end])
function getboundary(pts::PointVector, method::AbstractSurfaceInterp) 
    hull = spt.ConvexHull(collect(hcat(collect.(project(pts, method))...)') )
    if length(hull.vertices) == 3
        # If the convex hull is a triangle, just return the triangle 
        Triangle(Point.(pts[hull.vertices .+ 1])...)
    else 
        # If the the convex hull is not a triangle but a polygon, construct the 
        # boundary polygon, construct a mesh from the polygon and return the 
        # maximum triangle with the maximum area.
        polygon = Ngon(
            SVector{length(hull.vertices)}(Point.(pts[hull.vertices .+ 1]))
        )
        msh = GeometryBasics.mesh(polygon)
        idx = argmax([area(trig.points) for trig in msh])
        msh[idx]
    end 
end 


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

locate(pnt::AbstractPoint{1, T}, tess::LineString) where {T} = findfirst(((p1, p2),) -> p1[1] ≤ pnt[1] ≤ p2[1], tess)
locate(pnt::AbstractPoint{2, T}, tess::PyObject)   where {T} = tess.find_simplex(pnt)[1] + 1  

wrap(f0, tess::Tessellation, mappings::AbstractVector{<:Tuple{T, S}}, niter::Int) where {T, S} = 
    ((f0, tess, mappings)) |> ∘((wrapper for i in 1 : niter)...) 

function wrapper((f, tess, mappings))
    function fnext(x...) 
        pnt = Point(x...) 
        n = locate(pnt, tess)
        n == 0 && error("Point $pnt cannot be found.")
        linv, F = mappings[n]
        val = linv(x...) 
        F(val..., f(val...)...)
    end, tess, mappings
end
