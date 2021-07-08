# This file includes methods for data generation. 

export getdata, getpoint, uniformdomain, tomesh

"""
    getdata(f, vtx, npts) 

Returns a vector of (d + m)-dimensional points that are dispersed randomly over d-dimensioanl domain with vertex points
`vtx`. `f` is the vector-valued function used to map each point in the domain. `npts` is the number of points to returned.

    getdata(vtx, npts) 

Returns a vector of d-dimensional points that are dispersed randomly over d-dimensioanl domain with vertex points
`vtx`. `npts` is the number of points to be returned` 
"""
function getdata end 

getdata(f, vtx::AbstractVector{<:AbstractVector}, args::Int...) = getdata(f, ngon(vtx), args...)
getdata(vtx::AbstractVector{<:AbstractVector}, args::Int...) = getdata(ngon(vtx), args...)

getdata(domain::Line, nedgepoints::Int) = boundarypoints(domain, nedgepoints)
getdata(f, domain::Line, nedgepoints::Int) = [Point(pnt..., f(pnt...)...) for pnt in boundarypoints(domain, nedgepoints)] 

# nedgepoints=10 is for backwards compatibility
getdata(domain::Ngon, ninterpoints::Int, nedgepoints::Int=10) = disperse(domain, ninterpoints, nedgepoints)
getdata(f, domain::Ngon, ninterpoints::Int, nedgepoints::Int=10) = 
    [Point(pnt..., f(pnt...)...) for pnt in disperse(domain, ninterpoints, nedgepoints)] 

ngon(pts::AbstractVector{<:AbstractVector}) = Ngon(SVector{length(pts)}(map(item -> Point(item...), pts)))

"""
    getpoint(ngon::Ngon, maxiter::Int=100_000) 

Returns a valid that is inside of `ngon`.

    getpoint(line::Line) 

Returns a valid point that is inside `line`. 
"""
function getpoint end 

function getpoint(ngon::Ngon{Dim, T, N, P}; maxiter::Int=100_000) where {Dim, T, N, P}
    tess = spt.Delaunay(coordinates(ngon))      
    A, b = boundboxtransforms(ngon) 
    iter = 1 
    while iter ≤ maxiter 
        val = A * rand(T, Dim) + b 
        pnt = Point(val...) 
        isvalidpoint(pnt, tess) && return pnt
        iter += 1
    end 
    return Point(NaN, NaN) 
end 

function getpoint(line::Line{1, T}) where {T}
    p0, p1 = only(line[1]), only(line[2])
    pnt = p0 + (p1 - p0) * rand(T)
    Point(pnt)
end 

"""
    $SIGNATURES

Returns a uniform ngon whose vertex points are centered at `p0` with a radius `r`.
"""
uniformdomain(n::Int, T::Type{<:Real}=Float64, p0::AbstractVector{<:Real}=zeros(T, 2), r::Real=1.) =
    [p0 + T[r * cos(θ), r*sin(θ)] for θ in (0 : n - 1) / n * 2π] |> ngon 

"""
    $SIGNATURES

Returns a mesh whose points are `pnts`
"""
tomesh(pnts::AbstractVector) = GeometryBasics.Mesh(topoint.(pnts), tofaces(topoint.(pnts)))

tofaces(pnts::AbstractVector{<:AbstractPoint{3,T}}) where {T} = tofaces(project(pnts))
function tofaces(pnts::AbstractVector{<:AbstractPoint{2,T}}) where {T}
    tess = spt.Delaunay(pnts) 
    [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
end 

topoint(pnt::AbstractPoint) = pnt
topoint(vect::AbstractVector{<:Real}) = Point(vect...)
topoint(vect::Real) = Point(vect)

tovector(pnt::AbstractPoint) = [pnt...]
tovector(vect::AbstractVector{<:Real}) = vect


isvalidpoint(pnt::AbstractPoint, tess) = tess.find_simplex(pnt)[1] ≥ 0 && pnt !== Point(NaN, NaN)


function disperse(ngon::Ngon, ninterpoints::Int, nedgepoints::Int) 
    interpoints = interiorpoints(ngon, ninterpoints)
    boundpoints = boundarypoints(ngon, nedgepoints)
    vcat(boundpoints, interpoints) |> unique
end

function interiorpoints(ngon::Ngon, ninterpoints::Int)
    allpnts = [getpoint(ngon) for i in 1 : 10 * ninterpoints]
    [Point(val[1], val[2]) for val in eachcol(kmeans(hcat(collect.(allpnts)...), ninterpoints).centers)]
end 

boundarypoints(line::Line, nedgepoints::Int=10) = edgepoints(line[1], line[2], nedgepoints)
boundarypoints(ngon::Ngon, nedgepoints::Int=10) = vcat(map(((pnt1, pnt2),) -> edgepoints(pnt1, pnt2, nedgepoints), 
    TupleView{2, 1}(SVector([ngon.points; [ngon.points[1]]]...)))...)

function edgepoints(p1::AbstractPoint{2, T}, p2::AbstractPoint{2, T}, nedgepoints::Int=10) where {T}
    x = collect(LinRange(p1[1], p2[1], nedgepoints))
    y = collect(LinRange(p1[2], p2[2], nedgepoints))
    [Point(xi, yi) for (xi, yi) in zip(x, y)]
end
edgepoints(p1::AbstractPoint{1, T}, p2::AbstractPoint{1, T}, nedgepoints::Int=10) where {T} = 
    [Point(xi) for xi in collect(LinRange(p1[1], p2[1], nedgepoints))]

function boundboxtransforms(ngon::Ngon)
    coords = coordinates(ngon) 
    x, y = getindex.(coords, 1), getindex.(coords, 2) 
    xmin, xmax, ymin, ymax = minimum(x), maximum(x), minimum(y), maximum(y)
    xwidth, ywidth = xmax - xmin, ymax - ymin 
    A = [xwidth 0; 0 ywidth]
    b = [xmin, ymin] 
    A, b
end

# TODO: Complete function
function filtertriangle(quality, args...; kwargs...) end 
