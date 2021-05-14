# This file includes methods for data generation. 

export getdata, getpoint

"""
    getdata(f, vtx, npts) 

Returns a vector of (d + m)-dimensional points that are dispersed randomly over d-dimensioanl domain with vertex points
`vtx`. `f` is the vector-valued function used to map each point in the domain. `npts` is the number of points to returned.

    getdata(vtx, npts) 

Returns a vector of d-dimensional points that are dispersed randomly over d-dimensioanl domain with vertex points
`vtx`. `npts` is the number of points to be returned` 
"""
getdata(f, vtx::AbstractVector{<:AbstractVector}, npts::Int) = getdata(f, ngon(vtx), npts)
getdata(vtx::AbstractVector{<:AbstractVector}, npts::Int)    = getdata(ngon(vtx), npts)
getdata(f, domain::Ngon, npts::Int) = [Point(pnt..., f(pnt...)...) for pnt in disperse(domain, npts)] 
getdata(domain::Ngon, npts::Int)    = disperse(domain, npts)

ngon(pts::AbstractVector{<:AbstractVector}) = Ngon(SVector{length(pts)}(map(item -> Point(item...), pts)))

# TODO: Test `getdata` function with `Line` domains.
# FIXME: #37 Test `getdata` function returns correct number of points.

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
    return Point(NaN, NaN)  # For type-stability 
end 

function getpoint(line::Line{1, T}) where {T}
    p0, p1 = only(line[1]), only(line[2])
    pnt = p0 + (p1 - p0) * rand(T)
    Point(pnt)
end 

boundarypoints(line::Line; numpoints::Int=10) = linepoint(line[1], line[2])
boundarypoints(ngon::Ngon; numpoints::Int=10) = vcat(
        map(
            ((pnt1, pnt2),) -> linepoint(pnt1, pnt2), TupleView{2, 1}(SVector([ngon.points; [ngon.points[1]]]...))
            )...
        )

function linepoint(p1::AbstractPoint, p2::AbstractPoint; numpoints::Int=10) 
    x = collect(LinRange(p1[1], p2[1], numpoints))
    y = collect(LinRange(p1[2], p2[2], numpoints))
    [Point(xi, yi) for (xi, yi) in zip(x, y)]
end

isvalidpoint(pnt::AbstractPoint, tess) = tess.find_simplex(pnt)[1] ≥ 0 && pnt !== Point(NaN, NaN)

function disperse(ngon::Ngon, npoints::Int) 
    allpnts = [getpoint(ngon) for i in 1 : 10 * npoints]
    ctrpnts = [Point(val[1], val[2]) for val in eachcol(kmeans(hcat(collect.(allpnts)...), npoints).centers)]
    vcat(ngon.points, boundarypoints(ngon), ctrpnts)
end

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
