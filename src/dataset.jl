# This file includes methods for data generation. 

export Dataset, Tessellation, getdata, getpoint, project, uniformdomain, tomesh

# --------------------------------------------- Tessellation --------------------------------- # 

struct Tessellation{T} 
    tess::T 
end 

function Tessellation(points::AbstractVector, domain::AbstractVector=[], verbose::Bool=false) 
    if length(points[1]) == 1 
        tess = LineString([Point(point...) for point in points])
    elseif isempty(domain) 
        tess = spt.Delaunay(points)
    else 
        triin = Triangulate.TriangulateIO()
        triin.pointlist = hcat(points...)
        list = 1 : length(domain)
        triin.segmentmarkerlist = collect(list) 
        triin.segmentlist = hcat(collect.(zip(list, circshift(list, -1)))...)
        args = verbose ? "p" : "pQ"
        triout, vorout = Triangulate.triangulate(args, triin)
        tess = tri.Triangulation(triout.pointlist[1, :], triout.pointlist[2, :], triout.trianglelist' .- 1)
    end
    Tessellation{typeof(tess)}(tess)
end 

struct ScipyTriangulation end 
struct MatplotlibTriangulation end 
struct GeometryBasicsTriangulation end 

function triangulationtype(tess::Tessellation) 
    _tess = tess.tess 
    if typeof(_tess) <: LineString
        return GeometryBasicsTriangulation() 
    elseif _tess.__class__ == tri.triangulation.Triangulation 
        return MatplotlibTriangulation() 
    elseif _tess.__class__ == spt.qhull.Delaunay
        return ScipyTriangulation() 
    else
        error("Unknown triangulation type") 
    end 
end 


# ---------------------------------------- Dataset -------------------------- # 

struct Dataset{T1<:AbstractVector, T2<:AbstractVector}
    points::T1 
    domain::T2 
end 

Dataset(domain::AbstractVector, ninterpoints::Int, nedgepoints::Int) = 
    Dataset(getdata(domain, ninterpoints, nedgepoints), domain)

Dataset(f, domain::AbstractVector, ninterpoints::Int, nedgepoints::Int) = 
    Dataset(getdata(f, domain, ninterpoints, nedgepoints), domain)

function getdata(domain, ninterpoints, nedgepoints)
    if length(domain) == 2
        boundarypoints(domain, nedgepoints)
    else
        interpoints = interiorpoints(domain, ninterpoints)
        boundpoints = boundarypoints(domain, nedgepoints)
        # `domain` is inserted to the beginning of the points. this is very place of `domian` is very important because the 
        # constrains of the triangulation is performed with respect to the position of the domain. 
        vcat(domain, boundpoints, interpoints) |> unique
    end
end 
getdata(f, domain, ninterpoints, nedgepoints) = map(p -> [p; f(p...)], getdata(domain, ninterpoints, nedgepoints))

function interiorpoints(domain, ninterpoints=100)
    Meshes.sample(
        Meshes.Ngon(Meshes.Point.(domain)), 
        Meshes.HomogeneousSampling(ninterpoints)
        ) |> collect .|> Meshes.coordinates .|> collect
end 

function boundarypoints(domain, nedgepoints=10)
    if length(domain) == 2 # domain is a line 
        edgepoints(domain[1], domain[2], nedgepoints)
    else    # domain is an ngon 
        vcat(map(((pnt1, pnt2),) -> edgepoints(pnt1, pnt2, nedgepoints), 
            TupleView{2, 1}(SVector([domain; [domain[1]]]...)))...)
    end
end

function edgepoints(p1, p2, nedgepoints=10) 
    n = length(p1) 
    m = length(p2) 
    if m == n == 1   # p1 and p2 are one dimensional space
        [[xi] for xi in collect(LinRange(p1[1], p2[1], nedgepoints))]
    else    # p1 and p2 are one n-dimensional space
        x = collect(LinRange(p1[1], p2[1], nedgepoints))
        y = collect(LinRange(p1[2], p2[2], nedgepoints))
        [[xi, yi] for (xi, yi) in zip(x, y)]
    end 
end

uniformdomain(n, T=Float64, p0=zeros(T, 2), r=1.) = [p0 + T[r * cos(θ), r*sin(θ)] for θ in (0 : n - 1) / n * 2π]

Tessellation(dataset::Dataset) = Tessellation(dataset.points, dataset.domain)

ngon(pts) = Ngon(SVector{length(pts)}(map(item -> Point(item...), pts)))

function getpoint(domain) 
    if length(domain) == 2 
        p0, p1 = domain[1], line[2]
        p0 + (p1 - p0) * rand(T)
    else
        interiorpoints(domain, 1)
    end
end 

# ---------------------------------- Meshing --------------------------------------- # 

tomesh(dataset::Dataset) = tomesh(dataset.points, dataset.domain)
function tomesh(points::AbstractVector, domain=[])
    n = length(points[1]) 
    if n > 2 
        _points = project(points, n - 2)
    end 
    if isempty(domain)
        tess = Tessellation(domain, _points)
        faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.triangles .+ 1)]
    else 
        tess = spt.Delaunay(_points) 
        faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    end 
    GeometryBasics.Mesh([Point(point...) for point in points], faces)
end

# --------------------------------------- Projection ------------------------------------------------------# 

project(points::AbstractVector, n=1) = _project.(points, n)
project(dataset::Dataset, n=1) = Dataset(project(dataset.points, n), dataset.domain)
project(msh::GeometryBasics.Mesh, n=1) = GeometryBasics.Mesh(project(msh.position, n), faces(msh)) 

_project(point::AbstractVector, n=1) = point[1 : end - n]
_project(point::AbstractPoint{Dim, T}, n=1) where {Dim, T} = Point{Dim - n, T}(point[1 : end - n]...)
