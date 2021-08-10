# This file includes methods for data generation. 

export Dataset, Tessellation, getdata, boundarypoints, interiorpoints, getpoint, project, uniformdomain, tomesh, locate

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

function gettriangles(tess::Tessellation)
    T = typeof(triangulationtype(tess))
    if T == ScipyTriangulation 
        tess.tess.simplices .+ 1
    elseif T == MatplotlibTriangulation
        tess.tess.triangles .+ 1
    else 
        collect(1 : length(tess.tess))
    end 
end 

function locate(tess::Tessellation, point::AbstractVector)
    T = typeof(triangulationtype(tess))
    if T == ScipyTriangulation
        n = tess.tess.find_simplex(point)[1] + 1  
    elseif T == MatplotlibTriangulation  
        trf = tess.tess.get_trifinder()
        n = trf([point[1]], [point[2]])[1] + 1 
    else
        n = findfirst(((p1, p2),) -> p1[1] ≤ point[1] ≤ p2[1], tess.tess)
    end 
    n == 0 || return n 
    point = doubleprecision(point)
    locate(tess, point) 
end 


function doubleprecision(point)
    if eltype(point) == BigFloat
        prec = 2 * precision(BigFloat)
        if prec ≤ MAXPREC
            @info "BigFloat precision: $prec"
            setprecision(prec)
        else 
            error("Exceeded maximum allowed precision $MAXPREC for point $point") 
        end 
        BigFloat.(point)
    else
        prec = precision(BigFloat)
        @info "Passing to arbitrary precision arithmetic. BigFloat precision: $prec"
        BigFloat.(string.(point))
    end
end 


# ---------------------------------------- Dataset -------------------------- # 

struct Dataset{T1<:AbstractVector, T2<:AbstractVector}
    points::T1 
    domain::T2 
end 

Dataset(points::AbstractVector) = Dataset(points, [])

Dataset(domain::AbstractVector, ninterpoints::Int, nedgepoints::Int) = 
    Dataset(getdata(domain, ninterpoints, nedgepoints), domain)

Dataset(f, domain::AbstractVector, ninterpoints::Int, nedgepoints::Int) = 
    Dataset(getdata(f, domain, ninterpoints, nedgepoints), domain)

getdata(f, domain::AbstractVector, ninterpoints::Int, nedgepoints::Int) = 
    map(p -> [p; f(p...)], getdata(domain, ninterpoints, nedgepoints))

function getdata(domain::AbstractVector, ninterpoints::Int, nedgepoints::Int)
    n = length(domain[1]) 
    pdomain = project(domain, n - 2)
    if length(pdomain) == 2
        boundarypoints(pdomain, nedgepoints)
    else
        interpoints = interiorpoints(pdomain, ninterpoints)
        boundpoints = boundarypoints(pdomain, nedgepoints)
        # `domain` is inserted to the beginning of the points. this is very place of `domian` is very important because the 
        # constrains of the triangulation is performed with respect to the position of the domain. 
        vcat(pdomain, boundpoints, interpoints) |> unique
    end
end 


interiorpoints(f, domain::AbstractVector, ninterpoints::Int=100) = 
    map(p -> [p; f(p...)], interiorpoints(domain, ninterpoints))

function interiorpoints(domain::AbstractVector, ninterpoints::Int=100)
    Meshes.sample(
        Meshes.Ngon(Meshes.Point.(domain)), 
        Meshes.HomogeneousSampling(ninterpoints)
        ) |> collect .|> Meshes.coordinates .|> collect
end 

boundarypoints(f, domain::AbstractVector, nedgepoints::Int=10) = map(p -> [p; f(p...)], boundarypoints(domain, nedgepoints))

function boundarypoints(domain::AbstractVector, nedgepoints::Int=10)
    if length(domain) == 2 # domain is a line 
        pts = edgepoints(domain[1], domain[2], nedgepoints)
    else    # domain is an ngon 
        pts = vcat(map(((pnt1, pnt2),) -> edgepoints(pnt1, pnt2, nedgepoints), 
            TupleView{2, 1}(SVector([domain; [domain[1]]]...)))...)
    end
    unique(pts)
end

function edgepoints(p1::AbstractVector, p2::AbstractVector, nedgepoints::Int=10) 
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

uniformdomain(n::Int, T::Type=Float64, p0::AbstractVector=zeros(T, 2), r::Real=1.) = 
    [p0 + T[r * cos(θ), r*sin(θ)] for θ in (0 : n - 1) / n * 2π]
uniformdomain(f, n::Int, T::Type=Float64, p0::AbstractVector=zeros(T, 2), r::Real=1.) = 
    map(point -> [point; f(point...)], uniformdomain(n, T, p0, r))

Tessellation(dataset::Dataset) = Tessellation(dataset.points, dataset.domain)

ngon(pts) = Ngon(SVector{length(pts)}(map(item -> Point(item...), pts)))

function getpoint(domain) 
    if length(domain) == 2
        n = length(domain[1])
        pdomain = project(domain, n - 2) 
        p0, p1 = pdomain[1], pdomain[2]
        p0 + (p1 - p0) * rand(T)
    else
        interiorpoints(domain, 1)
    end
end 

# ---------------------------------- Meshing --------------------------------------- # 

tomesh(dataset::Dataset, tess=nothing) = tomesh(dataset.points, dataset.domain, tess)

function tomesh(points::AbstractVector, domain::AbstractVector=[], tess=nothing)
    if length(points[1]) == 1 
        return LineString(points)
    end 
    if tess === nothing
        tess = project_and_tessellate(points, domain)
    end 
    faces = [TriangleFace(idx[1], idx[2], idx[3]) for idx in eachrow(gettriangles(tess))]
    return GeometryBasics.Mesh([Point(point...) for point in points], faces)
end

function project_and_tessellate(points::AbstractVector, domain::AbstractVector=[])
    n = length(points[1]) 
    if 3 ≤ n ≤ 4  
        _points = project(points, n - 2)
    else
        error("Expected dimension of points is 3 or 4, got $n instead.")
    end 
    Tessellation(_points, domain)
end 

# --------------------------------------- Projection ------------------------------------------------------# 

project(points::AbstractVector, n=1) = _project.(points, n)
project(dataset::Dataset, n=1) = Dataset(project(dataset.points, n), dataset.domain)
project(msh::GeometryBasics.Mesh, n=1) = GeometryBasics.Mesh(project(msh.position, n), faces(msh)) 

_project(point::AbstractVector, n=1) = point[1 : end - n]
_project(point::AbstractPoint{Dim, T}, n=1) where {Dim, T} = Point{Dim - n, T}(point[1 : end - n]...)
