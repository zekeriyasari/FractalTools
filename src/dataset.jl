# This file includes methods for data generation. 

export Dataset, Tessellation, getdata, boundarypoints, interiorpoints, project, uniformdomain, tomesh, locate, 
    bigfloat, combine, generate, getpoints

# Methods arbitrary precision arithmetic
Meshes.atol(::Type{BigFloat}) = 100eps()
bigfloat(Ω::AbstractArray) = @. map(item -> BigFloat(string(item)), Ω)


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

function getpoints(tess::Tessellation)
    T = typeof(triangulationtype(tess))
    if T == ScipyTriangulation 
        tess.tess.points
    elseif T == MatplotlibTriangulation
        combine(tess.tess.x, tess.tess.y)
    else 
        coordinates(tess.tess) .|> collect
    end 
end 

function locate(tess::Tessellation, point::AbstractVector)
    # Construct a knn tree 
    pts = getpoints(tess) 
    tree = KDTree(hcat(pts...))

    # Search for the triangle 
    T = typeof(triangulationtype(tess))
    if T == ScipyTriangulation
        n = tess.tess.find_simplex(point)[1] + 1  
    elseif T == MatplotlibTriangulation  
        trf = tess.tess.get_trifinder()
        n = trf([point[1]], [point[2]])[1] + 1 
    else
        n = findfirst(((p1, p2),) -> p1[1] ≤ point[1] ≤ p2[1], tess.tess)
    end 

    # If point could not be located, locate its nearest neighbour
    if n == 0   
        point = pts[nn(tree, point) |> first]
        return locate(tess, point)
    end 
    return n
end 


# function doubleprecision(point)
#     if eltype(point) == BigFloat
#         prec = 2 * precision(BigFloat)
#         if prec ≤ MAXPREC
#             @info "BigFloat precision: $prec"
#             setprecision(prec)
#         else 
#             error("Exceeded maximum allowed precision $MAXPREC for point $point") 
#         end 
#         BigFloat.(point)
#     else
#         prec = precision(BigFloat)
#         @info "Passing to arbitrary precision arithmetic. BigFloat precision: $prec"
#         BigFloat.(string.(point))
#     end
# end 


# ---------------------------------------- Dataset -------------------------- # 

struct Dataset{T1<:AbstractVector, T2<:AbstractVector}
    points::T1 
    domain::T2 
end 

Dataset(points::AbstractVector) = Dataset(points, [])

Dataset(domain::AbstractVector, lc::Real) = Dataset(getdata(domain, lc)...)

Dataset(f, domain::AbstractVector, lc::Real=0.1; kwargs...) = Dataset(getdata(f, domain, lc; kwargs...)...)

function getdata(f, domain::AbstractVector, lc::Real=0.1; kwargs...) 
    apts, bpts = getdata(domain, lc; kwargs...)
    map(p -> [p; f(p...)], apts), map(p -> [p; f(p...)], bpts)
end

function getdata(domain::AbstractVector, lc::Real=0.1; writefile::Bool=false, filepath::String=randommeshpath()) 
    # Initialize gmsh 
    gmsh.initialize()

    # Add model 
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t1")
    geo = gmsh.model.geo 
    msh = gmsh.model.mesh 

    # Add points 
    for (k, pt) in enumerate(domain)
        geo.addPoint(pt..., 0, lc, k)
    end 

    # Add lines 
    idx = 1 : length(domain)
    for (k, (i, j)) in zip(idx, circshift(idx, -1)) |> enumerate
        geo.addLine(i, j, k)
    end 

    # Add surface
    geo.addCurveLoop(idx, 1)
    geo.addPlaneSurface([1], 1)

    # Mesh surface 
    geo.synchronize()
    msh.generate(2) # 2-dimensional meshing 


    # Extract bounary points 
    _, linenodetags = msh.getElementsByType(1) # Get line points 
    linenodetags = linenodetags .|> Int 
    nodes = msh.getNode.(linenodetags) .|> first 
    boundpts = combine(getindex.(nodes, 1), getindex.(nodes, 2)) |> unique

     # Extract all points 
    _, xyz = msh.getNodes()
    allx = xyz[1 : 3 : end] 
    ally = xyz[2 : 3 : end] 
    _allpts = combine(allx, ally)
    allpts = [boundpts; setdiff(_allpts, boundpts)]

    # Finalize gmsh 
    if writefile
        endswith(filepath, ".msh") ?  gmsh.write(filepath) : error("Expected `.msh` as file extension.")
    end 
    gmsh.finalize()

    # Return dataset 
    return (allpts, boundpts)
end 


# interiorpoints(f, pts::AbstractVector, lc::Real=0.01; kwargs...) = 
#     map(pnt -> [pnt; f(pnt...)], interiorpoints(pts, lc; kwargs...))

# function interiorpoints(pts::AbstractVector, lc::Real=0.01; writefile::Bool=false, filepath::String=randommeshpath()) 
#     # Initialize gmsh 
#     gmsh.initialize()

#     # Add model 
#     gmsh.option.setNumber("General.Terminal", 1)
#     gmsh.model.add("t1")
#     geo = gmsh.model.geo 
#     msh = gmsh.model.mesh 

#     # Add points 
#     for (k, pt) in enumerate(pts)
#         geo.addPoint(pt..., 0, lc, k)
#     end 

#     # Add lines 
#     idx = 1 : length(pts)
#     for (k, (i, j)) in zip(idx, circshift(idx, -1)) |> enumerate
#         geo.addLine(i, j, k)
#     end 

#     # Add surface
#     geo.addCurveLoop(idx, 1)
#     geo.addPlaneSurface([1], 1)

#     # Mesh surface 
#     geo.synchronize()
#     msh.generate(2) # 2-dimensional meshing 

#     # Collect triangulated coordinates
#     (_, gpts, _) = gmsh.model.mesh.getNodes()
#     x = gpts[1 : 3 : end]
#     y = gpts[2 : 3 : end]
#     rpts = combine(x, y)

#     # Finalize gmsh 
#     if writefile
#         endswith(filepath, ".msh") ?  gmsh.write(filepath) : error("Expected `.msh` as file extension.")
#     end 
#     gmsh.finalize()

#     # Return coordinates 
#     return rpts 
# end 


# boundarypoints(f, pts::AbstractVector, lc::Real=0.1; kwargs...) = map(p -> [p; f(p...)], boundarypoints(pts, lc; kwargs...))

# function boundarypoints(pts::AbstractVector, lc::Real=0.1; writefile::Bool=false, filepath::String=randommeshpath())
#     # Initialize gmsh 
#     gmsh.initialize()

#     # Initialize model 
#     gmsh.option.setNumber("General.Terminal", 1)
#     gmsh.model.add("t1")
#     geo = gmsh.model.geo 
#     msh = gmsh.model.mesh

#     # Add points 
#     for (k, pt) in enumerate(pts)
#         geo.addPoint(pt..., 0, lc, k)
#     end 

#     # Add lines 
#     idx = 1 : length(pts)
#     for (k, (i, j)) in zip(idx, circshift(idx, -1)) |> enumerate
#         geo.addLine(i, j, k)
#     end 

#     # Mesh surface 
#     geo.synchronize()
#     msh.generate(2)

#     # Extract boundarypoints 
#     _, _, _nodetags = msh.getElements(1) 
#     nodetags = _nodetags |> only .|> Int 
#     nodes = msh.getNode.(nodetags) .|> first 
#     rpts = combine(getindex.(nodes, 1), getindex.(nodes, 2)) |> unique

#     # Finalize gmsh 
#     if writefile
#         endswith(filepath, ".msh") ? gmsh.write(filepath) : error("Expected file extension is `.msh`")
#     end 
#     gmsh.finalize()

#     # Return points 
#     return rpts
# end 


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
