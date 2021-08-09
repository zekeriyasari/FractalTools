# This file includes methods for data generation. 

export Dataset, Tessellation, tessellate, getdata, getpoint, project, uniformdomain

# --------------------------------------------- Tessellation --------------------------------- # 

struct Tessellation{T} 
    tess::T 
end 

function Tessellation(points, domain=[]) 
    if isempty(domain) 
        tess = spt.Delaunay(points)
    elseif length(domain) == 2 
        tess = LineString([Point(point) for point in points])
    else 
        triin = Triangulate.TriangulateIO()
        triin.pointlist = hcat(points...)
        list = 1 : length(domain)
        triin.segmentmarkerlist = collect(list) 
        triin.segmentlist = hcat(collect.(zip(val, circshift(val, -1)))...)
        triout, vorout = Triangulate.triangulate("pQ", triin)
        tess = tri.Triangulation(triout.pointlist[1, :], triout.pointlist[2, :], triout.trianglelist' .- 1)
    end
    Tessellation(tess)
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
        vcat(boundpoints, interpoints) |> unique
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

tessellate(dataset::Dataset) = Tessellation(dataset.domain, dataset.points)

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

tomesh(dataset::Dataset) = tomesh(dataset.point, dataset.domain)
function tomesh(points::AbstractVector, domain=[])
    n = length(points[1]) 
    if n > 2 
        points = project(points, n - 2)
    end 
    if isempty(domain)
        tess = Tessellation(domain, points)
        faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.triangles .+ 1)]
    else 
        tess = spt.Delaunay(points) 
        faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    end 
    GeometryBasics.Mesh(points, faces)
end


_project(point, n=1) = point[1 : end - n]

project(points::AbstractVector, n=1) = _project.(points, n)
project(dataset::Dataset, n=1) = Dataset(project(dataset.points, n), dataset.domain)
