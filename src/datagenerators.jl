# This file includes methods for data generation. 

export Dataset, tessellate, mesh, uniformdomain, project 

# -------------------------------------- Dataset -------------------------------------- # 

struct Dataset{T1<:AbstractVector, T2<:AbstractVector}
    domain::T1 
    points::T2
end 

function Dataset(domain::AbstractVector, numpoints::Int) 
    if length(domain) == 2 
        points = range(domain[1], domain[end], length=numpoints) |> collect
    else
        sampledpoints = Meshes.sample(
            Meshes.Ngon(Meshes.Point2.(domain)), 
            Meshes.HomogeneousSampling(numpoints)
        ) |> collect .|> Meshes.coordinates .|> collect
        points = [domain; sampledpoints]
    end
    Dataset(domain, points)
end 

function Dataset(f, domain::AbstractVector, numpoints::Int) 
    ds = Dataset(domain, numpoints)
    points = map(pnt -> [pnt; f(pnt...)], ds.points) 
    Dataset(domain,  points)
end 

uniformdomain(n, T=Float64, p0=zeros(T, 2), r=1.) = [p0 + T[r * cos(θ), r*sin(θ)] for θ in (0 : n - 1) / n * 2π]

# -------------------------------------- Tessellation -------------------------------------- # 

struct Tessellation{T}
    tess::T
end 

tessellate(dataset::Dataset) = tessellate(dataset.points, dataset.domain)
function tessellate(pts::AbstractVector, domain::AbstractVector)
    triin = Triangulate.TriangulateIO()
    triin.pointlist = hcat(pts...)
    N = length(domain)
    triin.segmentlist = hcat([[i, i + 1] for i in 1 : N - 1]..., [N, 1])
    triin.segmentmarkerlist = collect(1 : N)
    triout, vorout = triangulate("pQ", triin)
    tess = tri.Triangulation(triout.pointlist[1, :], triout.pointlist[2, :], triout.trianglelist' .- 1)
    Tessellation(tess)
end

# -------------------------------------- Meshing -------------------------------------- # 

"""
    $SIGNATURES

Returns a mesh whose points are `pnts`
"""
function mesh(dataset::Dataset)
    points = dataset.points
    domain = dataset.domain
    GeometryBasics.Mesh(topoint.(points), tofaces(topoint.(points), domain))
end 

tofaces(pnts::AbstractVector{<:AbstractPoint{N,T}}, domain) where {N,T} = tofaces(project(pnts, N - 2), domain)
function tofaces(pnts::AbstractVector{<:AbstractPoint{2,T}}, domain) where {T}
    tessellation = tessellate(pnts, domain)
    [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tessellation.tess.triangles .+ 1)]
end 

topoint(pnt::AbstractPoint) = pnt
topoint(vect::AbstractVector{<:Real}) = Point(vect...)
topoint(vect::Real) = Point(vect)

tovector(pnt::AbstractPoint) = [pnt...]
tovector(vect::AbstractVector{<:Real}) = vect

# -------------------------------------- Projection -------------------------------------- # 

_project(pnt, n) = pnt[1 : end - n]
project(dataset::Dataset, n::Int) = Dataset(dataset.domain, _project.(dataset.points, n))
