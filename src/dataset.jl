# This file includes methods for data generation. 

export getdata, tomesh, uniformdomain

getdata(f, domain::AbstractVector, lc::Real=0.1; kwargs...) = map(pnt -> [pnt; f(pnt...)], getdata(domain, lc; kwargs...))

function getdata(domain::AbstractVector, lc::Real=0.1; kwargs...)
    if length(first(domain)) == 1 
        getdata1d(domain, lc; kwargs...) 
    else  
        getdata2d(domain, lc; kwargs...)
    end 
end 

function getdata1d(domain::AbstractVector, lc::Real=0.1; writefile::Bool=false, filepath::String=randommeshpath())
    # Initialize gmsh 
    gmsh.initialize()

    # Add model 
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t1")
    geo = gmsh.model.geo 
    msh = gmsh.model.mesh 

    # Add points 
    geo.addPoint(only(domain[1]), 0, 0, lc, 1) 
    geo.addPoint(only(domain[2]), 0, 0, lc, 2) 

    # Add lines 
    geo.addLine(1, 2, 1) 

    # Mesh surface 
    geo.synchronize()
    msh.generate(2) # 2-dimensional meshing 

     # Extract all points 
    _, xyz = msh.getNodes()
    allpts = xyz[1 : 3 : end] 

    # Finalize gmsh 
    if writefile
        endswith(filepath, ".msh") ?  gmsh.write(filepath) : error("Expected `.msh` as file extension.")
    end 
    gmsh.finalize()

    # Return dataset 
    return [[item] for item in sort(allpts)]
end 

function getdata2d(domain::AbstractVector, lc::Real=0.1; writefile::Bool=false, filepath::String=randommeshpath()) 
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

    # # Extract bounary points 
    # _, linenodetags = msh.getElementsByType(1) # Get line points 
    # linenodetags = linenodetags .|> Int 
    # nodes = msh.getNode.(linenodetags) .|> first 
    # boundpts = combine(getindex.(nodes, 1), getindex.(nodes, 2)) |> unique

    #  # Extract all points 
    # _, xyz = msh.getNodes()
    # allx = xyz[1 : 3 : end] 
    # ally = xyz[2 : 3 : end] 
    # _allpts = combine(allx, ally)
    # allpts = [boundpts; setdiff(_allpts, boundpts)]

     # Extract all points 
    _, xyz = msh.getNodes()
    allpts = combine(xyz[1 : 3 : end] ,  xyz[2 : 3 : end])

    # Finalize gmsh 
    if writefile
        endswith(filepath, ".msh") ?  gmsh.write(filepath) : error("Expected `.msh` as file extension.")
    end 
    gmsh.finalize()

    # Return dataset 
    return allpts 
end 

"""
    $SIGNATURES

Returns a uniform ngon whose vertex points are centered at `p0` with a radius `r`.
"""
uniformdomain(n::Int, T::Type{<:Real}=Float64, p0::AbstractVector{<:Real}=zeros(T, 2), r::Real=1.) =
    [p0 + T[r * cos(θ), r*sin(θ)] for θ in (0 : n - 1) / n * 2π] 

"""
    $SIGNATURES

Returns a mesh whose points are `pnts`
"""
tomesh(pnts::AbstractVector) = GeometryBasics.Mesh(topoint.(pnts), tofaces(topoint.(pnts)))

tofaces(pnts::AbstractVector{<:AbstractPoint{N,T}}) where {N,T} = tofaces(project(pnts, N - 2))
function tofaces(pnts::AbstractVector{<:AbstractPoint{2,T}}) where {T}
    tess = spt.Delaunay(pnts) 
    [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
end 

topoint(pnt::AbstractPoint) = pnt
topoint(vect::AbstractVector{<:Real}) = Point(vect...)
topoint(vect::Real) = Point(vect)

# tovector(pnt::AbstractPoint) = [pnt...]
tovector(vect::AbstractVector{<:Real}) = vect


