# This file includes plot recipes for 

export trisurf, topoint, tovector

# Configurations of Trisurf recipe 
@recipe(Trisurf, msh) do scene 
    Attributes(
        wireframe2 = false,
        wfcolor = :black,
        wflinewidth = 2,
        vmarkercolor = :red,
        wflinewidth3 = 3,
        vmarkercolor3 = :orange, 
        vmarkersize3 = 20,
        meshcolor3 = nothing,
        colormap = :viridis,
        visible = true
    )
end

# Main plotting function. Basically, we plot a mesh with a wireframe, with the configurations of trisurf recipe.
function plot!(plt::Trisurf) 
    msh3 = plt[1][]
    color = plt.meshcolor3[] === nothing ? last.(coordinates(msh3)) : plt.meshcolor3
    mesh!(plt, msh3, color=color, colormap=plt.colormap, visible=plt.visible)  # mesh plot 
    wireframe!(plt, msh3, linewidth=plt.wflinewidth3)   # wireframe plot 
    plt
end

# Argument conversion, 
# We convert arguments of `trisurf` function. We add a couple of methods of `trisurf` function
# Note: Always a tuple is returned from convert_arguments 

function convert_arguments(::Type{<:Trisurf}, msh2::GeometryBasics.Mesh, f)
    (GeometryBasics.Mesh([Point(pnt..., f(pnt...)) for pnt in msh2.position], faces(msh2)),)
end 

function convert_arguments(::Type{<:Trisurf}, pnts3d::AbstractVector{<:AbstractPoint{3, T}}) where {T}
    pnts2d = project(pnts3d)
    tess = spt.Delaunay(pnts2d) 
    fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    (GeometryBasics.Mesh(pnts3d, fcs),)
end 

convert_arguments(plt::Type{<:Trisurf}, pnts::AbstractVector{<:AbstractPoint{2,T}}, f::Function) where {T} = 
    convert_arguments(plt, pnts, map(pnt -> Point(f(pnt...)...), pnts))

function convert_arguments(plt::Type{<:Trisurf}, 
                           pnts2d::AbstractVector{<:AbstractPoint{2,T}}, 
                           pnts1d::AbstractVector{<:Point1}) where {T}
    convert_arguments(plt, [Point(pnt2[1], pnt2[2], pnt1[1]) for (pnt2, pnt1) in zip(pnts2d, pnts1d)])
end 

# Fallback method. If args is a vector (for example, a vector of reals), convert the it to a vector of points
convert_arguments(plt::Type{<:Trisurf}, args::AbstractVector...) = convert_arguments(plt, map(arg -> topoint.(arg), args)...)

topoint(pnt::AbstractPoint) = pnt
topoint(vect::AbstractVector{<:Real}) = Point(vect...)
topoint(vect::Real) = Point(vect)

tovector(pnt::AbstractPoint) = [pnt...]
tovector(vect::AbstractVector{<:Real}) = vect

