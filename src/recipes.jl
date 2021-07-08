# This file includes plot recipes for 

export trisurf, topoint, tovector, tomesh 

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

convert_arguments(::Type{<:Trisurf}, pnts::AbstractVector...) where {T} = (tomesh(combine.(pnts...)),)

convert_arguments(plt::Type{<:Trisurf}, pnts::AbstractVector, f::Function) =  
    convert_arguments(plt, map(pnt -> Point(pnt..., f(pnt...)...), pnts))

convert_arguments(::Type{<:Trisurf}, msh2::GeometryBasics.Mesh, f::Function) = 
    (GeometryBasics.Mesh([Point(pnt..., f(pnt...)) for pnt in msh2.position], faces(msh2)),)

combine(ps...) = vcat([[p...] for p in ps]...)
