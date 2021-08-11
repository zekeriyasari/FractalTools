# This file includes plot recipes for 

export trisurf

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
function Makie.plot!(plt::Trisurf) 
    msh3 = plt[1][]
    color = plt.meshcolor3[] === nothing ? last.(coordinates(msh3)) : plt.meshcolor3
    mesh!(plt, msh3, color=color, colormap=plt.colormap, visible=plt.visible)  # mesh plot 
    wireframe!(plt, msh3, linewidth=plt.wflinewidth3)   # wireframe plot 
    plt
end

# Argument conversion, 
# We convert arguments of `trisurf` function. We add a couple of methods of `trisurf` function
# Note: Always a tuple is returned from convert_arguments 

# Makie.convert_arguments(::Type{<:Trisurf}, pnts::AbstractVector...)  = (tomesh(combine(pnts...)),)
Makie.convert_arguments(::Type{<:Trisurf}, dataset::Dataset) = (tomesh(dataset),)
Makie.convert_arguments(::Type{<:Trisurf}, pnts1::AbstractVector, pnts2::AbstractVector, domain::AbstractVector) = 
    (tomesh(combine(pnts1, pnts2), domain),)

Makie.convert_arguments(plt::Type{<:Trisurf}, pnts::AbstractVector, f::Function) =  
    convert_arguments(plt, map(pnt -> Point(pnt..., f(pnt...)...), pnts))

Makie.convert_arguments(::Type{<:Trisurf}, msh2::GeometryBasics.Mesh, f::Function) = 
    (GeometryBasics.Mesh([Point(pnt..., f(pnt...)) for pnt in msh2.position], faces(msh2)),)
