# This file includes plot recipes for Makie

export trisurf

@recipe(Trisurf, msh, f) do scene 
    Makie.Attributes(
        wireframe2 = false,
        wfcolor = :black,
        wflinewidth = 2,
        vmarkercolor = :red,
        wflinewidth3 = 3,
        vmarkercolor3 = :orange, 
        vmarkersize3 = 20,
        meshcolor3 = :black,
        colormap = :viridis,
        visible = true
    )
end

function Makie.plot!(plt::Trisurf) 
    msh3 = plt[1][]
    Makie.mesh!(plt, msh3, color=plt.meshcolor3, colormap=plt.colormap, visible=plt.visible) 
    Makie.wireframe!(plt, msh3, linewidth=plt.wflinewidth3)
    plt
end

function Makie.convert_arguments(::Type{<:Trisurf}, msh2::GeometryBasics.Mesh, f)
    msh3 = GeometryBasics.Mesh([Point(pnt..., f(pnt...)) for pnt in msh2.position], faces(msh2))
    return (msh3, msh2, f)
end 

function Makie.convert_arguments(plt::Type{<:Trisurf}, 
                                            pnts::AbstractVector{<:AbstractPoint{2,T}}, 
                                            f::Function) where T
    Makie.convert_arguments(plt, pnts, map(pnt -> Point(f(pnt...)...), pnts))
end 

function Makie.convert_arguments(plt::Type{<:Trisurf}, 
                                            pnts2d::AbstractVector{<:AbstractPoint{2,T}}, 
                                            pnts1d::AbstractVector{<:Point1}) where {T} 
    Makie.convert_arguments(plt, [Point(pnt2[1], pnt2[2], pnt1[1]) for (pnt2, pnt1) in zip(pnts2d, pnts1d)])
end 

function Makie.convert_arguments(::Type{<:Trisurf}, 
                                            pnts3d::AbstractVector{<:AbstractPoint{3, T}}) where {T}
    pnts2d = project(pnts3d)
    tess = spt.Delaunay(pnts2d) 
    fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    msh3 = GeometryBasics.Mesh(pnts3d, fcs)
    return (msh3, pnts3d)
end 

