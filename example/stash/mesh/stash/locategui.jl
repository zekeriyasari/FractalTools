using FractalTools 
using Makie 

# Contruct a tesslation 
tridln = Node(TriDelaunay([0., 0.], [1., 0], [0.5, 1.], addboundarypoints=true))

# Definea an observables
msh = lift(tridln) do val 
    tomesh(val.delaunay)
end 
vtx = lift(msh) do val 
    val.position 
end 

# Test point 
tpnt = Node(getpoint(tridln[]))
tripts = lift(tpnt) do val 
    tidx = simplices(tridln[])[locate(tridln[], val), :]
    map(val -> Point(val...), eachrow(points(tridln[])[tidx, :]))
end 

# Plot tesselation 
fig = Figure(resolution=(500, 500), scale_plot = false, camera = AbstractPlotting.campixel!)
ax = fig[1, 1] =  Axis(fig) 
mesh!(ax, msh, color=:lightblue)
wireframe!(ax, msh, color=:black)
scatter!(ax, vtx, color=:black)
scatter!(ax, tripts, color=:red, markersize=10)
scatter!(ax, tpnt, color=:orange, markersize=10)

# Add mouse button interaction 
on(ax.scene.events.mousebuttons) do buttons
   if ispressed(ax.scene, Mouse.left)
        area = pixelarea(ax.scene)[]
        mp = ax.scene.events.mouseposition[]
        mp = Point2f0(mp) .- minimum(area)
        pos = to_world(ax.scene, mp)
        @info "$pos added"
        tpnt[] = pos
   end
end

# Add button interaction.
on(ax.scene.events.keyboardbuttons) do button 
    if ispressed(button, Keyboard.a)
        pnt = addpoint!(tridln[])
        @info "$pnt is added"
    end 
    if ispressed(button, Keyboard.f)
        finegrain!(tridln[], 10)
        @info "delaunay is fine grained"
    end 
    tridln[] = tridln[]   # Trigger change
end 

# Display scene 
display(fig)
