using Makie 
using GeometryBasics 
using FractalTools

# Construct 2d mesh 
msh2 = GeometryBasics.mesh(Tesselation(Rect2D(0., 0., 10., 10.), (50, 50)))

# Define functions, nodes and observables 
f1(x, y) = sin(x) * sin(y) 
f2(x, y) = sin(x) * cos(y) 
f3(x, y) = cos(x) * sin(y) 
f4(x, y) = cos(x) * cos(y) 
func = Node{Any}(f1)
msh3 = lift(func) do f
    pnts = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in msh2.position]
    fcs = faces(msh2)
    GeometryBasics.Mesh(pnts, fcs)
end 

# Consruct a figure 
fig = Figure() 

# Insert a menu 
menu = Menu(fig, options=zip(["1", "2", "3", "4"], [f1, f2, f3, f4]))
on(menu.selection) do s 
    func[] = s 
    update_cam!(ls.scene, Vec3(10.), Vec3(0.))
    center!(ls.scene)
end 
gd = fig[1, 1] = vgrid!(Label(fig, "Choose function"), menu,tellheight=false)

# Insert a mehs plot 
ls = fig[1, 2] = LScene(fig, scenekw=(camera=cam3d!, raw=false), tellheight=false, tellwidth=false)  
foo = mesh!(ls.scene, msh3, color=first.(msh3[].position))
wireframe!(ls.scene, msh3, linewidth=3)

# Display figure 
display(fig) 
