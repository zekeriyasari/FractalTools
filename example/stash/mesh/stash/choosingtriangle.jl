using FractalTools 
using Makie 
using GeometryBasics

# Construct a mesh 
pts = [
    [0., 0.], 
    [1., 0.], 
    [1., 1.] 
]
msh = trimesh(pts...)
trange = 1 : length(faces(msh))

# Plot mesh 
fig = Figure() 
ax, plt = mesh(fig[1, 1], msh, color=last.(msh.position)) 
wireframe!(ax, msh)
cb = Colorbar(fig[1, 2], plt, width=20)
button1 = Button(fig[2, 1], tellwidth=false, label="Increase")
button2 = Button(fig[2, 2], tellwidth=false, label="Decrease")
label = Label(fig[3, :], tellwidth=false)

# Plot a triangle 
coords = coordinates(msh)
trigs = faces(msh) 
idx = Node(1)
x = @lift(getindex.(coords[trigs[$idx]], 1) |> collect)
y = @lift(getindex.(coords[trigs[$idx]], 2) |> collect)
scatter!(ax, x, y)

on(button1.clicks) do n 
    idx[] += 1 
    label.text = "idx = $(idx[]), x = $(x[]), y = $(y[])"
end 

on(button2.clicks) do n 
    idx[] -= 1 
    label.text = "idx = $(idx[]), x = $(x[]), y = $(y[])"
end 

# Bind keybord interacton 
scene = fig.scene 
on(scene.events.keyboardbuttons) do button
    if ispressed(button, Keyboard.up)
        idx[] += 1
        label.text = "idx = $(idx[]), x = $(x[]), y = $(y[])"
    end 
    if ispressed(button, Keyboard.down)
        idx[] -= 1
        label.text = "idx = $(idx[]), x = $(x[]), y = $(y[])"
    end 
end

display(fig) 
