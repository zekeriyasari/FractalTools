# Import Meshes.atol to package 
# Meshes.atol(::Type{BigFloat}) = 100 * eps()

using FractalTools 
using GLMakie 

# Dataset 
f(x, y) = x^2 + y^2 + 1
Ω = [[0., 0.], [1., 1.], [2., 0.], [3., 1.], 
    [3., 0.], [2., -1.], [1., 0.], [0., -1.]]
bpts = generate(f, Ω, 0.1, interior=false)
bpts = bpts[length(Ω) + 1 : end]
pts = generate(f, Ω, 0.1) 

fig, ax, plt = scatter(getindex.(bpts, 1), getindex.(bpts, 2))

fig = Figure() 
ax = Axis(fig[1, 1]) 
for bpt in bpts 
    scatter!(ax, [bpt[1]], [bpt[2]])
    display(fig) 
    sleep(0.25) 
end 
