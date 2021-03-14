using FractalTools
using PyPlot
using Triangulation

close("all")

datapath = joinpath(@__DIR__,  "../../data/gmsh_data")
# datapath = joinpath(pathof(FractalTools), "data/gmsh_data")
filename = "t3.msh"

# Read mesh file.
mesh = read_mesh(joinpath(datapath, filename))
meshes = triangular_partition(mesh)
regions = gmsh2matplotlibtri.(meshes) 
# region = regions[3]


# Compute z-values on mesh points.
f(x, y) = x.^2 .+ y.^2
z = map(region ->  f(region.x,region.y), regions)



# Calculate integral
I = map(i ->  integrate(regions[i], z[i], α=0.001), 1 : length(z))
@show I
#------------------------ Plots --------------------#

figure("interpolation")
plot_trisurf(region.x ,region.y, z)
triplot(region)
gcf()


