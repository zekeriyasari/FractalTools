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
region = regions[3]


# Compute z-values on mesh points.
# Compute z-values on mesh points.
func(x, y) = [
    x^2 + y^2,
    x^2 - y^2
]
number_of_regions = length(regions)
z = Vector{Vector}(undef, number_of_regions)
t = Vector{Vector}(undef, number_of_regions)
for k = 1 : number_of_regions
    z[k] = getindex.(func.(regions[k].x, regions[k].y), 1)
    t[k] = getindex.(func.(regions[k].x, regions[k].y), 2)
end


# Calculate integral
I = map(i ->  integrate(regions[i], z[i], t[i], α7=0.001, α8=0.001, α11=0.001, α12=0.001), 1 : length(z))
@show I
#------------------------ Plots --------------------#

figure("interpolation")
plot_trisurf(region.x ,region.y, z)
triplot(region)
gcf()


