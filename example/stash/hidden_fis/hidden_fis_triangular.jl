# This file includes an example for two dimensional interpolation over triangular interpolation domain.

using FractalTools
using Triangulation
using PyPlot 
using Statistics

close("all")

datapath = joinpath(@__DIR__,  "data/gmsh_data")
filename = "t3.msh"

# Read mesh file.
mesh = read_mesh(joinpath(datapath, filename))
scale_factor = 1 / 10  # Rescale the mesh.
for node in mesh.nodes
    node.x /= scale_factor
    node.y /= scale_factor
end
meshes = triangular_partition(mesh)
regions = gmsh2matplotlibtri.(meshes)
# region = regions[1]

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

# Refine region to evaluate interpolated function
subdiv = 3
x, y = refine(gmsh2matplotlibtri(mesh), subdiv, true)

# Compute initial and real values.
func0(xi, yi) = [0., 0.]

# Perform interpolation.
α7, α8, α11, α12 = fill(0.001, 4)
interpolant, coeffs = hiddenfis(regions, z, t, α7=α7, α8=α8, α11=α11, α12=α12, func0=func0, num_iter=10, get_coeffs=true)

# Interpolation and error
real_values = func.(x, y)
real_values_z = getindex.(real_values, 1)
real_values_t = getindex.(real_values, 2)
interpolated_values = interpolant.(x, y)
interpolated_values_z = getindex.(interpolated_values, 1)
interpolated_values_t = getindex.(interpolated_values, 2)
absolute_error_z = abs.(real_values_z - interpolated_values_z)
absolute_error_t = abs.(real_values_t - interpolated_values_t)
relative_error_z = absolute_error_z ./ (abs.(real_values_z .+ eps())) * 100
relative_error_t = absolute_error_t ./ (abs.(real_values_t .+ eps())) * 100
frerr_z = mean(absolute_error_z.^2)
frerr_t = mean(absolute_error_t.^2)

# Plot the results.
xd = xcoords(mesh)
yd = ycoords(mesh)
ztd = func.(xd, yd)
zd = getindex.(ztd, 1)
td = getindex.(ztd, 2)
minx, meanx, maxx = minimum(x), mean(x), maximum(x)
miny, meany, maxy = minimum(y), mean(y), maximum(y)

# Plot the results
figure("Real and Interpolated Values z")
plot_trisurf(x, y, real_values_z)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(real_values_z), maximum(real_values_z), length=5))
tight_layout()

plot_trisurf(x, y, interpolated_values_z)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(interpolated_values_z), maximum(interpolated_values_z), length=5))
tight_layout()
display(gcf())

figure("Absolute Values z")
plot_trisurf(x, y, absolute_error_z)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(absolute_error_z), maximum(absolute_error_z), length=5))
tight_layout()
display(gcf())

figure("Relative Error z")
plot_trisurf(x, y, relative_error_z)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(relative_error_z), maximum(relative_error_z), length=5))
tight_layout()
display(PyPlot.gcf())

figure("Real and Interpolated Values t")
plot_trisurf(x, y, real_values_t)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(real_values_t), maximum(real_values_t), length=5))
tight_layout()

plot_trisurf(x, y, interpolated_values_t)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(interpolated_values_t), maximum(interpolated_values_t), length=5))
tight_layout()
display(gcf())

figure("Absolute Values t")
plot_trisurf(x, y, absolute_error_t)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(absolute_error_t), maximum(absolute_error_t), length=5))
tight_layout()
display(gcf())

figure("Relative Error t")
plot_trisurf(x, y, relative_error_t)
colors = ["black", "red", "blue", "green", "purple", "magenta"]
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(relative_error_t), maximum(relative_error_t), length=5))
tight_layout()
display(PyPlot.gcf())
