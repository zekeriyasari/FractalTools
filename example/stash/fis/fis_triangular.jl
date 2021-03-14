# This file includes an example for two dimensional interpolation over triangular interpolation domain.

# TODO: Relative error must decrease with decreasing vertical scaling factor and increasing interpolation point. 
# TODO: Plot error levels with colors in 2D surfaces.
# TODO: Increase number of points in the mesh.
# TODO: Relative error calculation for zero. Look at MSE with scaling factor.
# TODO: Construct the interpolation domain while constructing the interpolant.

using FractalTools
using Triangulation
using PyPlot 
using Statistics

close("all")

datapath = joinpath(@__DIR__,  "../../data/gmsh_data")
# datapath = joinpath(pathof(FractalTools), "data/gmsh_data")
filename = "t3.msh"

# Read mesh file.
mesh = read_mesh(joinpath(datapath, filename))
scale_factor = 1   # Rescale the mesh.
for node in mesh.nodes
    node.x /= scale_factor
    node.y /= scale_factor
end
meshes = triangular_partition(mesh)
regions = gmsh2matplotlibtri.(meshes)

# Compute z-values on mesh points.
func(x, y) = x.^2 .+ y.^2 .+ 1.
# func(x, y) = -20 * exp.(-0.2 * sqrt.(0.5 * (x.^2 .+ y.^2))) - exp.(0.5 * (cos.(2 * pi * x) + cos.(2 * pi * y))) .+ MathConstants.e .+ 20
# func(x, y) = -(y .+ 47) .* sin.(sqrt.(abs.(x ./ 2 .+ (y .+ 47)))) .- x .* sin.(sqrt.(abs.(x .- (y .+ 47))))
# func(x, y) = 100 * sqrt.(abs.(y - 0.01 * x.^2)) .+ 0.01 * abs.(x .+ 10)
number_of_regions = length(regions)
z = Vector{Vector}(undef, number_of_regions)
for k = 1 : number_of_regions
    z[k] = func(regions[k].x, regions[k].y)
end

# Refine region to evaluate interpolated function
subdiv = 3
x, y = refine(gmsh2matplotlibtri(mesh), subdiv, true)

# Compute initial and real values.
func0(xi, yi) = 0.

# Perform interpolation.
# α = 0.025
# α = 0.151
α = 0.001
interpolant, coeffs = fis(regions, z, α=α, func0=func0, num_iter=10, get_coeffs=true)

# Interpolation and error
real_values = func.(x, y)
interpolated_values = interpolant.(x, y)
absolute_error = abs.(interpolated_values - real_values)
relative_error = absolute_error ./ (abs.(real_values)) * 100
frerr = mean(absolute_error.^2)

# Filtered relative error 
threshold = 0.05
mask = x.^2 .+ y.^2 .< threshold
idx = findall(mask)
relative_error_fltr = copy(relative_error)
relative_error_fltr[idx] .= 0.# TODO: NAN value

# # Plot the results.
xd = xcoords(mesh)
yd = ycoords(mesh)
zd = func.(xd, yd)
minx, meanx, maxx = minimum(x), mean(x), maximum(x)
miny, meany, maxy = minimum(y), mean(y), maximum(y)

# Plot the results
figure("Real Values")
plot_trisurf(x, y, real_values)

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
zticks(range(minimum(real_values), maximum(real_values), length=5))
tight_layout()
display(PyPlot.gcf())

figure("Interpolated Values")
plot_trisurf(x, y, interpolated_values)
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
zticks(range(minimum(interpolated_values), maximum(interpolated_values), length=5))
tight_layout()
display(PyPlot.gcf())

figure("Absolute Values")
plot_trisurf(x, y, absolute_error)
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
zticks(range(minimum(absolute_error), maximum(absolute_error), length=5))
tight_layout()
display(PyPlot.gcf())

figure("Relative Error")
plot_trisurf(x, y, relative_error)
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(relative_error), maximum(relative_error), length=5))
tight_layout()
display(PyPlot.gcf())

figure("Relative Error Filtered")
plot_trisurf(x, y, relative_error_fltr)
for (color, region) in zip(colors, regions)
    triplot(region, color=color)
    plot(region.x, region.y, ".", color="orange")
end
xlabel("x", labelpad=10, fontsize=12)
ylabel("y", labelpad=10, fontsize=12)
zlabel("z", labelpad=10, fontsize=12)
xticks(range(minx, maxx, length=5))
yticks(range(miny, maxy, length=5))
# zticks(range(minimum(relative_error_fltr), maximum(relative_error_fltr), length=5))
tight_layout()
display(PyPlot.gcf())

# Plot interpolation domain.
fig, ax = subplots(1)
for region in regions
    ax.triplot(region)
    ax.plot(region.x, region.y, ".", color="orange")
end
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("y", fontsize=12)
tight_layout()
display(PyPlot.gcf())
