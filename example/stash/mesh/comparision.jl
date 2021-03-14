# This file is used for comparison of fractal interpolation and spline intepolation 

using Makie 
using GeometryBasics 
using FractalTools
import Interpolations
import Interpolations: Gridded, Linear

# Define test functions 
f(x, y) = x^2 + y^2 + 1

# Construct gridded data to construct and test interpolant 
ndata = 10 
ntest = 20
data_mesh_2d = GeometryBasics.mesh(Tesselation(Rect2D(BigFloat.([-1., -1, 2., 2.])...), (ndata, ndata)))
data_points_3d = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in data_mesh_2d.position]     # points to construct interpolants 
test_mesh_2d = GeometryBasics.mesh(Tesselation(Rect2D(BigFloat.([-1., -1, 2., 2.])...), (ntest, ntest)))
test_points_3d = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in test_mesh_2d.position]    # points to test interpolants

# Construct interpolants 
fractal_interpolant = FractalTools.interpolate(data_points_3d)
spline_interpolant = Interpolations.interpolate(
    (unique(getindex.(data_points_3d, 1)), unique(getindex.(data_points_3d, 2))),       # x-y gridded data 
    collect(reshape(getindex.(data_points_3d, 3), ndata, ndata)),                       # z values 
    Gridded(Linear())                                                                   # Interpolation type 
)

# Evaluations 
true_data_points                        = map(pnt -> Point(f(pnt[1], pnt[2])), data_points_3d)
true_test_points                        = map(pnt -> Point(f(pnt[1], pnt[2])), test_points_3d)

fractal_interpolant_data_points         = map(pnt -> Point(fractal_interpolant(pnt[1], pnt[2])), data_points_3d)
fractal_interpolant_test_points         = map(pnt -> Point(fractal_interpolant(pnt[1], pnt[2])), test_points_3d)
fractal_interpolant_data_points_error   = Point.( 
    abs.(last.(true_data_points - fractal_interpolant_data_points)) ./ abs.(last.(true_data_points)) * 100
    )
fractal_interpolant_test_points_error   = Point.( 
    abs.(last.(true_test_points - fractal_interpolant_test_points)) ./ abs.(last.(true_test_points)) * 100
    )

spline_interpolant_data_points          = map(pnt -> Point(spline_interpolant(pnt[1], pnt[2])), data_points_3d)
spline_interpolant_test_points          = map(pnt -> Point(fractal_interpolant(pnt[1], pnt[2])), test_points_3d)
spline_interpolant_data_points_error   = Point.( 
    abs.(last.(true_data_points - spline_interpolant_data_points)) ./ abs.(last.(true_data_points)) * 100
    )
spline_interpolant_test_points_error    = Point.( 
    abs.(last.(true_test_points - spline_interpolant_test_points)) ./ abs.(last.(true_test_points)) * 100
    )

# Plots 
fig = Figure() 
ls11 = LScene(fig[1, 1])
ls12 = LScene(fig[1, 2])
ls13 = LScene(fig[1, 3])
ls21 = LScene(fig[2, 1])
ls22 = LScene(fig[2, 2])
ls23 = LScene(fig[2, 3])

trisurf!(ls11, data_mesh_2d.position, true_data_points, meshcolor3=last.(true_data_points))
trisurf!(ls11, data_mesh_2d.position, fractal_interpolant_data_points, meshcolor3=last.(fractal_interpolant_data_points), colormap=:heat)
wireframe!(ls11, data_mesh_2d)

trisurf!(ls12, data_mesh_2d.position, fractal_interpolant_data_points_error, meshcolor3=last.(fractal_interpolant_data_points_error))
wireframe!(ls12, data_mesh_2d)

trisurf!(ls13, test_mesh_2d.position, fractal_interpolant_test_points_error, meshcolor3=last.(fractal_interpolant_test_points_error))
wireframe!(ls13, test_mesh_2d)

trisurf!(ls21, data_mesh_2d.position, true_data_points, meshcolor3=last.(true_data_points))
trisurf!(ls21, data_mesh_2d.position, spline_interpolant_data_points, meshcolor3=last.(spline_interpolant_data_points), colormap=:heat)
wireframe!(ls21, data_mesh_2d)
trisurf!(ls22, data_mesh_2d.position, spline_interpolant_data_points_error, meshcolor3=last.(spline_interpolant_data_points_error))
wireframe!(ls22, data_mesh_2d)
trisurf!(ls23, test_mesh_2d.position, spline_interpolant_test_points_error, meshcolor3=last.(spline_interpolant_test_points_error))
wireframe!(ls23, test_mesh_2d)

display(fig)
