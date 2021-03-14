using Makie 
using GeometryBasics 
using PyCall 
using FractalTools
using StaticArrays 
using Random 

# Construct interpolation data 
hexagon = GeometryBasics.Ngon(@SVector [Point(BigFloat(cos(θ)), BigFloat(sin(θ))) for θ in 0 : π / 3 : 2π - π / 3])
f(x, y) = x^2 + y^2
pnts3d = getdata(f, hexagon, 100)

# Construct interpolant 
itp = interpolate(pnts3d)

function erf(xd, yd)
    fval = f(xd, yd) 
    ival = itp(xd, yd) 
    abs(fval - ival) / abs(fval) * 100
end

# Calculate error values 
pnts2d = project(pnts3d)
ivals = map(pnt -> Point(pnt[1], pnt[2], itp(pnt...)), pnts2d)
fvals = map(pnt -> Point(pnt[1], pnt[2], f(pnt...)), pnts2d)
evals = map(pnt -> Point(pnt[1], pnt[2], erf(pnt...)), pnts2d)

# Construct GUI
fig         = Figure() 
lb          = Label(fig, "Show Interpolant") 
tog         = Toggle(fig, active=false) 
ax1         = LScene(fig, scenekw=(camema=cam3d!, raw=false), tellwidth=false)
fig[1, 1]   = vgrid!(ax1, hgrid!(lb, tog, tellwidth=false), tellwidth=false)
ax2         = fig[1, 2] = LScene(fig, scenekw=(camema=cam3d!, raw=false), tellwidth=false)

        trisurf!(ax1, fvals, meshcolor3=last.(fvals), colormap=:viridis)
plt =   trisurf!(ax1, ivals, meshcolor3=last.(ivals), colormap=:heat)
        wireframe!(ax1, project(triangulate(pnts3d)[2]))  
        trisurf!(ax2, evals, meshcolor3=last.(evals), colormap=:viridis)

on(tog.active) do val 
    plt.visible[] = val
end 

tog.active[] = false

display(fig)

# # Record video 
# fps = 24
# record(fig.scene, "errorplot.mp4"; framerate = fps) do io
#     for i = 1 : 10 * fps
#         sleep(1/fps)
#         recordframe!(io)
#     end
# end
