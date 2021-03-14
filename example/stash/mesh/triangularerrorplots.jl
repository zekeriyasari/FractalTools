# This script includes a demo for the errot plots. 

using Makie 
using FractalTools 
using GeometryBasics
using Random 

Random.seed!(0)

# Generate data 
trig = Triangle(
    Point(BigFloat(-1.), BigFloat(-1.)), 
    Point(BigFloat(1.), BigFloat(-1)),
    Point(BigFloat(0.), BigFloat(2.)) 
)
f(x, y) = x^2 + y^2
pnts3d = getdata(f, trig, 100)
 
# Construct interpolant 
itp = interpolate(pnts3d)
function erf(x, y)
    ival = itp(x, y) 
    fval = f(x, y) 
    abs(fval - ival) / abs(fval) * 100
end

# Construct gui
_, msh3 = triangulate(pnts3d)
msh2 = project(msh3)
pnts2d = project(pnts3d)   

fig         = Figure() 
lb          = Label(fig, "Show Interpolant") 
tog         = Toggle(fig, active=false) 
ax1         = LScene(fig, scenekw=(camema=cam3d!, raw=false), tellwidth=false)
fig[1, 1]   = vgrid!(ax1, hgrid!(lb, tog, tellwidth=false), tellwidth=false)
ax2         = fig[1, 2] = LScene(fig, scenekw=(camema=cam3d!, raw=false), tellwidth=false)

        trisurf!(ax1, msh2, f, meshcolor3=first.(msh3.position), colormap=:viridis)
plt =   trisurf!(ax1, msh2, itp, meshcolor3=first.(msh3.position), colormap=:heat)
        wireframe!(ax1, msh2)  
        trisurf!(ax2, msh2, erf, meshcolor3=first.(msh2.position))

on(tog.active) do val 
    plt.visible[] = val
end 

tog.active[] = false

display(fig)

# fps = 60
# record(fig.scene, "errorplot.mp4"; framerate = fps) do io
#     for i = 1 : 10 * fps
#         sleep(1/fps)
#         recordframe!(io)
#     end
# end