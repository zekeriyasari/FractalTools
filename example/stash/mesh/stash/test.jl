using Makie 
using GeometryBasics

outtrig = Triangle(
    Point(0., 0.), 
    Point(1., 0.),
    Point(0.5, 1.)
)

intrig = Triangle(
    Point(rand(2)...), 
    Point(rand(2)...), 
    Point(rand(2)...)
)

fig = Figure() 
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
mesh!(ax1, outtrig.points)
mesh!(ax2, intrig.points)
using GeometryBasics 
using Makie 
using FractalTools 

using Random 
Random.seed!(0)


# Define the triangle 
trig = Triangle(
    Point(-1., -1.), 
    Point(1., -1.), 
    Point(0., 2)
)

# Get valid points in the domain 
pnts2d = disperse(trig, 10)

# Define function 
f(x, y) = x^2 + y^2 

# Evaluate function 
pnts3d = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in pnts2d]  # Interpolation data 

# Check transformations 
_, msh3 = triangulate(pnts3d)
outtrig = findouttriangle(pnts3d)
outtrigpnts2d = project(outtrig.points)

i = Node(39) 
ttrig = @lift(msh3[$i])

ttrigpnts2d = lift(ttrig) do ttrig 
    project(ttrig.points)
end 

mtrig = lift(ttrig) do ttrig 
    tfm = gettransform(outtrig, ttrig, 0.01)
    L(x) = tfm.A[1:2, 1:2] * x + tfm.b[1:2]
    Triangle(map(pnt -> Point(L(pnt)...), outtrigpnts2d)...)
end 

fig = Figure() 
ls = LScene(fig[1, 1])
buttongrid = fig[1, 2] =  GridLayout(tellheight=false)
button1 = buttongrid[1, 1] = Button(fig, label = "Inc")
button2 = buttongrid[2, 1] = Button(fig, label = "Dec")
label = buttongrid[3, 1] = Label(fig, "Tri idx:$(string(i[]))")
mesh!(ls, outtrig)
wireframe!(outtrigpnts2d, color=:black, linewidth=3)
mesh!(ttrig, color=:blue)
wireframe!(ttrigpnts2d, color=:blue, linewidth=3)
mesh!(mtrig, color=:green)

on(button1.clicks) do n 
    i[] += 1
    label.text[] = "Tri idx: " * string(i[])
end 

on(button2.clicks) do n 
    i[] -= 1
    label.text[] = "Tri idx: " * string(i[])
end 

display(fig) 

# # Interpolate the data
# itp = interpolate(pnts3d, maxiters=4)

# # Evaluat the itp 
# tpnt = pnts2d
# # tpnt = [getpoint(trig) for i in 1 : 500]
# val = map(pnt -> itp(pnt...), tpnt) 
# rval =  map(pnt -> f(pnt...), tpnt) 
# relerr = abs.(val - rval) ./ abs.(rval) * 100

# @show findall(isnan, val) 
# @show sum(relerr .≥ 1)


# # Plots 
# tess, mesh3 = triangulate(pnts3d)
# mesh2 = GeometryBasics.Mesh(project(coordinates(mesh3)), faces(mesh3))
# fig, ax, plt = scatter(pnts3d, color=:red, markersize=20)
# wireframe!(trig.points, linewidth=3)
# wireframe!(mesh2) 
# scatter!(pnts2d, markersize=8, color=:black)
# scatter!(tpnt, markersize=8, color=:green)
# # mesh!(mesh3, color=first.(mesh3.position))
# wireframe!(mesh3, color=:black)

# # Test point 
# tpnt = getpoint(trig)      
# idx = tess.find_simplex(tpnt)[1] + 1  
# scatter!(tpnt, markersize=30, color=:orange)
# mesh!(mesh2[idx])

# tpnt2 = (trig.points[1] + trig.points[2]) / 2 
# idx = tess.find_simplex(tpnt2)[1] + 1  
# scatter!(tpnt2, markersize=30, color=:blue)
# mesh!(mesh2[idx], color=:green)

# tpnt3 = pnts2d[40]
# idx = tess.find_simplex(tpnt3)[1] + 1  
# scatter!(tpnt3, markersize=30, color=:green)
# mesh!(mesh2[idx], color=:black)

# tpnt4 = pnts2d[39]
# scatter!(tpnt4, markersize=30, color=:blue)

# tpnt5 = (pnts2d[40] + pnts2d[39]) / 2 
# scatter!(tpnt5, markersize=30, color=:blue)
# idx = tess.find_simplex(tpnt5)[1] + 1  
# mesh!(mesh2[idx])

# display(fig) 