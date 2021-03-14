using GeometryBasics 
using Makie 
using FractalTools 
using Random 
Random.seed!(0)

# Define the triangle 
trig = Triangle(
    Point(BigFloat(-1.), BigFloat(-1.)), 
    Point(BigFloat(1.), BigFloat(-1.)), 
    Point(BigFloat(0.), BigFloat(2))
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

# Add observables for interactivity
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

# Plot figures 
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

# Button intrections
on(button1.clicks) do n 
    i[] += 1
    label.text[] = "Tri idx: " * string(i[])
end 

on(button2.clicks) do n 
    i[] -= 1
    label.text[] = "Tri idx: " * string(i[])
end 

display(fig) 
