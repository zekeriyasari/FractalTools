# Fractal Surface Interpolation demo 

using Makie 
using Random 
using FractalTools
using GeometryBasics 
using LinearAlgebra

# Get interpolation data 
trig  = Triangle(
    Point(BigFloat(-1.), BigFloat(-1)), 
    Point(BigFloat(1.), BigFloat(-1.)), 
    Point(BigFloat(0.), BigFloat(2.))
)
f(x, y) = x^2 + y^2 + 1
pnts3d = getdata(f, trig, 100)

# Construct test points 
tpnts3d = getdata(f, trig, 100)
tpnts = project(tpnts3d)

α = Node(0.001)
niter = Node(1) 
err = lift(niter, α) do n, α 
    itp = interpolate(pnts3d, α=α, maxiters=n)
    fval = map(pnt -> f(pnt...), tpnts)  
    ival = map(pnt -> itp(pnt...), tpnts) 
    val = abs.(fval - ival) ./ abs.(fval) * 100
    val
end 

fig = Figure(title="Fractal Interpolation") 
ax = Axis(fig[1,1], xlabel="npoints", ylabel="Relative Error (%)", title="Fractal Surface Interpolation Error Demo") 
plt = lines!(ax, err) 

gl = GridLayout()
gl2 = GridLayout()
gl[1, 1] = slider = Slider(fig, range=0.001 : 0.0001 : 0.01, startvalue=0.001, tellwidth=false)
gl2[1, 1] = button1 = Button(fig, label="Increment", tellwidth=false)
gl2[1, 2] = button2 = Button(fig, label="Decrement", tellwidth=false)
gl[2, 1] = gl2
gl[1, 2] = label1 = Label(fig, "alpha = $(slider.value[])")
gl[2, 2] = label2 = Label(fig, "niter = $(niter[])")
fig[2, 1] = gl

on(button1.clicks) do n 
    niter[] += 1
    label2.text[] = "niter = " * string(niter[])
    autolimits!(ax)
end 

on(button2.clicks) do n 
    niter[] -= 1
    label2.text[] = "niter = " * string(niter[])
    autolimits!(ax)
end 

on(slider.value) do val 
    α[] = val
    label1.text[] = "alpha = " * string(val)
    autolimits!(ax)
end

display(fig)

# # For video recording.
# fps = 60
# record(fig.scene, "demo.mp4"; framerate = fps) do io
#     for i = 1 : 50 * fps
#         sleep(1/fps)
#         recordframe!(io)
#     end
# end
