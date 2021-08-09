using FractalTools 
import GLMakie 

# Dataset 
f(x, y) = x^2 + y^2
numinterior = 100
numboundary = 10 
Ω = uniformdomain(3) 
pts = getdata(f, Ω, numinterior, numboundary)

msh = tomesh(pts, Ω)

# Plots 
fig, ax, plt = GLMakie.wireframe(msh) 
GLMakie.wireframe!(project(msh, 1))
