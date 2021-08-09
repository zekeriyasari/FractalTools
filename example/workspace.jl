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

pts1 = getdata([[0., 0.], [1., 0], [1., 1.]], numinterior, numboundary)
tess1 = Tessellation(pts1, Ω)
tess2 = Tessellation(pts1)
pts13 = getdata([[0.], [1.]], numinterior, numboundary)
tess3 = Tessellation(pts3)

FractalTools.triangulationtype(tess1) 
FractalTools.triangulationtype(tess2) 
FractalTools.triangulationtype(tess3) 
