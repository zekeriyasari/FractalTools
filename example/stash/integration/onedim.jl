using FractalTools
using Plots; plotlyjs()
using Cubature 

# Definitions 
ti, ts, tf = 0, 0.1, 1.
t = ti : ts : tf 
f(t) = sin(2π * t) 
g(t) = cos(2π * t)
# f(t) = t * (1 - t)

# Construct interpolation 
interp = interpolate(t, f.(t), d=0.01)

# Construct hidden interpolation 
hiddeninterp = hiddenfif(t, f.(t), g.(t)) 

# Quadratute integration
hquad(tk) = hquadrature(f, ti, tk, abstol=1e-9)[1]

# Fractal integration 
hfrac(tk) = integrate(ti:ts:tk, f.(ti:ts:tk), d=0.01)

# Fractal hidden integration 
hhiddenfrac(tk) = integrate(ti:ts:tk, f.(ti:ts:tk), g.(ti:ts:tk))

# Plots 
td = ts : 0.01 * ts : tf 
plot(td, f.(td), label="f")
plot!(td, hquad.(td), label="intfquad")
plot!(td, hfrac.(td), label="intffrac")
plot!(td, hhiddenfrac.(td), label="intfhiddenfrac")
plot!(td, getindex.(hiddeninterp.(td), 1), label="fhiddenfrac")
scatter!(t, f.(t), marker=(:circle, 2), label="data")
