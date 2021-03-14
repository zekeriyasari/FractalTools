using FractalTools
using Plots 

# Define function 
f(x) = 4 - (x - 2)^2 + 1

# Interpolation points 
xi, dx, xf = 0, 1, 4 
x = collect(xi : dx : xf) 
y = f.(x) 
d = 0.5 * ones(length(x) - 1)

# Construct interpolation 
interp = interpolate(x, y, d)
ws = interp.ifs.ws

# End points 
xc = [0, 4, 4, 0, 0]
yc = [0, 0, 4, 4, 0]
pc = [xc yc]

for niter in 1 : niters
    xt = collect(xi : 0.1dx : xf)
    plt = plot(xt, f.(xt))
    pc = vcat([mapslices(wt, pc, dims=2) for wt in ws]...)
end


function transform_plot()
    # plt = plot(getindex.(pc, 1), getindex.(pc,2), marker=(:circle, 3))
    # xt = collect(xi : 0.1dx : xf)
    # plot!(xt, f.(xt))
    for wt in ws 
        wp = []
        for (xi, yi) in zip(xc, yc)
            push!(wp, wt([xi, yi]))
        end 
        plot!(getindex.(wp, 1), getindex.(wp,2), marker=(:circle, 3), subplot=1)
    end
    plt
end 

transform_plot()
