using FractalTools
using Debugger

ifs = Sierpinski() 
generator = ifs.generator
# generator = randalg_sequential_generator(ifs.ws, ifs.probs)

function indicator_measure_2D(x, x0, ϵ)
    return (1/3) * is_in_Ball_2D(x,x0,ϵ) 
end

function deneme_2D(x,k)
    1 * ones(size(x,2))
end

elton_integral_2D(deneme_2D, indicator_measure_2D, generator,(2,), ([0.25,0.25],[0.25,0.25]), 1e-7, 1e7)

# elton_integral(deneme_2D, generator,(1,), 1e-8, 1e7)
