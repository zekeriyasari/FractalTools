using FractalTools
using Debugger

ifs = Sierpinski() 
generator = ifs.generator
# generator = randalg_sequential_generator(ifs.ws, ifs.probs)

function indicator_measure(x, x0, ϵ)
    return 1 * is_in_Ball(x,x0,ϵ,norm) 
end

function deneme_2D(x,k)
    1 * ones(size(x,2))
end

elton_integral(deneme_2D, indicator_measure, generator,(2,), ([0,0],1), 1e-10, 1e7)

# elton_integral(deneme_2D, generator,(1,), 1e-8, 1e7)
