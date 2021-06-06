using Base: is_interactive

using FractalTools

function deneme(x,k)
    k * x
end

# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [0]
A2 = reshape([1/2],1,1)
b2 = [1/2]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)
generator = randalg_sequential_generator(ifs.ws, ifs.probs)

function indicator_measure_1D(x, x0, ϵ)
    return ϵ * is_in_Ball(x,x0,ϵ) 
end


elton_integral(deneme, indicator_measure_1D, generator,(1,), (0, 0.5), 1e-11, 2e7)

# ToDo : plot elton_integral vs max_iter

