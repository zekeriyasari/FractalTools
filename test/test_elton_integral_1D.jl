using FractalTools
using Debugger

function deneme(x,k)
    # k * x
    1
end

function deneme2(x,w)
    sin(w*x)
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
generator = ifs.generator
# generator = randalg_sequential_generator(ifs.ws, ifs.probs)

function indicator_measure_1D(x, x0, ϵ)
    return 2*ϵ * is_in_Ball(x,x0,ϵ) 
end


elton_integral_1D(deneme, indicator_measure_1D, generator, (1,), (0.5, 0.25), 1e-11, 1e7)

elton_integral_1D(deneme2, indicator_measure_1D, generator, (2,), (0.6,0.2), 1e-11, 1e7)



# elton_integral(dipole_pocklington, is_in_Ball, generator,(0.9,0.005,2π), (0.9, 0.05))


