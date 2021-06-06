using FractalTools

A1 = reshape([1/2],1,1)
b1 = [0]
A2 = reshape([1/2],1,1)
b2 = [1/2]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)

generator = randalg_sequential_generator(ifs.ws, ifs.probs)

take!(generator)
