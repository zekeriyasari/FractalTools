using FractalTools

# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [-1/8]
A2 = reshape([1/2],1,1)
b2 = [1/8]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)


take!(ifs.generator)