using FractalTools

# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [0]
A2 = reshape([1/2],1,1)
b2 = [1/2]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs1 = IFS(w)
ifs1.ws
ifs1.probs

# An example of two dimensional ifs

A = [1/2 0;0 1/2]
b1 = [0; 0] 
b2 = [0; 1/2] 
b3 = [1/2; 1/2] 
b4 = [1/2; 0] 

w1 = Transformation(A,b1)
w2 = Transformation(A,b2)
w3 = Transformation(A,b3)
w4 = Transformation(A,b4)

w = [w1, w2, w3, w4]
ifs2 = IFS(w)
ifs2.ws
ifs2.probs

# An example of different probabilities

p = [0.2, 0.4, 0.1, 0.3]

ifs3 = IFS(w,p)
ifs3.ws
ifs3.probs

