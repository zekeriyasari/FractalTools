using FractalTools
using Plots

function example_attractor_of_ifs(N_iter)
    # An example of one dimensional ifs
    A1 = reshape([1/2],1,1)
    b1 = [0]
    A2 = reshape([1/2],1,1)
    b2 = [1/2]

    w1 = Transformation(A1,b1)
    w2 = Transformation(A2,b2)

    w = [w1, w2]
    ifs = IFS(w)
    initset = [[0.1]]
    atr = attractor(ifs, initset; alg=RandAlg(), numiter=N_iter,allocated=true) 
    attrset = vcat(atr.set...)
    t = ones(length(attrset))
    return atr, attrset, t
end

atr, attrset, t = example_attractor_of_ifs(10)
p1 = scatter(t, attrset, title="Attractor of the IFS", label="Attractor Points", markersize = 4)

atr, attrset, t = example_attractor_of_ifs(100)
p2 = scatter(t, attrset, title="Attractor of the IFS", label="Attractor Points", markersize = 4)

atr, attrset, t = example_attractor_of_ifs(10000)
p3 = scatter(t, attrset, title="Attractor of the IFS", label="Attractor Points", markersize = 4)


plot(p1, p2, p3, layout=(1,3))
xlabel!("x")
ylabel!("z")
xlims!((0.9,1.1))
# savefig("myplot.png")


# An example of atr with Channel
atr = attractor(ifs1, initset; alg=RandAlg(), numiter=100, allocated=false) 
collect(atr.set)

# An example of atr with randalg_sequential_generator
generator = randalg_sequential_generator(ifs1.ws, ifs1.probs)
take!(generator)


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

