using FractalTools
using Debugger

function antenna_pocklington(zp, z, ra, k=2π)
    R = sqrt.(ra^2 .+ (z .- zp).^2)
    exp.(-im*k .*R ) .*  ((1 .+ im*k .*R) .* (2 .* R.^2 .- 3*ra^2) .+ (k*ra .*R).^2) ./ (4π .*R.^5)
end


function indicator_measure_1D(x, x0, ϵ)
    return ϵ * is_in_Ball(x,x0,ϵ) 
end


# An example of one dimensional ifs
A1 = reshape([1/2],1,1)
b1 = [-1/8]
A2 = reshape([1/2],1,1)
b2 = [1/8]

w1 = Transformation(A1,b1)
w2 = Transformation(A2,b2)

w = [w1, w2]
ifs = IFS(w)
generator = randalg_sequential_generator(ifs.ws, ifs.probs)

elton_integral(antenna_pocklington, indicator_measure_1D, generator, (-0.2, 0.005), (-0.25, -0.15), 1e-7, 2e7)

elton_integral(antenna_pocklington, indicator_measure_1D, generator, (0.1, 0.005), (0.2, 0.05), 1e-7, 2e7)

