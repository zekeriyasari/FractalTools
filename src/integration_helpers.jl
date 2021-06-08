# using FractalAntennaTools
# using FractalTools
using LinearAlgebra
export is_in_Ball, heaviside, interval

function is_in_Ball_v2(x::AbstractVector, x0::AbstractVector, ϵ, norm_func = norm, args...)
    return norm_func(x - x0, args...) < ϵ
end
# is_in_Ball_v2([0 1 ; 0  0.1], [0; 0.2], 0.05)

# Check norm.()

# function is_in_Ball(x, x0, ϵ, norm_func = norm, args...)
#     return norm_func.(x - x0, args...) < ϵ
# end

function heaviside(t)
    0.5 .* (sign.(t) .+ 1)
end

function interval(t, a, b)
    heaviside.(t-a) .- heaviside.(t-b)
end 


# function getatrset1D(N_iter)
#     A1 = reshape([1/2],1,1)
#     b1 = [-1/8]
#     A2 = reshape([1/2],1,1)
#     b2 = [1/8]

#     w1 = FractalTools.Transformation(A1,b1)
#     w2 = FractalTools.Transformation(A2,b2)
#     w = [w1, w2]
#     ifs = IFS(w)

#     initset = [[0.1]]
#     # atr = random_sequential_generator(...)
#     # atr = ifs_sequential_generator(...)
#     atr = attractor(ifs, initset, alg=RandAlg(), numiter=N_iter, allocated=false) 
# end 

# function getatrset2D(N_iter)
#     ifs = Square() 
#     initset = [rand(2)]
#     atr2D = attractor(ifs, initset, alg=RandAlg(), numiter=Int(N_iter), allocated=false)
# end


# function deneme1D(point, c)
#     2 * getindex(point,1) - c
# end 

# function deneme2D(point, c)
#     getindex(point,1) .+ getindex(point,2) .- c 
# end

# function antenna_pocklington(zp, z, ra, k=2π)
#     R = sqrt.(ra^2 .+ (z .- zp).^2)
#     exp.(-im*k .*R ) .*  ((1 .+ im*k .*R) .* (2 .* R.^2 .- 3*ra^2) .+ (k*ra .*R).^2) ./ (4π .*R.^5)
# end


# # function elton_general(func, ind_func, atr, a, b, args...)
# #     number = 0
# #     total = 0
# #     measure = prod(b-a)
# #     @show measure
# #     for point in atr.set
# #         @show point
# #         indicator  = ind_func.(point, a, b)
# #         if all(indicator .> 0 )
# #             number += 1
# #             @show point
# #             @show number
# #             val =  func(point, args...)
# #             @show val
# #             total += sum(val)
# #             @show total 
# #         end
# #     end
# #     elton = total ./ (number + 1) * measure
# # end

# # atr = getatrset1D(10^1)
# # elton_general(deneme1D, interval, atr, 0, 0.1, 0)

# # atr = getatrset2D(1e1)
# # elton_general(deneme2D, interval, atr, [50,30], [100,100], 0)


# # function elton_general(func, ind_func, atr, a, b, N_iter, args...)
# #     number = 0
# #     total = 0
# #     measure = prod(b-a)
# #     for i = 1: N_iter
# #         point = take!(atr.set)
# #         indicator  = ind_func.(point, a, b)
# #         if all(indicator .> 0 )
# #             number += 1
# #             val =  func(point, args...)
# #             total += sum(val)
# #         end
# #     end
# #     elton_val = total ./ (number + 1) * measure
# # end


