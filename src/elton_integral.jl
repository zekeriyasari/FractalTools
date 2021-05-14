# using FractalAntennaTools
# using FractalTools

# d \mü 
# TODO if there not exist points channel, use a channel generated random points. 
# TODO ::AbstractFunction vs.
function elton_integral(func, ind_func, points::Channel; func_params::NamedTuple=(;), ind_params::NamedTuple=(;), ϵ::Real=1e-6, max_iter::Int=Int(1e9))
    number = 0
    total = 0
    Δ = 0
    # measure = prod(b - a)
    while Δ ≥ ϵ || number ≤ max_iter
        old_total = total 
        set_of_points = take!(points)
        indicator = ind_func.(set_of_points, ind_params...)
        filtered_set_of_points = set_of_points .* indicator
        number += sum(indicator)
        val =  func.(filtered_set_of_points, func_params...) 
        total += sum(val)
        Δ = abs(total / (number + 1) - old_total / (number + 1 + chunk_size))
    end
    # elton_val = total ./ (number + 1) 
    # elton_val = total ./ (number + 1) * measure
    return number > max_iter & Δ > ϵ,  elton_val, Δ
end


# function elton_general(func, ind_func, atr, a, b, args...)
#     number = 0
#     total = 0
#     measure = prod(b-a)
#     for point in atr.set
#         indicator  = ind_func.(point, a, b)
#         if all(indicator .> 0 )
#             number += 1
#             val =  func(point, args...)
#             total += sum(val)
#         end
#     end
#     elton_val = total ./ (number + 1) * measure
# end
