# using FractalAntennaTools
# using FractalTools

function elton_integral(func, ind_func, generator::Channel; func_params::NamedTuple=(;), 
                        ind_params::NamedTuple=(;), ϵ::Real=1e-6, max_iter::Int=Int(1e9))
    number = 0
    total = 0
    Δ = 0
    # measure = prod(b - a)
    while  Δ ≥ ϵ || number ≤ max_iter
        old_total = total 
        point = take!(generator)
        # point = next(generator)
        indicator = ind_func.(point; ind_params...)
        point_ind = point .* indicator
        number += sum(indicator)
        val =  func(point_ind; func_params...) 
        total += sum(val) .*  indicator
        Δ = abs(total - old_total)
    end
    elton_val = total ./ (number + 1) 
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
