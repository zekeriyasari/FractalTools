using FractalTools
using LinearAlgebra

export elton_integral_1D, mean
#ToDo: Assume that indicator is 1 for one version
# ToDo : Elton integral function and norm function must run vector.
function mean(x)
    sum(x) / length(x)
end

function elton_integral_1D(func, ind_func, points::Channel, func_params=tuple(), ind_params=tuple(), ϵ::Real=1e-10, max_iter=1e6, chunk_size = 10)
    number_of_points = 0
    old_number_of_points = 0
    result = 0
    Δ = 1e15 * ones(chunk_size)
    k = 0
    # measure = prod(b - a)
    while mean(Δ) >= ϵ && number_of_points ≤ max_iter
    # while number_of_points ≤ max_iter
        old_number_of_points = number_of_points
        old_result = result
        number_of_indicator = 0
        set_of_points = []
        indicator = []
        while number_of_indicator < (chunk_size / 2)
            set_of_points = take!(points)
            indicator = ind_func.(set_of_points, ind_params...)
            number_of_indicator = sum(indicator .> 0)
        end
        filtered_set_of_points = set_of_points[indicator.>0] 
        number_of_points += number_of_indicator
        val =  func.(filtered_set_of_points, func_params...) .* indicator[indicator .> 0]
        result = (1 / (number_of_points +1)) * ((old_number_of_points + 1) * old_result + sum(val))
        Δ[mod(k,chunk_size)+1] = abs(result - old_result) 
        k += 1
    end
    return  result , mean(Δ), mean(Δ) >= ϵ && number_of_points > max_iter, number_of_points
end


function elton_integral_1D(func, points::Channel, func_params=tuple(), ϵ::Real=1e-10, max_iter=1e6, chunk_size = 10)
    number_of_points = 0
    old_number_of_points = 0
    result = 0
    Δ = 1e15 * ones(chunk_size)
    k = 0
    # measure = prod(b - a)
    while mean(Δ) >= ϵ && number_of_points ≤ max_iter
    # while number_of_points ≤ max_iter
        old_number_of_points = number_of_points
        old_result = result
        set_of_points = take!(points)
        val =  func.(set_of_points, func_params...) 
        number_of_points += chunk_size
        result = (1 / (number_of_points +1)) * ((old_number_of_points + 1) * old_result + sum(val))
        Δ[mod(k,chunk_size)+1] = abs(result - old_result) 
        k += 1
    end
    return  result , mean(Δ), mean(Δ) >= ϵ && number_of_points > max_iter, number_of_points
end






