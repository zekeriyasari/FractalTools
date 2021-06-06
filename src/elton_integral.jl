using FractalTools
using LinearAlgebra

export elton_integral

function elton_integral(func, ind_func, points::Channel, func_params=tuple(), ind_params=tuple(), ϵ::Real=1e-10, max_iter=1e9, chunk_size = 1024)
    number_of_points = 0
    old_number_of_points = 0
    result = 0
    Δ = 1e15
    # measure = prod(b - a)
    # while Δ >= ϵ && number_of_points ≤ max_iter
    while number_of_points ≤ max_iter
        old_number_of_points = number_of_points
        old_result = result
        number_of_indicator = 0
        set_of_points = []
        indicator = []
        while number_of_indicator < chunk_size 
            set_of_points = take!(points)
            indicator = ind_func.(set_of_points, ind_params...)
            number_of_indicator = sum(indicator .> 0)
        end
        filtered_set_of_points = set_of_points[indicator.>0] .* indicator[indicator .> 0]
        number_of_points += number_of_indicator
        val =  func.(filtered_set_of_points, func_params...) 
        result = (1 / (number_of_points +1)) * ((old_number_of_points + 1) * old_result + sum(val))
        Δ = abs(result - old_result) 
    end
    return  result , Δ, Δ > ϵ && number_of_points > max_iter, number_of_points
end




# An example of one dimensional ifs
# A1 = reshape([1/2],1,1)
# b1 = [0]
# A2 = reshape([1/2],1,1)
# b2 = [1/2]

# w1 = Transformation(A1,b1)
# w2 = Transformation(A2,b2)

# w = [w1, w2]
# ifs = IFS(w)
# generator = randalg_sequential_generator(ifs.ws, ifs.probs)

# elton_integral(dipole_pocklington, is_in_Ball, generator,(0.9,0.005,2π), (0.9, 0.05))

