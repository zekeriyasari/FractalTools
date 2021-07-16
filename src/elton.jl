# This file includes methods to compute Elton integral 

export elton_integral, indicator_measure

function elton_integral(func, ind_func, points::Attractor, 
    func_params=tuple(), ind_params=tuple(), ϵ::Real=1e-10, max_iter=1e6, windowsize=10)
    number_of_points = 0
    old_number_of_points = 0
    result = 0
    Δ = 1e15 * ones(windowsize)
    chunksize = points.chunksize
    k = 0
    while mean(Δ) >= ϵ && number_of_points ≤ max_iter
        old_number_of_points = number_of_points
        old_result = result
        number_of_indicator = 0
        set_of_points = []
        indicator = []
        index_indicator = []
        while number_of_indicator < (chunksize / 2)
            set_of_points = take!(points)
            indicator = ind_func(set_of_points, ind_params...)
            index_indicator = indicator .> 0
            number_of_indicator = sum(index_indicator)
        end
        filtered_set_of_points = set_of_points[index_indicator] 
        number_of_points += number_of_indicator
        val =  func(filtered_set_of_points, func_params...) .* indicator[index_indicator]
        result = (1 / (number_of_points +1)) * ((old_number_of_points + 1) * old_result + sum(val))
        Δ[k % chunksize + 1] = abs(result - old_result) 
        k += 1
    end
    return result, mean(Δ), mean(Δ) >= ϵ && number_of_points > max_iter, number_of_points
end


function elton_integral(func, points::Attractor, func_params=tuple(), ϵ::Real=1e-10, 
    max_iter=1e6, windowsize=10)
    number_of_points = 0
    old_number_of_points = 0
    result = 0
    Δ = 1e15 * ones(windowsize)
    chunksize = points.chunksize
    k = 0
    while mean(Δ) >= ϵ && number_of_points ≤ max_iter
        old_number_of_points = number_of_points
        old_result = result
        set_of_points = take!(points)
        val =  func(set_of_points, func_params...) 
        number_of_points += chunksize
        result = (1 / (number_of_points +1)) * ((old_number_of_points + 1) * old_result + sum(val))
        Δ[k % chunksize + 1] = abs(result - old_result) 
        k += 1
    end
    return result, mean(Δ), mean(Δ) >= ϵ && number_of_points > max_iter, number_of_points
end

mean(x) = sum(x) / length(x)

function is_in_Ball(x, x0, ϵ, norm_func=norm, args...)
    norm_func.(map(xi -> xi - x0, x)) .≤ ϵ
    # norm_func.(eachcol(x .- x0), args...) .<= ϵ
end

# function is_in_Ball(x, x0, ϵ, norm_func=norm, args...)
#     return norm_func.(x - x0, args...) < ϵ
# end

# Indicator measure 
indicator_measure(x, x0, ϵ) = 2ϵ * is_in_Ball(x, x0, ϵ) 

# """
#     $SIGNATURES 

# Computes Elton integral of `f` over the attracttor `Ω` for `n` points 
# """
# function elton(f, Ω , n)
#     pts = vcat([take!(Ω) for i in 1 : n]...)
#     sum(map(p -> f(p...), pts)) / length(pts) 
# end 

# function elton(f, Ω, n, chunksize)
#     # Update  chunksize of Ω
#     Ω.chunksize = chunksize 
    
#     # Evaluate integral value in chunks 
#     total = 0 
#     npts = 0 
#     for i in 1 : n 
#         chunk = take!(Ω) 
#         total += sum(map(pnt -> f(pnt...), chunk))
#         npts  += length(chunk)
#     end 
#     total / npts
# end
