include("src/integration_helpers.jl")

using Test

# @testset ".jl" begin
    x0 = 2 * rand(10) .- 1
    x = 2 * rand(10) .- 1
    ϵ = rand(10)
    for _x0 in x0 
        for _x in x
            for _ϵ in ϵ
                result = _x0 - _ϵ < _x < _x0 + _ϵ
                 @show is_in_Ball(_x, _x0, _ϵ) == result
            end
        end
    end
# end

