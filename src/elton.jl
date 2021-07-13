# This file includes methods to compute Elton integral 

export elton 

"""
    $SIGNATURES 

Computes Elton integral of `f` over the attracttor `Ω` for `n` points 
"""
function elton(f, Ω , n)
    pts = vcat([take!(Ω) for i in 1 : n]...)
    sum(map(p -> f(p...), pts)) / length(pts) 
end 
