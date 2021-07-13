# This file includes methods to compute Elton integral 

export elton 

# TODO: The method below compute Elton integral after drawing n points in one shot. Update the method such that the points
# are drawn from the attractor stap by step. 

"""
    $SIGNATURES 

Computes Elton integral of `f` over the attracttor `Ω` for `n` points 
"""
function elton(f, Ω , n)
    pts = vcat([take!(Ω) for i in 1 : n]...)
    sum(map(p -> f(p...), pts)) / length(pts) 
end 

function elton(f, Ω, n, chunksize)
    # Update  chunksize of Ω
    Ω.chunksize = chunksize 
    
    # Evaluate integral value in chunks 
    total = 0 
    npts = 0 
    for i in 1 : n 
        chunk = take!(Ω) 
        total += sum(map(pnt -> f(pnt...), chunk))
        npts  += length(chunk)
    end 
    total / npts
end
