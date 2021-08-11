# This file includes utility functions 

export combine, randomstring, randommeshpath

"""
    $SIGNATURES 

Combines `vals`.

# Example 
```julia 

julia> foo = [rand(2) for i in 1 : 3]
3-element Vector{Vector{Float64}}:
 [0.07304969216674784, 0.13922697837562814]
 [0.17391711119297937, 0.01724058232513692]
 [0.6244409832793123, 0.9377684091200582]

julia> bar = [rand(1) for i in 1 : 3]
3-element Vector{Vector{Float64}}:
 [0.694956323617361]
 [0.5868716964238885]
 [0.5625412213942036]

julia> combine(foo, bar) 
3-element Vector{Vector{Float64}}:
 [0.07304969216674784, 0.13922697837562814, 0.694956323617361]
 [0.17391711119297937, 0.01724058232513692, 0.5868716964238885]
 [0.6244409832793123, 0.9377684091200582, 0.5625412213942036]
```
"""
combine(vals::AbstractVector...) = [vcat(val...) for val in zip(vals...)]

"""
    $SIGNATURES

Returns a random string whose characters are between `A` and `z` .`n` is the number of characters. 

# Example
```julia 
julia> randomstring(10) 
"EydQuv\\txF"
```
"""
randomstring(n::Int=10) = String(rand('A' : 'z', n))

""" 
    $SIGNATURES 

Returns a random path for `.msh` file. `n` is the number of characters. Default directory for the file is temp directory of
the system 
    
# Example 
```julia 

```
"""
randommeshpath(n::Int=10) = joinpath(tempdir(), randomstring(n) * ".msh")
