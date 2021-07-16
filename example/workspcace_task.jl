# This file is for parallelizatin in task-level 

getatr(dim=2, chunksize=10) = Channel(0) do ch 
    while true
        chunk = [rand(dim) for i in 1 : chunksize]
        put!(ch, chunk)
    end
end 

f(x,y) = x^2 + y^2
function worker(f, chunkchannel, resultchannel) 
    while true 
        chunk = take!(chunkchannel) 
        all(chunk .=== NaN) && break 
        f_val = map(x -> f(x...), chunk)
        res = sum(f_val) 
        put!(resultchannel, res) 
    end 
end 

chunkchannel = Channel(0) 
resultchannel = Channel(0) 

task = @async worker(f, chunkchannel, resultchannel)

atr = getatr()

for i in 1 : 10 
    put!(chunkchannel, take!(atr))
    @show i, take!(resultchannel)
end 

