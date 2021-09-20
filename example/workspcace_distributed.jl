# This file is for parallelizatin in task-level 

using Distributed 

# Load process 
numprocs = nprocs() 
numcores = Sys.CPU_THREADS
numprocs == numcores - 1 || addprocs(numcores - 1 - numprocs)

# Define worker function 
@everywhere begin 
   
    function workerfunc(f, chunkchannel, resultchannel) 
        while true 
            chunk = take!(chunkchannel) 
            all(chunk .=== NaN) && break 
            f_val = map(x -> f(x...), chunk)
            res = sum(f_val) 
            put!(resultchannel, (myid(), res)) 
        end 
    end 
    f(x, y) = x^2 + y^2
end 

# Launch communcition chnannel 
chunkchannel = RemoteChannel(() -> Channel(Inf))
resultchannel = RemoteChannel(() -> Channel(Inf))

# Launch worker tasks 
for worker in workers()
    @spawnat worker workerfunc(f, chunkchannel, resultchannel)
end 

# Define attractor 
getatr(dim=2, chunksize=10) = Channel(0) do ch 
    while true
        chunk = [rand(dim) for i in 1 : chunksize]
        put!(ch, chunk)
    end
end 

# Construct attractor 
atr = getatr()

# Feed workers 
for i in 1 : 100 
    put!(chunkchannel, take!(atr))
    @show i, take!(resultchannel)
end 

