
function f1(x1, y1=1, z1=2; t1=3, u1=5)
    @show x1, y1, z1, t1, u1
    nothing 
end 

function f2(x2, y2=1, z2=2; t2=3, u2=5)
    @show x2, y2, z2, t2, u2
    nothing 
end 

function f3(x3, f1args = tuple(), f2args = tuple(); f1kwargs = (;), f2kwargs=(;))
    @info "Calling f1"
    f1(x3, f1args...; f1kwargs...)
    @info "Calling f2"
    f2(x3, f2args...; f2kwargs...)
end

f3(5) 
f3(5, (2, 3)) 
f3(5, (4, 10), (6, 0))
f3(5, (4, 10), (6, 0), f1kwargs=(t1 = 5,))
f3(5, (4, 10), (6, 0), f1kwargs=(t1 = 5, u1 = 10))
f3(5, (4, 10), (6, 0), f1kwargs=(t1 = 5, u1 = 10), f2kwargs=(t2 = 5,))
f3(5, (4, 10), (6, 0), f1kwargs=(t1 = 5, u1 = 10), f2kwargs=(t2 = 5, u2 = 10))