
function foo(ch; xi=0, dx=1, xf=1)
    for _x in xi : dx : xf 
        put!(ch, _x)
    end 
end 

ch = Channel(0)
task = @async foo(ch, xf = 10)
bind(ch, task)

for i in 1 : 9 
    @show take!(ch) 
end 
task
