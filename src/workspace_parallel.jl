using Distributed
using SpecialFunctions
using BenchmarkTools
using Base.Threads
import Base.Threads.@spawn


x = range(0,100, length=10000)
results = zeros(length(x))
results .= besselj1.(x)
@btime results .= besselj1.(x)

@btime for i in 1:length(x)
    results[i] = besselj1(x[i])
end

@btime @threads for i in 1:length(x)
    results[i] = besselj1(x[i])
end

function fib(n::Int)
    if n < 2
        return n
    end
    t = fib(n - 2)
    return fib(n - 1) + t
end

function fib_threads(n::Int)
    if n < 2
        return n
    end
    t = @spawn fib_threads(n - 2)
    return fib(n - 1) + fetch(t)
end

function slow_func(x)
    sleep(0.005) #sleep for 5ms
    return x
end

@btime let
    a = @spawn slow_func(2)
    b = @spawn slow_func(4)
    c = @spawn slow_func(42)
    d = @spawn slow_func(12)
    res = fetch(a) .+ fetch(b) .* fetch(c) ./ fetch(d)
end

@btime let
    a = slow_func(2)
    b = slow_func(4)
    c = slow_func(42)
    d = slow_func(12)
    res = a .+ b .* c ./ d
end

@btime let
    x = 1:100
    a = @spawn sin(2)
    b = @spawn sin(4)
    c = @spawn sin(42)
    d = @spawn sin(12)
    res = fetch(a) .+ fetch(b) .* fetch(c) ./ fetch(d)
end

@btime let
    x = 1:100
    a = sin(2)
    b = sin(4)
    c = sin(42)
    d = sin(12)
    res = a .+ b .* c ./ d
end

fetch(@spawn myid())
fetch(@spawnat 3 myid())

@everywhere function my_func(x)
    return x^3*cos(x)
end
fetch(@spawnat 2 my_func(0.42))

@everywhere using SharedArrays
res = SharedArray(zeros(10))
@distributed for x in 1:10
    res[x] = my_func(x)
end
pmap(my_func, 1:10)
