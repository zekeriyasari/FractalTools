# This file includes some testing functions

# ----------------------------------- Univariate test functions ----------------------------------- #

parabol(x) = 1 - (2x - 1)^2

sinusoid(x, f=1) = sin(2π * f * x) 

weierstrass(x, a=0.5, b=3, N=500) = sum([a^n * cos(b^n *  π * x) for n = 1 : N])

cellerier(x, a=2, N=500) = sum([1 / (a^k) * sin(a^k * x) for k = 1 : N])

riemann(x, N=500) = sum([1 / k^2 * sin(k^2 * x) for k = 1 : N])

darboux(x, N=10) = sum([1 / factorial(k) * sin(factorial(k + 1) * x) for k = 1 : N])

function p(x) 
    xi = abs(x) % 2 
    0       ≤   xi  ≤ 1 / 3 ?   0           : 
    1 / 3   ≤   xi  ≤ 2 / 3 ?   3 * xi - 1  : 
    2 / 3   ≤   xi  ≤ 4 / 3 ?   1           : 
    4 / 3   ≤   xi  ≤ 5 / 3 ?   5 - 3 * xi  : 
    5 / 3   ≤   xi  ≤ 2     ?   0           :
    throw(DomainError("$xi is not in the expected domian")) 
end 
schoenberg1(x, N=500) = sum([(1 / 2)^k * p(3^(2k) * x) for k = 0 : N])
schoenberg2(x, N=500) = sum([(1 / 2)^k * p(3^(2k + 1) * x) for k = 0 : N])

function g(x) 
    xi = abs(x) % 4 
    0   ≤   xi  ≤ 2 ?  1 - xi : 
    2   ≤   xi  ≤ 4 ? -3 + xi : 
    throw(DomainError("$xi is not in the expected domian")) 
end 
mccarthy(x, N=500) = sum([(1 / 2)^k * g(4^k * x) for k = 1 : N])

wen(x, N=1000) = prod([1 + (1 / 2)^n * sin(6^n * π * x) for n = 1 : N])

# --------------------------------------- Bivariate test functions --------------------------------------- #

paraboloid(x, y) = x^2 + y^2

crossintray(x, y) = -1e-4 * (abs(sin(x) * sin(y) * exp(abs(100 - sqrt(x^2 + y^2) / π))) + 1)^(0.1)

dropwave(x, y) = -(1 + cos(12 * sqrt(x^2 + y^2))) / (0.5 * (x^2 + y^2) + 2)

griewank(x, y) = (x^2 + y^2) / 4000 - cos(x) * cos(y / sqrt(2)) + 1

holder(x, y) = -abs(sin(x) * cos(y) * exp(abs(1 - sqrt(x^2 + y^2) / π)))

schewefel(x, y) = 837.9658 - x * sin(sqrt(abs(x))) - y * sin(sqrt(abs(y))) 

rastrigin(x, y) = 20 + (x^2 - 10 * cos(2π * x)) + (y^2 - 10 * cos(2π * y))

levy(x, y) = sin(3π * x)^2 + (x - 1)^2 * (1 + sin(3π * y)^2) + (y - 1)^2 * (1 + sin(2π * y)^2)

eggholder(x, y) = -(y + 47) * sin(sqrt(abs(y + x / 2 + 47))) - x * sin(sqrt(abs(x - (y + 47))))

ackley(x, y, a=20, b=0.2, c=2π) = 
    -a * exp(-b * sqrt(1 / 2 * (x^2 + y^2))) - exp(1 / 2 * (cos(c * x) + cos(x * y))) + a + MathConstants.e

schaffer1(x, y) = 0.5 + (sin(x^2 - y^2)^2 - 0.5) / (1 + 0.001 * (x^2 + y^2))^2

schaffer2(x, y) = 0.5 + (cos(sin(x^2 - y^2))^2 - 0.5) / (1 + 0.001 * (x^2 + y^2))^2

schubert(x, y) = sum([i * cos((i + 1) * x + i) for i in 1 : 5]) * sum([i * cos((i + 1) * y + i) for i in 1 : 5])
