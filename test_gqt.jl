include("minpack2.jl")
using Minpack2

n = 5
center = ones(n)
center[:] = randn(n)
#fill!(center, 0.0)
#center[1] = 1.0
a = eye(n)
a[:,:] = diagm(rand(n))
lda = n
b = zeros(n)
b[:] = - a * center
delta = 1.0
rtol = 1e-4
atol = 1e-6
itmax = 1000
par = 0.0
f = 0.0
x = Array{Float64}(n)
info = 6987234598757 # random integer not interned
z = Array{Float64}(n)
wa1 = Array{Float64}(n)
wa2 = Array{Float64}(n)

@show a

Minpack2.gqt(n, a, lda, b, delta, rtol ,atol, itmax,
  Ptr{Float64}(pointer_from_objref(par)),
  Ptr{Float64}(pointer_from_objref(f)),
  x, Ptr{Int}(pointer_from_objref(info)), z, wa1, wa2)
