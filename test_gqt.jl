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
par = Ref{Float64}()
f = Ref{Float64}()
x = Array{Float64}(n)
info = Ref{Int}()
z = Array{Float64}(n)
wa1 = Array{Float64}(n)
wa2 = Array{Float64}(n)

@show a

Minpack2.gqt(n, a, lda, b, delta, rtol ,atol, itmax,
  par, f, x, info, z, wa1, wa2)
