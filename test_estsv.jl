include("minpack2.jl")
using Minpack2

r = rand(3,3)
R = full(UpperTriangular(r))

(ldr, n) = size(r); z = Array{Float64}(n);
svmin = 0.0

@time Minpack2.estsv(n, r, ldr, convert(Ptr{Float64}, pointer_from_objref(svmin)), z)

@time _, S, V = svd(R)
svmin_j = S[end]
z_j = -V[:,end]

@show svmin
@show svmin_j
@show z
@show z_j
