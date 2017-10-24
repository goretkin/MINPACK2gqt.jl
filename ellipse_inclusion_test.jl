include("EllipsoidQueries.jl")
include("ell_lib.jl")
# no random ellipses inside random ellipse happening!

#check inclusion
# generate a random ellipse, and then generate random ellipses within that ellipse.
# L = randn(3,3)
# Σ = L * L'
L = eye(2) + .2 * randn(2,2)
L *= .8
c = randn(2) * 0.0
Σ = L' * L
LL = cholfact(Σ)[:L]

i = 0
j = 0
K = 10

bad = []
Σin = nothing
for j = 1:100
    L = 0.0 * eye(2) + 1.0 * randn(2,2)
    Σt = L' * L
    LLt = cholfact(Σt)[:L]
    ct = randn(2)

    rmaxf = ELL_LIB.ell_pair_cover_query(LL, c, LLt, ct)
    rmaxj = EllipsoidQueries.ellipsoid_pair_cover(Σ, c, Σt, ct)

    if abs(rmaxf - rmaxj) > 1e-8
      push!(bad, (LL, LLt, Σ, Σt, c, ct, rmaxf, rmaxj))
    end
end

f(e) = abs(e[8] - e[7])
(LL, LLt, Σ, Σt, c, ct, rmaxf, rmaxj) = sort(bad, by=f)[99]
