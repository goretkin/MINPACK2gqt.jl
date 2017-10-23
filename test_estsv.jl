include("minpack2.jl")

r = rand(300,300)
R = full(UpperTriangular(r))

(ldr, n) = size(r); z = Array{Float64}(n);
svmin = Ref{Float64}()

Profile.clear_malloc_data()

@time Minpack2.estsv(n, r, ldr, svmin, z)
@time Minpack2.estsv(n, r, ldr, svmin, z)
@profile Minpack2.estsv(n, r, ldr, svmin, z)


@time _, S, V = svd(R)
@time _, S, V = svd(R)
svmin_j = S[end]
z_j = -V[:,end]

@show svmin[]
@show svmin_j
@show z
@show z_j

include("ell_lib.jl")

bad = []
for n = 1:100
  for rc = 1:n
    A = 5 * (rand(n,n,3) - .5)
    A[rc,:,2] = 0.0
    A[:,rc,3] = 0.0

    for (M,info) in zip([A[:,:,i] for i=1:3], ["rand", "-row", "-col"])
      v1, sv1 = ELL_LIB.estsv(M)
      v2, sv2 = Minpack2.estsv(M)

      ve = (v1-v2)
      ve_flipped = (v1+v2) # sometimes the vector is negated (which may be unacceptable) but let's know that it's jsut a sign flip not totally wrong
      vemax = maximum(ve)
      vemax_overall = min(vemax, maximum(ve_flipped))
      sve = sv1 - sv2
      if vemax_overall > 1e3eps()
        println("Quite bad!!")
        @show vemax_overall
      end
      if !(vemax <= 1e3*eps()) || !(sve <= 1e3*eps())
        println()
        push!(bad, M)
        @show info
        @show n
        @show vemax
        @show sve
        #@show M
      end
    end
  end
end
