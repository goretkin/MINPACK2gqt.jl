using GQTPAR
using Base.Test

module Foreign
  using Base.Test
  jp(x) = joinpath(dirname(@__FILE__), x)

  # must be const for `ccall`
  const gqtpar_lib_path = jp("gqtpar.$(Libdl.dlext)")
  run(`gfortran -shared -fPIC -lblas -llapack -o$(gqtpar_lib_path) $(jp("destsv.f")) $(jp("dgqt.f"))`)
  @test isfile(gqtpar_lib_path)

  using GQTPAR
  include("estsv.jl")
  include("gqt.jl")
end

srand(0)
for n = 1:10:200
  R = zeros(n,n)
  for i=1:100
    randn!(R) # only using upper triangular part
    svec_j, sval_j = GQTPAR.estsv(R)
    svec_f, sval_f = Foreign.estsv(R)
    @test (sval_j - sval_f) <= eps()

    # singular vectors may be off by a sign
    err = minimum(maximum(abs.(svec_j - c*svec_f)) for c in [-1, 1])
    @test err <= 1e5*eps()  #hmph
  end
end

srand(0)
atol = 1e-12
rtol = 1e-6
itmax = 50
for n = 1:10:200
  ws_j = GQTPAR.GQTWorkspace{Float64}(n)
  ws_f = GQTPAR.GQTWorkspace{Float64}(n)

  for i=1:10
    fill!(ws_j, 0.0)
    fill!(ws_f, 0.0)

    randn!(ws_j.a)
    randn!(ws_j.b)
    delta = 5 * rand()
    ws_j.par = 0.0

    copy!(ws_f.a, ws_j.a)
    copy!(ws_f.b, ws_j.b)
    ws_f.par = ws_j.par

    GQTPAR.solve!(ws_j,  delta, itmax, atol, rtol)
    Foreign.solve!(ws_f,  delta, itmax, atol, rtol)

    @test ws_j.x ≈ ws_f.x
    @test abs(ws_j.z[1]- ws_f.z[1]) <= 1.0 # iterations taken
    @test ws_j.info == ws_f.info
  end
end
