using GQTPAR.Minpack2Interface
using Base.Test

jp(x) = joinpath(dirname(@__FILE__), x)

# must be const for `ccall`
const gqtpar_lib_path = jp("gqtpar.$(Libdl.dlext)")
run(`gfortran -shared -fPIC -lblas -llapack -o$(gqtpar_lib_path) $(jp("destsv.f")) $(jp("dgqt.f"))`)
@test isfile(gqtpar_lib_path)

include("estsv.jl")
include("gqt.jl")
