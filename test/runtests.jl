include("../minpack2_interface.jl")

using .Minpack2Interface
using Base.Test

# must be const for `ccall`
const gqtpar_lib_path = joinpath(dirname(@__FILE__), "gqtpar.$(Libdl.dlext)")
run(`gfortran -shared -fPIC -lblas -llapack -o$(gqtpar_lib_path) destsv.f dgqt.f`)
@test isfile(gqtpar_lib_path)

include("estsv.jl")
include("gqt.jl")
