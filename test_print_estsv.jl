using JLD
R = load("bad5.jld", "a") * 1e0

ENV["GFORTRAN_UNBUFFERED_ALL"] = "1"

include("ell_lib.jl")
include("minpack2.jl")

if length(ARGS) == 0 || contains(ARGS[1], "F")
println("FORTRAN********************"); flush(STDOUT)
println("Final: ", ELL_LIB.estsv(R)); flush(STDOUT)
end

if length(ARGS) == 0 || contains(ARGS[1], "J")
println("JULIA********************"); flush(STDOUT)
println("Final: ", Minpack2.estsv(R)); flush(STDOUT)
end
