A = -(ones(2,2)  + eye(2))
b = [1.0, -1.0]

ENV["GFORTRAN_UNBUFFERED_ALL"] = "1"

include("ell_lib.jl")
include("minpack2.jl")

if length(ARGS) == 0 || contains(ARGS[1], "F")
println("FORTRAN********************"); flush(STDOUT)
println("Final: ", ELL_LIB.solve_gqt(A, b, 1.0)); flush(STDOUT)
end

if length(ARGS) == 0 || contains(ARGS[1], "J")
println("JULIA********************"); flush(STDOUT)
println("Final: ", Minpack2.solve_gqt(A, b, 1.0)); flush(STDOUT)
end
