module Minpack2
  include("fortran_help.jl")
  using .FortranHelp
  #import FortranHelp: @aref
  include("estsv.jl")
  include("gqt.jl")
end
