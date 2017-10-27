module Minpack2
  include("fortran_help.jl")
  using .FortranHelp
  #import FortranHelp: @aref

  include("minpack2_interface.jl")
  using .Minpack2Interface

  include("ell_lib.jl")

  include("estsv.jl")
  include("gqt_workspace.jl")
  include("gqt.jl")

  function estsv(R)
    n = size(R, 1)
    ldr = stride(R,2)
    svmin = Ref{Float64}()
    z = Array{Float64}(n)

    estsv(n, R, ldr, svmin, z)
    return (z, svmin[])
  end
end
