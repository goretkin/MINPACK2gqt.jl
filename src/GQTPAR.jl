module GQTPAR
  export GQTWorkspace, solve!
  include("fortran_help.jl")
  using .FortranHelp
  #import FortranHelp: @aref

  include("minpack2_interface.jl")
  using .Minpack2Interface

  include("estsv.jl")
  include("gqt_workspace.jl")
  include("gqt.jl")

  function estsv(R)
    n = size(R, 1)
    z = Array{Float64}(n)
    return estsv!(R, z)
  end

  function estsv!(R, z)
    svmin = Ref{Float64}()  # will heap-allocate on 0.6
    return estsv!(R, z, svmin)
  end

  function estsv!(R, z, svmin)
    n = size(R, 1)
    ldr = stride(R,2)
    estsv(n, R, ldr, svmin, z)
    return (z, svmin[])
  end

end
