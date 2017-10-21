module Minpack2
  import Base.show
  import Base.showerror

  include("fortran_help.jl")
  using .FortranHelp
  #import FortranHelp: @aref

  struct MINPACK2Info
    code::Int
  end

  struct MINPACK2Exception
    info::MINPACK2Info
  end

  gqt_status = Dict(
    1=>"The function value f(x) has the relative accuracy specified by rtol.",
    2=>"The function value f(x) has the absolute accuracy specified by atol.",
    3=>"Rounding errors prevent further progress. On exit x is the best available approximation.",
    4=>"Failure to converge after itmax iterations. On exit x is the best available approximation."
    )

  function show(io::IO, v::MINPACK2Info)
    print("MINPACK2Info: ", v.code, " ")
    println(io, get(gqt_status, v.code, "Unknown Error Code"))
  end

  function showerror(io::IO, ex::MINPACK2Exception)
    print(io, "MINPACK2Exception: ")
    show(io, ex.info)
  end

  function throw_if_error(v::MINPACK2Info)
    if v.code != 1 && v.code != 2
      throw(Minpack2.MINPACK2Exception(v))
    end
  end

  include("estsv.jl")
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
