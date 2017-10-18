module Minpack2
  import Base.show
  import Base.showerror

  include("fortran_help.jl")
  using .FortranHelp
  #import FortranHelp: @aref
  include("estsv.jl")
  include("gqt.jl")

  struct MINPACK2Info
    info::Int
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
    print("MINPACK2Info: ", v.info, " ")
    println(io, get(gqt_status, v.info, "Unknown Error Code"))
  end

  function showerror(io::IO, ex::MINPACK2Exception)
    print(io, "MINPACK2Exception: ")
    show(io, ex.info)
  end

  function throw_if_error(v::MINPACK2Info)
    if v.info != 1 && info != 2
      throw(Minpack2.MINPACK2Exception(v))
    end
  end

end
