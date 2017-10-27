module Minpack2Interface
  import Base.show
  import Base.showerror

  export MINPACK2Info, MINPACK2Exception, GQTWorkspace, throw_if_error

  struct MINPACK2Info
    code::Int
  end

  struct MINPACK2Exception
    info::MINPACK2Info
  end

  include("gqt_workspace.jl")

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
      throw(MINPACK2Exception(v))
    end
  end
end
