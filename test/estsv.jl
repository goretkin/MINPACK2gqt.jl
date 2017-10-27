function estsv(R)
  @assert 1 == stride(R, 1)
  n = size(R, 1)
  ldr = stride(R, 2)
  svmin = Ref{Float64}()
  z = Array{Float64}(n)

  ccall(
    (:destsv_, gqtpar_lib_path),
    Void,
    (Ref{Cint}, Ptr{Float64}, Ref{Cint}, Ref{Float64}, Ptr{Float64}),
    n, R, ldr, svmin, z
  )
  return (z, svmin[])
end
