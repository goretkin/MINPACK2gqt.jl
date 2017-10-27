function solve!(ws::GQTWorkspace, delta, itmax, atol, rtol)
  @assert 1 == stride(ws.a, 1)
  n = size(ws.a, 1)
  lda = stride(ws.a, 2)

  info_ = Ref{Cint}()
  f_ = Ref{Float64}()

  ccall(
    ("dgqt_", gqtpar_lib_path),
    Void,
    (Ref{Cint}, Ptr{Float64}, Ref{Cint}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Cint}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Cint},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
    n, ws.a, lda, ws.b, delta, rtol, atol,
    itmax, ws.par, f_, ws.x, info_,
    ws.z, ws.wa1, ws.wa2
  )
  ws.info = MINPACK2Info(info_[])
  return (ws, f_[])
end
