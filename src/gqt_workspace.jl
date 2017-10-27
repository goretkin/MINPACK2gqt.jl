type GQTWorkspace{T}
  a::DenseArray{T,2}
  b::DenseArray{T,1}
  x::DenseArray{T,1}
  z::DenseArray{T,1}
  wa1::DenseArray{T,1}
  wa2::DenseArray{T,1}
  par::T
  info::MINPACK2Info
  function GQTWorkspace{T}(n::Int) where {T}
    new{T}(
      Matrix{T}(n,n),
      Vector{T}(n),
      Vector{T}(n),
      Vector{T}(n),
      Vector{T}(n),
      Vector{T}(n),
      zero(T),
      MINPACK2Info(0))
    end
end

import Base.fill!
function fill!(ws::GQTWorkspace, v)
  fill!(ws.a, v)
  fill!(ws.b, v)
  fill!(ws.x, v)
  fill!(ws.z, v)
  fill!(ws.wa1, v)
  fill!(ws.wa2, v)
  ws.par = v
end
