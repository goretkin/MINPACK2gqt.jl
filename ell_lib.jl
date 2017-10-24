module ELL_LIB
# ELL_LIB uses a packed format for triangular matrices

# here is some code, from ell_pair_cover_query.f90 that unpacks gg1 into g1
# g1 will be lower triangular
#=
k  = 0
g1 = 0.d0
g2 = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g1(i,j) = gg1(k)
   end do
end do
=#

function pack!(l::AbstractArray, L::LowerTriangular)
  n = LinAlg.checksquare(L)
  k = 0
  for j=1:n, i=j:n
    k += 1
    l[k] = L[i,j]
  end
  return nothing
end

function pack{T}(L::LowerTriangular{T})
  n = LinAlg.checksquare(L)
  np = n*(n+1)÷2
  l = Array{T,1}(np)
  pack!(l, L)
  return l
end


function ell_pair_cover_query(L1, c1, L2, c2)
  n = LinAlg.checksquare(L1)
  assert(n == LinAlg.checksquare(L2))
  assert((n,) == size(c1))
  assert((n,) == size(c2))

  rmax = Ref{Float64}()
  ccall(
    ("ell_pair_cover_query_", "ELL_LIB_R1/ell.so"),
    Void,
    (Ref{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}),
    n, c1, pack(L1), c2, pack(L2), rmax
  )
  return rmax[]
end

function estsv(R)
  n = size(R, 1)
  ldr = stride(R,2)
  svmin = Ref{Float64}()
  z = Array{Float64}(n)

  ccall(
    ("destsv_", "ELL_LIB_R1/ell.so"),
    Void,
    (Ref{Int}, Ptr{Float64}, Ref{Int}, Ref{Float64}, Ptr{Float64}),
    n, R, ldr, svmin, z
  )
  return (z, svmin[])
end

include("minpack2_interface.jl")
using .Minpack2Interface

function solve!(ws::GQTWorkspace, delta, itmax, atol, rtol)
  n = size(ws.a, 1)
  lda = stride(ws.a, 2)

  info_ = Ref{Int}()
  f_ = Ref{Float64}()

  ccall(
    ("dgqt_", "ELL_LIB_R1/ell.so"),
    Void,
    (Ref{Int}, Ptr{Float64}, Ref{Int}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Int}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Int},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
    n, ws.a, lda, ws.b, delta, rtol, atol,
    itmax, ws.par, f_, ws.x, info_,
    ws.z, ws.wa1, ws.wa2
  )
  ws.info = MINPACK2Info(info_[])
  return (ws, f_[])
end

function solve_gqt{T}(A::Matrix{T}, b::Vector{T}, delta::T, itmax::Int=100, atol::T=5e3*eps(T), rtol::T=5e11*eps(T))
  n = size(A,1)
  ws = GQTWorkspace{T}(n)
  fill!(ws, 9.87654321)
  ws.a = A
  ws.b = b
  solve!(ws, delta, itmax, atol, rtol)
  throw_if_error(ws.info)
  return ws.x
end

end
