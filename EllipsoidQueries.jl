module EllipsoidQueries

include("minpack2.jl")
using .Minpack2

export ellipsoid_extremal_point, ellipsoid_pair_cover

function ellipsoid_extremal_point(E, c, p, sense=1)
  # minimize sense*|x-p| s.t. |x-c|^2_E ≦ 1
  # return (x, info)
  n = size(E, 1)

  U = chol(E) # E = U'*U
  iU = inv(U)
  # y = U * (x-c)
  # x = c + U^-1 * y

  # χ = 1/2 y' * U'^-1 * U^-1 * y + (c-p)' * y
  # min sense*χ
  # s.t. |y| ≦ 1

  a = sense * iU' * iU
  b = sense * (c - p)
  Δ = 1.0

  y = Minpack2.solve_gqt(a, b, Δ, 100, 1e-12, 1e-4)

  x = c + iU * y
  return x
end

function ellipsoid_pair_cover(U1::UpperTriangular, c1, U2::UpperTriangular, c2)
  n = size(U1, 1)
  yU2 = U2 * inv(U1)

  yS2 = yU2' * yU2
  yc2 = U1 * (c2-c1)
  # implicit: yS1 = eye(n)
  # implicit: yc1 = zeros(n)

  x = ellipsoid_extremal_point(yS2, yc2, zeros(n), -1)
  return norm(x) # if < 1, E1 covers E2
end

function ellipsoid_pair_cover(S1, c1, S2, c2)
  # Does E1 cover E2

  # transform space so that E1 is unit ball at origin
  # E1 = {x | |x-c1|^2_S1 ≦ 1}
  # E2 = {x | |x-c2|^2_S2 ≦ 1}
  # y = U1*(x-c1)
  # x = c1 + U1^-1 * y

  # yE1 = {y | |(c1 + U1^-1 * y)-c1|^2_S1 ≦ 1}
  # yE1 = {y | |y|^2 ≦ 1}
  # yE2 = {y | |(c1 + U1^-1 * y)-c2|^2_S2 ≦ 1}
  # yE2 = {y | |(c1-c2) + U1^-1 * y|^2_(U2'*U2) ≦ 1}
  # yE2 = {y | |U1*(c1-c2) + y|^2_(U2'^-1*U2'*U2*U1^-1) ≦ 1}

  U1 = chol(S1)
  U2 = chol(S2)
  return ellipsoid_pair_cover(U1, c1, U2, c2)
end

function ellipsoid_extremal_point_work(E, c, p, sense=1)
  # minimize sense*|x-p| s.t. |x-c|^2_E ≦ 1
  # return (x, info)
  n = size(E, 1)

  U = chol(E) # E = U'*U
  iU = inv(U)
  # y = U * (x-c)
  # x = c + U^-1 * y

  # χ = 1/2 y' * U'^-1 * U^-1 * y + (c-p)' * y
  # min sense*χ
  # s.t. |y| ≦ 1

  a = sense * iU' * iU
  b = sense * (c - p)
  Δ = 1.0

  return (a, b)
end

function ellipsoid_pair_cover_work(U1::UpperTriangular, c1, U2::UpperTriangular, c2)
  n = size(U1, 1)
  yU2 = U2 * inv(U1)

  yS2 = yU2' * yU2
  yc2 = U1 * (c2-c1)
  # implicit: yS1 = eye(n)
  # implicit: yc1 = zeros(n)

  return ellipsoid_extremal_point_work(yS2, yc2, zeros(n), -1)

end

function ellipsoid_pair_cover_work(S1, c1, S2, c2)
  U1 = chol(S1)
  U2 = chol(S2)
  return ellipsoid_pair_cover_work(U1, c1, U2, c2)

end

end
