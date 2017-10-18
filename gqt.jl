function gqt(n::Integer,a::DenseArray{Float64,2},lda::Integer,b::DenseArray{Float64,1},delta::Float64,rtol::Float64,atol::Float64,itmax::Integer,par::Float64,f::Float64,x::DenseArray{Float64,1},info::Integer,z::DenseArray{Float64,1},wa1::DenseArray{Float64,1},wa2::DenseArray{Float64,1})

#
# Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
# see CopyrightMINPACK.txt
#
#     ***********
#
#     Subroutine dgqt
#
#     Given an n by n symmetric matrix A, an n-vector b, and a
#     positive number delta, this subroutine determines a vector
#     x which approximately minimizes the quadratic function
#
#           f(x) = (1/2)*x'*A*x + b'*x
#
#     subject to the Euclidean norm constraint
#
#           norm(x) <= delta.
#
#     This subroutine computes an approximation x and a Lagrange
#     multiplier par such that either par is zero and
#
#            norm(x) <= (1+rtol)*delta,
#
#     or par is positive and
#
#            abs(norm(x) - delta) <= rtol*delta.
#
#     If xsol is the solution to the problem, the approximation x
#     satisfies
#
#            f(x) <= ((1 - rtol)**2)*f(xsol)
#
#     The subroutine statement is
#
#       subroutine dgqt(n,a,lda,b,delta,rtol,atol,itmax,
#                        par,f,x,info,z,wa1,wa2)
#
#     where
#
#       n is an integer variable.
#         On entry n is the order of A.
#         On exit n is unchanged.
#
#       a is a double precision array of dimension (lda,n).
#         On entry the full upper triangle of a must contain the
#            full upper triangle of the symmetric matrix A.
#         On exit the array contains the matrix A.
#
#       lda is an integer variable.
#         On entry lda is the leading dimension of the array a.
#         On exit lda is unchanged.
#
#       b is an double precision array of dimension n.
#         On entry b specifies the linear term in the quadratic.
#         On exit b is unchanged.
#
#       delta is a double precision variable.
#         On entry delta is a bound on the Euclidean norm of x.
#         On exit delta is unchanged.
#
#       rtol is a double precision variable.
#         On entry rtol is the relative accuracy desired in the
#            solution. Convergence occurs if
#
#              f(x) <= ((1 - rtol)**2)*f(xsol)
#
#         On exit rtol is unchanged.
#
#       atol is a double precision variable.
#         On entry atol is the absolute accuracy desired in the
#            solution. Convergence occurs when
#
#              norm(x) <= (1 + rtol)*delta
#
#              max(-f(x),-f(xsol)) <= atol
#
#         On exit atol is unchanged.
#
#       itmax is an integer variable.
#         On entry itmax specifies the maximum number of iterations.
#         On exit itmax is unchanged.
#
#       par is a double precision variable.
#         On entry par is an initial estimate of the Lagrange
#            multiplier for the constraint norm(x) <= delta.
#         On exit par contains the final estimate of the multiplier.
#
#       f is a double precision variable.
#         On entry f need not be specified.
#         On exit f is set to f(x) at the output x.
#
#       x is a double precision array of dimension n.
#         On entry x need not be specified.
#         On exit x is set to the final estimate of the solution.
#
#       info is an integer variable.
#         On entry info need not be specified.
#         On exit info is set as follows:
#
#            info = 1  The function value f(x) has the relative
#                      accuracy specified by rtol.
#
#            info = 2  The function value f(x) has the absolute
#                      accuracy specified by atol.
#
#            info = 3  Rounding errors prevent further progress.
#                      On exit x is the best available approximation.
#
#            info = 4  Failure to converge after itmax iterations.
#                      On exit x is the best available approximation.
#
#       z is a double precision work array of dimension n.
#          code modified to return number of iterations in z(1)  (SBP 12/1/2006)
#
#       wa1 is a double precision work array of dimension n.
#
#       wa2 is a double precision work array of dimension n.
#
#     Subprograms called
#
#       MINPACK-2  ......  destsv
#
#       LAPACK  .........  dpotrf
#
#       Level 1 BLAS  ...  daxpy, dcopy, ddot, dnrm2, dscal
#
#       Level 2 BLAS  ...  dtrmv, dtrsv
#
#     MINPACK-2 Project. October 1993.
#     Argonne National Laboratory and University of Minnesota.
#     Brett M. Averick, Richard Carter, and Jorge J. More'
#
#     ***********
const p001 = 1.0
const p5 = 0.5

iter = 0

indef = 0
rznorm = 0

# double precision dasum, ddot, dnrm2
# external destsv, daxpy, dcopy, ddot, dnrm2, dscal, dtrmv, dtrsv

#     Initialization.

parf = zero(Float64)
xnorm = zero(Float64)
rxnorm = zero(Float64)
rednc = false
for j = 1:n
   x(j) = zero(Float64)
   z(j) = zero(Float64)
end

#     Copy the diagonal and save A in its lower triangle.

#all dcopy(n,a,lda+1,wa1,1)
for j = 1:(n - 1)
   call dcopy(n-j,a(j,j+1),lda,a(j+1,j),1)
end

#     Calculate the l1-norm of A, the Gershgorin row sums,
#     and the l2-norm of b.

anorm = zero(Float64)
for j = 1:n
   wa2(j) = dasum(n,a(1,j),1)
   anorm = max(anorm,wa2(j))
end
for j = 1:n
   wa2(j) = wa2(j) - abs(wa1(j))
end
bnorm = dnrm2(n,b,1)

#     Calculate a lower bound, pars, for the domain of the problem.
#     Also calculate an upper bound, paru, and a lower bound, parl,
#     for the Lagrange multiplier.

pars = -anorm
parl = -anorm
paru = -anorm
for j = 1:n
   pars = max(pars,-wa1(j))
   parl = max(parl,wa1(j)+wa2(j))
   paru = max(paru,-wa1(j)+wa2(j))
end
parl = max(zero(Float64),bnorm/delta-parl,pars)
paru = max(zero(Float64),bnorm/delta+paru)

#     If the input par lies outside of the interval (parl,paru),
#     set par to the closer endpoint.

par = max(par,parl)
par = min(par,paru)

#     Special case: parl = paru.

paru = max(paru,(one(Float64)+rtol)*parl)

#     Beginning of an iteration.

info = 0
for iter = 1:itmax

#        Safeguard par.

   if (par <= pars && paru > zero(Float64)) par = max(p001, sqrt(parl/paru))*paru end

#        Copy the lower triangle of A into its upper triangle and
#        compute A + par*I.

   for j = 1:(n - 1)
      call dcopy(n-j,a(j+1,j),1,a(j,j+1),lda)
   end
   for j = 1:n
      a(j,j) = wa1(j) + par
   end

#        Attempt the  Cholesky factorization of A without referencing
#        the lower triangular part.

   call dpotrf('U',n,a,lda,indef)

#        Case 1: A + par*I is positive definite.

   if (indef == 0)

#           Compute an approximate solution x and save the
#           last value of par with A + par*I positive definite.

      parf = par
      call dcopy(n,b,1,wa2,1)
      call dtrsv('U','T','N',n,a,lda,wa2,1)
      rxnorm = dnrm2(n,wa2,1)
      call dtrsv('U','N','N',n,a,lda,wa2,1)
      call dcopy(n,wa2,1,x,1)
      call dscal(n,-one(Float64),x,1)
      xnorm = dnrm2(n,x,1)

#           Test for convergence.

      if (abs(xnorm-delta) <= rtol*delta || (par == zero(Float64) && xnorm <= (one(Float64)+rtol)*delta)) info = 1 end

#           Compute a direction of negative curvature and use this
#           information to improve pars.

      call destsv(n,a,lda,rznorm,z)
      pars = max(pars,par-rznorm**2)

#           Compute a negative curvature solution of the form
#           x + alpha*z where norm(x+alpha*z) = delta.

      rednc = false
      if (xnorm < delta)

#              Compute alpha

         prod = ddot(n,z,1,x,1)/delta
         temp = (delta-xnorm)*((delta+xnorm)/delta)
         alpha = temp/(abs(prod)+sqrt(prod**2+temp/delta))
         alpha = sign(alpha,prod)

#              Test to decide if the negative curvature step
#              produces a larger reduction than with z = 0.

         rznorm = abs(alpha)*rznorm
         if ((rznorm/delta)**2+par*(xnorm/delta)**2 <= par) rednc = true

#              Test for convergence.

         if (p5*(rznorm/delta)**2 <= rtol*(one(Float64)-p5*rtol)*(par+(rxnorm/delta)**2))
            info = 1
         else if (p5*(par+(rxnorm/delta)**2) <= (atol/delta)/delta && info == 0)
            info = 2
         end
      end

#           Compute the Newton correction parc to par.

      if (xnorm == zero(Float64))
         parc = -par
      else
         call dcopy(n,x,1,wa2,1)
         temp = one(Float64)/xnorm
         call dscal(n,temp,wa2,1)
         call dtrsv('U','T','N',n,a,lda,wa2,1)
         temp = dnrm2(n,wa2,1)
         parc = (((xnorm-delta)/delta)/temp)/temp
      end

#           Update parl or paru.

      if (xnorm > delta) parl = max(parl,par) end
      if (xnorm < delta) paru = min(paru,par) end
   else

#           Case 2: A + par*I is not positive definite.

#           Use the rank information from the Cholesky
#           decomposition to update par.

      if (indef > 1)

#              Restore column indef to A + par*I.

         call dcopy(indef-1,a(indef,1),lda,a(1,indef),1)
         a(indef,indef) = wa1(indef) + par

#              Compute parc.

         call dcopy(indef-1,a(1,indef),1,wa2,1)
         call dtrsv('U','T','N',indef-1,a,lda,wa2,1)
         call dcopy(indef-1,wa2,1,a(1,indef),1)
         temp = dnrm2(indef-1,a(1,indef),1)
         a(indef,indef) = a(indef,indef) - temp**2
         call dtrsv('U','N','N',indef-1,a,lda,wa2,1)
      end
      wa2(indef) = -one(Float64)
      temp = dnrm2(indef,wa2,1)
      parc = -(a(indef,indef)/temp)/temp
      pars = max(pars,par+parc)

#           If necessary, increase paru slightly.
#           This is needed because in some exceptional situations
#           paru is the optimal value of par.

      paru = max(paru,(one(Float64)+rtol)*pars)
   end

#        Use pars to update parl.

   parl = max(parl,pars)

#        Test for termination.

   if (info == 0)
      if (iter == itmax) info = 4 end
      if (paru <= (one(Float64)+p5*rtol)*pars) info = 3 end
      if (paru == zero(Float64)) info = 2 end
   end

#        If exiting, store the best approximation and restore
#        the upper triangle of A.

   if (info != 0)

#           Compute the best current estimates for x and f.

      par = parf
      f = -p5*(rxnorm**2+par*xnorm**2)
      if (rednc)
         f = -p5*((rxnorm**2+par*delta**2)-rznorm**2)
         call daxpy(n,alpha,z,1,x,1)
      end

#           Restore the upper triangle of A.

      for j = 1:(n - 1)
         call dcopy(n-j,a(j+1,j),1,a(j,j+1),lda)
      end
      call dcopy(n,wa1,1,a,lda+1)
      z(1) = iter ! SBP: modification to return number of iterations
      return
   end

#        Compute an improved estimate for par.

   par = max(parl,par+parc)

#        End of an iteration.

end

end
