function estsv(n::Integer,r::AbstractArray{Float64,1},ldr::Integer,svmin::Float64,z::AbstractArray{Float64,1})
      #TODO assert size(r) == (ldr,n), size(z) == (n)

# Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
# see CopyrightMINPACK.txt
#
#     **********
#
#     Subroutine destsv
#
#     Given an n by n upper triangular matrix R, this subroutine
#     estimates the smallest singular value and the associated
#     singular vector of R.
#
#     In the algorithm a vector e is selected so that the solution
#     y to the system R'*y = e is large. The choice of sign for the
#     components of e cause maximal local growth in the components
#     of y as the forward substitution proceeds. The vector z is
#     the solution of the system R*z = y, and the estimate svmin
#     is norm(y)/norm(z) in the Euclidean norm.
#
#     The subroutine statement is
#
#       subroutine estsv(n,r,ldr,svmin,z)
#
#     where
#
#       n is an integer variable.
#         On entry n is the order of R.
#         On exit n is unchanged.
#
#       r is a double precision array of dimension (ldr,n)
#         On entry the full upper triangle must contain the full
#            upper triangle of the matrix R.
#         On exit r is unchanged.
#
#       ldr is an integer variable.
#         On entry ldr is the leading dimension of r.
#         On exit ldr is unchanged.
#
#       svmin is a double precision variable.
#         On entry svmin need not be specified.
#         On exit svmin contains an estimate for the smallest
#            singular value of R.
#
#       z is a double precision array of dimension n.
#         On entry z need not be specified.
#         On exit z contains a singular vector associated with the
#            estimate svmin such that norm(R*z) = svmin and
#            norm(z) = 1 in the Euclidean norm.
#
#     Subprograms called
#
#       Level 1 BLAS ... dasum, daxpy, dnrm2, dscal
#
#     MINPACK-2 Project. October 1993.
#     Argonne National Laboratory
#     Brett M. Averick and Jorge J. More'.
#
#     **********
      double precision one, p01, zero
      parameter (zero=0.d0,p01=1.0d-2,one=1.0d0)

      integer i, j
      double precision e, s, sm, temp, w, wm, ynorm, znorm

      double precision dasum, dnrm2
      external dasum, daxpy, dnrm2, dscal

      do i = 1, n
         z(i) = zero
      end do

#     This choice of e makes the algorithm scale invariant.

      e = abs(r(1,1))
      if (e .eq. zero) then
         svmin = zero
         z(1) = one
         return
      end if

#     Solve R'*y = e.

      do i = 1, n

#        Scale y. The factor of 0.01 reduces the number of scalings.

         e = sign(e,-z(i))
         if (abs(e-z(i)) .gt. abs(r(i,i))) then
            temp = min(p01,abs(r(i,i))/abs(e-z(i)))
            call dscal(n,temp,z,1)
            e = temp*e
         end if

#        Determine the two possible choices of y(i).

         if (r(i,i) .eq. zero) then
            w = one
            wm = one
         else
            w = (e-z(i))/r(i,i)
            wm = -(e+z(i))/r(i,i)
         end if

#        Choose y(i) based on the predicted value of y(j) for j > i.

         s = abs(e-z(i))
         sm = abs(e+z(i))
         do j = i + 1, n
            sm = sm + abs(z(j)+wm*r(i,j))
         end do
         if (i .lt. n) then
            call daxpy(n-i,w,r(i,i+1),ldr,z(i+1),1)
            s = s + dasum(n-i,z(i+1),1)
         end if
         if (s .lt. sm) then
            temp = wm - w
            w = wm
            if (i .lt. n) call daxpy(n-i,temp,r(i,i+1),ldr,z(i+1),1)
         end if
         z(i) = w

      end do

      ynorm = dnrm2(n,z,1)

#     Solve R*z = y.

      do j = n, 1, -1

#        Scale z.

         if (abs(z(j)) .gt. abs(r(j,j))) then
            temp = min(p01,abs(r(j,j))/abs(z(j)))
            call dscal(n,temp,z,1)
            ynorm = temp*ynorm
         end if
         if (r(j,j) .eq. zero) then
            z(j) = one
         else
            z(j) = z(j)/r(j,j)
         end if
         temp = -z(j)
         call daxpy(j-1,temp,r(1,j),1,z,1)

      end do

#     Compute svmin and normalize z.

      znorm = one/dnrm2(n,z,1)
      svmin = ynorm*znorm
      call dscal(n,znorm,z,1)

      end
