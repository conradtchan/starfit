module norm

  use typedef, only: &
       real64, int64

  implicit none

  private

  ! these are constants relating to double precision arithmetic */

  real(kind=real64), parameter :: LNNORM_MAX_X = 38.d0
  real(kind=real64), parameter :: LNNORM_MIN_X = -1.d9
  real(kind=real64), parameter :: LNANORM_MIN_X = -37.519d0

  ! Some more constants.  If upping precision from double, these next
  ! three will need to be made more precise.

  real(kind=real64), parameter :: pi = 4.d0*atan(1.d0)
  real(kind=real64), parameter :: sqr2pi = sqrt(2.d0*pi)
  real(kind=real64), parameter :: sqrpi = 1.d0 / sqr2pi
  real(kind=real64), parameter :: sqrt2 = sqrt(2.d0)

  ! these next two are used for comparison, not calculation, they do not
  ! need to be as precise as the above

  real(kind=real64), parameter :: thrsh = 0.67448975d0
  real(kind=real64), parameter :: root32 = sqrt(32.d0)

  ! some more simple constants

  real(kind=real64), parameter :: sqrt2i = sqrt(0.5d0)
  real(kind=real64), parameter :: sqrt2im = -sqrt(0.5d0)

  public :: &
       logcdf, logcdfp, &
       logcdf2, logcdf2a, cdf, lnenorm, log1p, lnanorm, lnnorm, enorm

contains

!#######################################################################
!     This algorithm for the log of the standard normal distribution is
!     written by
!
!     (C) Alexander Heger
!     Monash University
!     20140320
!
!     This function evaluates the log of the normal distribution function:
!
!                                 / x
!                        1       |       -t*t/2
!      ln(P(x)) =ln( ----------- |      e       dt )
!                    sqrt(2 pi)  |
!                                /-oo
!

  elemental pure function logcdf(x) result(y)

    implicit None

    real(kind=real64), intent(in) :: &
         x
    real(kind=real64) :: &
         y

    real(kind=real64) :: &
         t

    t = sqrt2i * x
    if (t > 1.d0) then
       y = log1p(-0.5d0*erfc(t))
    else
       y = log(erfc_scaled(-t)*0.5d0) - t**2
    endif

  end function logcdf

!=======================================================================

!
!                 -0.5*x*x       __oo
!                e               \       n            -2*n
!     \Phi(x) = -------------- ( /_  (-1)  (2*i-1)!! x     )
!                x*sqrt(2*pi)    i=0
!
!
!     the first summand (i=0) is just 1
!     fastest method for x < -20 and -Ofast

  elemental pure function logcdf2(x) result(y)

    implicit None

    real(kind=real64), parameter :: &
         c = -0.5d0 * log(8.d0 * atan(1.d0))

    real(kind=real64), intent(in) :: &
         x

    real(kind=real64) :: &
         y

    real(kind=real64) :: &
         r, i, f, x2, x2i, last

    if (x > -20.d0) then
       y = logcdf(x)
       return
    endif

    x2 = x**2
    if (x < -1.d9) then
       y = -0.5d0 * x2
       return
    endif

    r = 1.d0
    i = 0.d0
    f = 1.d0
    x2i = 1.d0 / x2
    last = 0.d0
    DO WHILE (last /= r)
       i = i + 1.d0
       f = f * x2i * (1.d0 - 2.d0 * i)
       last = r
       r = last + f
    enddo
    y = -0.5d0 * x2 + log(-r/x) + c

  end function logcdf2


!=======================================================================
!     a slightly optimized version
  elemental pure function logcdf2a(x) result(y)

    implicit None

    real(kind=real64), intent(in) :: &
         x

    real(kind=real64) :: &
         y

    real(kind=real64), parameter :: &
         c = -0.5d0 * log(8.d0 * atan(1.d0))
    real(kind=real64), parameter :: &
         eps = epsilon(x)
    real(kind=real64), parameter :: &
         lim1 = 0.5d0*log(eps) !-18.03 for real(kind=real64)
    real(kind=real64), parameter :: &
         lim2 = lim1/sqrt(eps) !-1.d9 for real(kind=real64)

    real(kind=real64) :: &
         r, f, x2, x2i
    integer(kind=int64) :: &
         i

    integer(kind=int64), parameter :: &
         n = ceiling(-lim1 * 0.75d0)
    real(kind=real64), parameter, dimension(n) :: &
         a = (/(1.d0 - 2.d0 * i,i=1,n)/)

     if (x > lim1) then
        y = logcdf(x)
        return
     endif

     x2 = x**2
     if (x < lim2) then
        y = -0.5d0 * x2
        return
     endif

     r = 1.d0
     f = 1.d0
     x2i = 1.d0 / x2
     do i = 1, n
        f = f * x2i * a(i)
        if (abs(f) < eps) exit
        r = r + f
     enddo
     y = -0.5d0 * x2 + log(-r/x) + c

   end function logcdf2a


!#######################################################################
!     This algorithm for the standard normal distribution is written by
!
!     (C) Alexander Heger
!     Monash University
!     20140320
!
!     This function evaluates the log of the normal distribution function:
!
!                                 / x
!                        1       |       -t*t/2
!      ln(P(x)) =    ----------- |      e       dt
!                    sqrt(2 pi)  |
!                                /-oo

   elemental pure function cdf(x) result(y)

     implicit None

     real(kind=real64), intent(in) :: &
          x

     real(kind=real64) :: &
          y

     y = 0.5d0*erfc(x * sqrt2im)

   end function cdf

!#######################################################################
!     This algorithm for the first derivative log of the standard normal
!     distribution is written by
!
!     (C) Alexander Heger
!     Monash University
!     20140329
!
!     This function evaluates the log of the normal distribution function:
!
!                                  / x
!      d                 -x*x/2   |       -t*t/2
!      -- ln(P(x)) =    e       / |      e       dt
!      dx                         |
!                                 /-oo
!

  elemental pure function logcdfp(x) result(y)

    implicit None

    real(kind=real64), intent(in) :: &
         x

    real(kind=real64) :: &
         y

    ! sqrt(2/pi)
    real(kind=real64), parameter :: &
         fac = 1.d0 / sqrt(2.d0*atan(1.d0))

    y = fac / erfc_scaled(x * sqrt2im)

  end function logcdfp

!#######################################################################

!  This algorithm for the log of the standard normal distribution is
!    written by
!	Jean Marie Linhart
!	StataCorp LP
!
!    Last modified: October 26, 2007
!    converted to FORTRAN on 20140323 by Alexander Heger

  elemental pure function lnenorm(xi) result(y)

    implicit none

    real(kind=real64), intent(in) :: &
         xi

    real(kind=real64) :: &
         y

    real(kind=real64) :: &
         x, ex

    logical :: upper

    ! Next two lines represent limits on double precision arithmetic for
    ! this problem.  One may want to test on one's own machine and
    ! change -- the chips hold data on a register with more precision
    ! than double, and this can effect results obtained.

    if (xi > LNNORM_MAX_X) then
       y = 0.d0
       return
    endif

    ! Call through to lnnorm to get correct behavior in the lower tail

    if (xi < LNANORM_MIN_X) then
       y = lnnorm(xi)
       return
    endif

    if (xi < 0.d0) then
       upper = .FALSE.
       x = -xi
    else
       upper = .TRUE.
       x = +xi
    endif

    if (x <= 0.5d0) then
       ex = 0.5d0 * erf(x/sqrt2)
       if (upper) then
          y = log(0.5d0 + ex)
       else
          y = log(0.5d0 - ex)
       endif
    else
       if (upper) then
          ex = 0.5d0 * erfc(x/sqrt2)
          if (ex < 0.1d0) then
             y = LOG1P(-ex)
          else
             y = log(1.d0 - ex)
          endif
       else
          ex = 0.5d0 * erfc(x/sqrt2)
          y = log(ex)
       endif
    endif

  end function lnenorm


!#######################################################################

!     This little function for log(1+x) is based on the Taylor series.
!     It is meant to be called only with small x.
!     Written by:
!	Jean Marie Linhart
!	StataCorp LP
!
!     Last modified October 26, 2007 (C version)
!
!     Adopted to FORTRAN 20140318 by Alexander Heger, Monash University
!
!     log1px takes a double and returns a double.
!     It is a Taylor series expansion of log(1+x).
!     x is presumed to be < 1.  As I have called it, x < .1,
!     and so I know the algorithm will terminate quickly.
!     The closer x is to 1, the slower this will be.

  elemental pure function log1p(x) result(y)

    implicit none

    real(kind=real64), intent(in) :: &
         x

    real(kind=real64) :: &
         y

    integer(kind=int64) :: &
         n, sn

    real(kind=real64) :: &
         xn, oans, term

    if (abs(x) >= 1.d0) then
       y = 0.d0 !  NaN
       return
    endif

    term = x
    y = x
    oans = x
    oans = y + 1.d0
    n = 1
    sn = 1
    xn = x

    ! Comparing ans!=oans is done here to insure that this calculation
    ! continues until the accuracy of the machine is reached.  At some
    ! point, the value is not being updated in successive iterations,
    ! that is time to quit.

    do while (y /= oans )
       oans = y
       sn = -sn
       xn = xn * x
       n = n + 1
       term= (xn*sn)/n
       y = y + term
    enddo

  end function log1p

!#######################################################################

!     This implementation for the log of the standard normal distribution
!     is by Jean Marie Linhart
!     StataCorp LP
!
!     It is based on the algorithm for the standard normal distribution
!     by W. J. Cody published in
!     Algorithm 715, Collected Algorithms from ACM.
!     This work published in Transactions on Mathematical Software,
!     Vol. 19, No. 1, March, 1993, pp. 22-32.
!
!     Includes changes given in remark by Price, TOMS 22 (2)
!
!------------------------------------------------------------------
!
!     This function evaluates the log of the normal distribution function:
!
!                                 / x
!                        1       |       -t*t/2
!      ln(P(x)) =ln( ----------- |      e       dt )
!                    sqrt(2 pi)  |
!                                /-oo
!
!     The main computation evaluates near-minimax approximations
!     derived from those in "Rational Chebyshev approximations for
!     the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!     This transportable program uses rational functions that
!     theoretically approximate the normal distribution function to
!     at least 18 significant decimal digits.  The accuracy achieved
!     depends on the arithmetic system, the compiler, the intrinsic
!     functions, and proper selection of the machine-dependent
!     constants.
!
! *******************************************************************
! *******************************************************************
!
! Explanation of machine-dependent constants.  Let
!
!   xmin  = the smallest positive floating-point number.
!
!   Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   eps   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number x
!           such that   1.0 + x = 1.0   to machine precision.
!   Approximate values for some important machines are:
!
!                          xmin        eps
!
!  CDC 7600      (S.P.)  3.13E-294   7.11E-15
!  CRAY-1        (S.P.)  4.58E-246   7.11E-157
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16
!  IBM 195       (D.P.)  5.40D-79    1.39D-17
!  VAX D-Format  (D.P.)  2.94D-39    1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   1.11D-16
!
! *******************************************************************
! *******************************************************************
!
! Error returns
!
!  The program returns  anorm = 0     for  arg > LNNORM_MAX_X
!
!  Original normal distribution Algorithm Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!  Author for the log of the normal distribution:
!	Jean Marie Linhart
!	StataCorp LP
!
!  Latest modification: January 3, 2008
!
! -----------------------------------------------------------------

  elemental pure function lnanorm(arg) result(y)

    implicit None

    real(kind=real64), intent(in) :: &
         arg

    real(kind=real64) :: &
         y

    integer(kind=int64) :: &
         i
    real(kind=real64) :: &
         del, x, xden, xnum, xa, xsq

    real(kind=real64), parameter :: &
         eps = 1.11d-16, &
         half = 0.5d0, &
         one = 1.d0, &
         sixten = 16.d0, &
         zero = 0.d0

! ------------------------------------------------------------------
!     Mathematical constants
!
!     sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), and
!     thrsh is the argument for which anorm = 0.75.
!
! ------------------------------------------------------------------
!     Next two lines represent limits on double precision arithmetic
!     for this problem.  One may want to test on one's own machine
!     and change -- the chips hold data on a register with more
!     precision than double, and this can effect results obtained.
!
!------------------------------------------------------------------
! Coefficients for approximation in first interval
!------------------------------------------------------------------
    real*8, parameter, dimension(0:4) :: a = (/ &
         2.2352520354606839287d+0, &
         1.6102823106855587881d+2, &
         1.0676894854603709582d+3, &
         1.8154981253343561249d+4, &
         6.5682337918207449113d-2 /)

    real*8, parameter, dimension(0:3) :: b = (/ &
         4.7202581904688241870d+1, &
         9.7609855173777669322d+2, &
         1.0260932208618978205d+4, &
         4.5507789335026729956d+4  /)

!------------------------------------------------------------------
! Coefficients for approximation in second interval
!------------------------------------------------------------------
    real*8, parameter, dimension(0:8) :: c = (/ &
         3.9894151208813466764d-1, &
         8.8831497943883759412d+0, &
         9.3506656132177855979d+1, &
         5.9727027639480026226d+2, &
         2.4945375852903726711d+3, &
         6.8481904505362823326d+3, &
         1.1602651437647350124d+4, &
         9.8427148383839780218d+3, &
         1.0765576773720192317d-8 /)

    real*8, parameter, dimension(0:7) :: d = (/ &
         2.2266688044328115691d+1, &
         2.3538790178262499861d+2, &
         1.5193775994075548050d+3, &
         6.4855582982667607550d+3, &
         1.8615571640885098091d+4, &
         3.4900952721145977266d+4, &
         3.8912003286093271411d+4, &
         1.9685429676859990727d+4 /)

!------------------------------------------------------------------
! Coefficients for approximation in third interval
!------------------------------------------------------------------
    real*8, parameter, dimension(0:5) :: p = (/ &
         2.158985340579569900d-1, &
         1.274011611602473639d-1, &
         2.223527787064980700d-2, &
         1.421619193227893466d-3, &
         2.911287495116879200d-5, &
         2.307344176494017303d-2 /)

    real*8, parameter, dimension(0:4) :: q = (/ &
         1.28426009614491121d+0, &
         4.68238212480865118d-1, &
         6.59881378689285515d-2, &
         3.78239633202758244d-3, &
         7.29751555083966205d-5  /)

!------------------------------------------------------------------
    if (arg > LNNORM_MAX_X) then
       y = 0.d0
       return
    endif

    ! Call through to lnnorm to get correct behavior in the lower tail

    if (arg < LNANORM_MIN_X) then
       y = lnnorm(arg)
       return
    endif
!------------------------------------------------------------------
    x = arg
    xa = abs(x)
    if (xa <= thrsh) then
!------------------------------------------------------------------
!     Evaluate  lnanorm  for  |x| <= 0.67448975
!------------------------------------------------------------------
       xsq = zero
       if (xa > eps) xsq = x**2
       xnum = a(4)*xsq
       xden = xsq
       do i=0,2
          xnum = (xnum + a(i)) * xsq
          xden = (xden + b(i)) * xsq
       enddo
       y = x * (xnum + a(3)) / (xden + b(3))
       y = log(half + y)
!------------------------------------------------------------------
!     Evaluate  lnanorm  for 0.67448975 <= |x| <= sqrt(32)
!------------------------------------------------------------------
    else if (xa <= root32) then
       xnum = c(8) * xa
       xden = xa
       do i=0,6
          xnum = (xnum + c(i)) * xa
          xden = (xden + d(i)) * xa
       enddo
       y = (xnum + c(7)) / (xden + d(7))
       xsq = aint(xa * sixten) / sixten
       del = (xa - xsq)*(y + xsq)
       y = exp(-xsq*xsq*half)*exp(-del*half)*y
       if (x > zero) then
          if (y < 0.1d0) then ! close to 0
             y = LOG1P(-y)
          else ! far from 0
             y = log(one-y)
          endif
       else
          y = log(y)
       endif
!------------------------------------------------------------------
!     Evaluate  lnanorm  for |x| > sqrt(32)
!------------------------------------------------------------------
    else
       y = zero
       xsq = one / x ** 2
       xnum = p(5)*xsq
       xden = xsq
       do i=0,3
          xnum = (xnum + p(i)) * xsq
          xden = (xden + q(i)) * xsq
       enddo
       y = xsq *(xnum + p(4)) / (xden + q(4))
       y = (sqrpi -  y) / xa
       xsq = aint(x * sixten) / sixten
       del = (x - xsq) * (x + xsq)
       y = exp(-xsq*xsq*half)*exp(-del*half)*y
       if (x > zero) then
          y = LOG1P(-y)
       else
          y = log(y)
       endif
    endif

  end function lnanorm

!#######################################################################

!     code for the lnnorm(x) function to return the logarithm of the normal
!     distribution.
!
!     It is based on the attached code for the normal distribution, originally
!     written at StataCorp LP, and based on Algorithm 304 and the remarks by
!     Adams and Holmgren.
!
!     This code is written by Jean Marie Linhart
!                             StataCorp LP
!                             jlinhart@stata.com
!
!     Last modified January 4, 2008
!
!     Ported to FORTRAN 20140323 by Alexander Heger
!
!     The following is a modification of Algorithm 304 to give the log of the
!     normal distribution function at x.
!
!     It takes a double precision argument and returns double precision
!

  elemental pure function lnnorm(zi) result(y)

    implicit None

    real(kind=real64), parameter :: &
         q2lim = 1.d30, &
         q2limi = 1.d0 / q2lim

    real(kind=real64), intent(IN) :: &
         zi

    real(kind=real64) :: &
         y

    logical :: &
         lower
    real(kind=real64) :: &
         z, z2, f, s, p1, q1, p2, q2, t, a1, a2
    real(kind=real64) :: &
         n, m

    z = zi

    if (z == 0.d0) then
       y = log(0.5d0)
       return
    endif

    ! Next two lines represent limits on double precision arithmetic
    ! for this problem.  One may want to test on one's own machine and
    ! change -- the chips hold data on a register with more precision
    ! than double, and this can effect results obtained.

    if (z > LNNORM_MAX_X) then
       y = 0.d0
       return
    endif
    if (z <= LNNORM_MIN_X) then
       y = -0.5d0 * z**2
       return
    endif

    ! In the original algorithm, the logical variable upper served a
    ! dual purpose.  On input, it indicated whether the user wanted
    ! F(z) or 1-F(z).  Inside the routine, it combined this with
    ! information on the positivity of z.  In this version, the
    ! algorithm always returns F(z)

    if (z < 0.d0) then
       z = -z
       lower = .True.
    else
       lower = .False.
    endif

    z2 = z**2

    ! f is the standard normal density function evaluated at z.

    f = exp(-0.5d0 * z2) / sqr2pi
    n = f / z

    ! The original Algorithm 304 had z < 2.32 or z < 3.5 depending on
    ! whether you were doing upper or lower.  This was changed for
    ! (slightly) greater accuracy

    if ( z <= 2.d0 ) then

       ! This is the series representation for F(z) where z lies near
       ! center of the distribution

       z = z * f
       s = z
       t = 0.d0

       ! Comparing s!=t is done here to insure that this calculation
       ! continues until the accuracy of the machine is reached.  At
       ! some point, the value is not being updated in successive
       ! iterations, that is time to quit.

       n = 3.d0
       do while (s /= t)
          t = s
          z = z * (z2 / n)
          s = s + z
          n = n + 2.d0
       enddo
       if (lower) then
          y = log(0.5d0 - s)
       else
          y = log(0.5d0 + s)
       endif
       return
    endif

    ! This is the continued fraction representation for Q(z) which
    ! is used to evaluate F(z) when z lies in one of the tails
    ! of this distribution.  This section contains the modifications
    ! suggested by Homgren and Adams.

    a1 = 2.d0
    a2 = 0.d0
    n = z2 + 3.d0
    p1 = 1.d0
    q1 = z
    p2 = n - 1.d0
    q2 = n * z
    m = p1 / q1
    t = p2 / q2
    s = m ! dummy assignment to not stop while

    ! Comparing m!=t and s!=t is done here to insure that this calculation
    ! continues until the accuracy of the machine is reached.  At some
    ! point, the value is not being updated in successive iterations, that
    ! is time to quit.

    do while ((m /= t) .and. (s /= t))
       n = n + 4.d0
       a1 = a1 - 8.d0
       a2 = a2 + a1
       s = a2*p1 + n*p2
       p1 = p2
       p2 = s
       s = a2 * q1 + n * q2
       q1 = q2
       q2 = s
       s = m
       m = t
       if (q2 > q2lim) then
          p1 = p1 * q2limi
          p2 = p2 * q2limi
          q1 = q1 * q2limi
          q2 = q2 * q2limi
       endif
       t = p2 / q2
    enddo
    if (lower) then
       y = log(t) - 0.5d0 * z2 - log(sqr2pi)
    else
       y = log1p(-f*t)
    endif

  end function lnnorm

!#######################################################################

!     This algorithm for the standard normal distribution is written by
!     Jean Marie Linhart
!     StataCorp LP
!
!     Last modified: October 16, 2007
!     converted to Fortran 20140324 by Alexander Heger

  elemental pure function enorm(xi) result(y)

    implicit none

    real(kind=real64), parameter :: &
         sqrt8i = sqrt(0.125d0)

    real(kind=real64), intent (IN) :: &
         xi

    real(kind=real64) :: &
         y

    real(kind=real64) :: &
         ex, x

    logical :: &
         upper

    if (xi < 0.d0) then
       upper = .FALSE.
       x = -xi * sqrt2i
    else
       upper = .TRUE.
       x = +xi * sqrt2i
    endif

    if (x <= sqrt8i) then
       ex = 0.5d0 * erf(x)
       if (upper) then
          y = 0.5d0 + ex
       else
          y = 0.5d0 - ex
       endif
    else
       ex = 0.5d0 * erfc(x)
       if (upper) then
          y = 1.d0 - ex
       else
          y = ex
       endif
    endif
  end function enorm

end module norm
