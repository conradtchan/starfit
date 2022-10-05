c$$$      program test
c$$$
c$$$      end


c#######################################################################
c     This algorithm for the log of the standard normal distribution is
c     written by
c
c     (C) Alexander Heger
c     Monash University
c     20140320
c
c     This function evaluates the log of the normal distribution function:
c
c                                 / x
c                        1       |       -t*t/2
c      ln(P(x)) =ln( ----------- |      e       dt )
c                    sqrt(2 pi)  |
c                                /-oo
c

      function logcdf(x)
      implicit None
      real*8 :: logcdf
      real*8, intent(in) :: x
      real*8 :: t
      real*8, external :: log1p
      real*8, parameter :: sqrt2i = sqrt(0.5d0)

      t = sqrt2i * x
      if (t > 1.d0) then
         logcdf = log1p(-0.5d0*erfc(t))
      else
         logcdf = log(erfc_scaled(-t)*0.5d0) - t**2
      endif
      end

c=======================================================================

c
c                 -0.5*x*x       __oo
c                e               \       n            -2*n
c     \Phi(x) = -------------- ( /_  (-1)  (2*i-1)!! x     )
c                x*sqrt(2*pi)    i=0
c
c
c     the first summand (i=0) is just 1
c     fastest method for x < -20 and -Ofast

      function logcdf2(x)
      implicit None
      real*8 :: logcdf2
      real*8, intent(in) :: x
      real*8, external :: logcdf
      real*8, parameter :: c = -0.5d0 * log(8.d0 * atan(1.d0))
      real*8, parameter :: sqrt2i = sqrt(0.5d0)

      real*8 :: r, i, f, x2, x2i, last

      if (x > -20.d0) then
         logcdf2 = logcdf(x)
         return
      endif

      x2 = x**2
      if (x < -1.d9) then
         logcdf2 = -0.5d0 * x2
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
      logcdf2 = -0.5d0 * x2 + log(-r/x) + c
      end
c=======================================================================
c     a slightly optimized version
c      function logcdf3(x)
c      implicit None
c      real*8 :: logcdf3
c      real*8, intent(in) :: x
c      real*8, external :: logcdf
c      real*8, parameter :: c = -0.5d0 * log(8.d0 * atan(1.d0))
c      real*8, parameter :: sqrt2i = sqrt(0.5d0)
c      real*8, parameter :: eps = epsilon(x)
c      real*8, parameter :: lim1 = 0.5d0*log(eps) !-18.03 for real*8
c      real*8, parameter :: lim2 = lim1/sqrt(eps) !-1.d9 for real*8
c
c      real*8 :: r, f, x2, x2i
c      integer :: i
c
c      integer, parameter :: n = -lim1*0.75
c      real*8, parameter, dimension(n) :: a = (/(1.d0 - 2.d0 * i,i=1,n)/)
c
c      if (x > lim1) then
c         logcdf3 = logcdf(x)
c         return
c      endif
c
c      x2 = x**2
c      if (x < lim2) then
c         logcdf3 = -0.5d0 * x2
c         return
c      endif
c
c      r = 1.d0
c      f = 1.d0
c      x2i = 1.d0 / x2
c      do i = 1, n
c         f = f * x2i * a(i)
c         if (abs(f) < eps) exit
c         r = r + f
c      enddo
c      logcdf3 = -0.5d0 * x2 + log(-r/x) + c
c      end



c#######################################################################
c     This algorithm for the standard normal distribution is written by
c
c     (C) Alexander Heger
c     Monash University
c     20140320
c
c     This function evaluates the log of the normal distribution function:
c
c                                 / x
c                        1       |       -t*t/2
c      ln(P(x)) =    ----------- |      e       dt
c                    sqrt(2 pi)  |
c                                /-oo
c

      function cdf(x)
      implicit None
      real*8 :: cdf
      real*8, intent(in) :: x
      real*8, parameter :: sqrt2i = -sqrt(0.5d0)

      cdf = 0.5d0*erfc(x * sqrt2i)
      end

c#######################################################################
c     This algorithm for the first derivative log of the standard normal
c     distribution is written by
c
c     (C) Alexander Heger
c     Monash University
c     20140329
c
c     This function evaluates the log of the normal distribution function:
c
c                                  / x
c      d                 -x*x/2   |       -t*t/2
c      -- ln(P(x)) =    e       / |      e       dt
c      dx                         |
c                                 /-oo
c

      function logcdfp(x)
      implicit None
      real*8 :: logcdfp
      real*8, intent(in) :: x
      real*8, parameter :: sqrt2i = -sqrt(0.5d0)
c     sqrt(2/pi)
      real*8, parameter :: fac = 1.d0 / sqrt(2.d0*atan(1.d0))

      logcdfp = fac / erfc_scaled(x * sqrt2i)

      end

c#######################################################################

c  This algorithm for the log of the standard normal distribution is
c    written by
c	Jean Marie Linhart
c	StataCorp LP
c
c    Last modified: October 26, 2007
c    converted to FORTRAN on 20140323 by Alexander Heger

      function lnenorm(xi)
      implicit None

      real*8 :: lnenorm
      real*8, intent(IN) :: xi

      real*8 :: x, ex, ans
      logical :: upper
      real*8, external :: lnnorm, log1p

      include 'norm.inc'

      ans = 0.d0

c     Next two lines represent limits on double precision arithmetic for
c     this problem.  One may want to test on one's own machine and
c     change -- the chips hold data on a register with more precision
c     than double, and this can effect results obtained.
c
c   LNNORM_MAX_X
c   LNNORM_MIN_X
c	defined in norminc.h

      if (xi > LNNORM_MAX_X) then
         lnenorm = 0.d0
         return
      endif
c     Call through to lnnorm to get correct behavior in the lower tail */
      if (xi < LNANORM_MIN_X) then
         lnenorm = lnnorm(xi)
         return
      endif

      if (x < 0.d0) then
         upper = .FALSE.
         x = -xi
      else
         upper = .TRUE.
         x = +xi
      endif

      if (x <= 0.5d0) then
         ex = 0.5d0 * erf(x/sqrt2)
         if (upper) then
            ans = log(0.5d0 + ex)
         else
            ans = log(0.5d0 - ex)
         endif
      else
         if (upper) then
            ex = 0.5d0 * erfc(x/sqrt2)
            if (ex < 0.1d0) then
               ans = LOG1P(-ex)
            else
               ans = log(1.d0 - ex)
            endif
         else
            ex = 0.5d0 * erfc(x/sqrt2)
            ans = log(ex)
         endif
      endif
      lnenorm = ans
      return
      end




c#######################################################################

c     This little function for log(1+x) is based on the Taylor series.
c     It is meant to be called only with small x.
c     Written by:
c	Jean Marie Linhart
c	StataCorp LP
c
c     Last modified October 26, 2007 (C version)
c
c     Adopted to FORTRAN 20140318 by Alexander Heger, Monash University
c
c     log1px takes a double and returns a double.
c     It is a Taylor series expansion of log(1+x).
c     x is presumed to be < 1.  As I have called it, x < .1,
c     and so I know the algorithm will terminate quickly.
c     The closer x is to 1, the slower this will be.

      pure function log1p(x)
      real*8 log1p
      real*8, intent(in) :: x

      integer :: n, sn
      real*8 xn, ans, oans, term, eps;

      if (abs(x) >= 1.d0) then
           log1p = 0.d0 !  NaN
           return
        endif
        term = x
        ans = x
        oans = x
        oans = ans + 1.d0
        n = 1
        sn = 1
        xn = x
c     Comparing ans!=oans is done here to insure that this calculation
c     continues until the accuracy of the machine is reached.  At some point,
c     the value is not being updated in successive iterations, that is time to
c     quit.
        do while (ans /= oans )
           oans = ans
           sn = -sn
           xn = xn * x
           n = n + 1
           term= (xn*sn)/n
           ans = ans + term
        enddo
        log1p = ans
        return
        end

c#######################################################################

c     This implementation for the log of the standard normal distribution
c     is by Jean Marie Linhart
c     StataCorp LP
c
c     It is based on the algorithm for the standard normal distribution
c     by W. J. Cody published in
c     Algorithm 715, Collected Algorithms from ACM.
c     This work published in Transactions on Mathematical Software,
c     Vol. 19, No. 1, March, 1993, pp. 22-32.
c
c     Includes changes given in remark by Price, TOMS 22 (2)
c
c------------------------------------------------------------------
c
c     This function evaluates the log of the normal distribution function:
c
c                                 / x
c                        1       |       -t*t/2
c      ln(P(x)) =ln( ----------- |      e       dt )
c                    sqrt(2 pi)  |
c                                /-oo
c
c     The main computation evaluates near-minimax approximations
c     derived from those in "Rational Chebyshev approximations for
c     the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
c     This transportable program uses rational functions that
c     theoretically approximate the normal distribution function to
c     at least 18 significant decimal digits.  The accuracy achieved
c     depends on the arithmetic system, the compiler, the intrinsic
c     functions, and proper selection of the machine-dependent
c     constants.
c
c *******************************************************************
c *******************************************************************
c
c Explanation of machine-dependent constants.  Let
c
c   xmin  = the smallest positive floating-point number.
c
c   Then the following machine-dependent constants must be declared
c   in DATA statements.  IEEE values are provided as a default.
c
c   eps   = argument below which anorm(x) may be represented by
c           0.5  and above which  x*x  will not underflow.
c           A conservative value is the largest machine number x
c           such that   1.0 + x = 1.0   to machine precision.
c   Approximate values for some important machines are:
c
c                          xmin        eps
c
c  CDC 7600      (S.P.)  3.13E-294   7.11E-15
c  CRAY-1        (S.P.)  4.58E-246   7.11E-157
c  IEEE (IBM/XT,
c    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8
c  IEEE (IBM/XT,
c    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16
c  IBM 195       (D.P.)  5.40D-79    1.39D-17
c  VAX D-Format  (D.P.)  2.94D-39    1.39D-17
c  VAX G-Format  (D.P.)  5.56D-309   1.11D-16
c
c *******************************************************************
c *******************************************************************
c
c Error returns
c
c  The program returns  anorm = 0     for  arg > LNNORM_MAX_X
c
c  Original normal distribution Algorithm Author: W. J. Cody
c          Mathematics and Computer Science Division
c          Argonne National Laboratory
c          Argonne, IL 60439
c  Author for the log of the normal distribution:
c	Jean Marie Linhart
c	StataCorp LP
c
c  Latest modification: January 3, 2008
c
c -----------------------------------------------------------------
c#include "norminc.h"

      function lnanorm(arg)
      implicit None

      real*8 :: lnanorm
      real*8, intent(in) :: arg

      integer :: i
      real*8 :: del, result, x, xden, xnum, y, xsq
      real*8, parameter :: eps = 1.11d-16
      real*8, parameter :: half = 0.5d0
      real*8, parameter :: one = 1.d0
      real*8, parameter :: sixten = 16.d0
      real*8, parameter :: zero = 0.d0

      real*8 :: log1p
      real*8 :: lnnorm

      include 'norm.inc'

c ------------------------------------------------------------------
c     Mathematical constants
c
c     sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), and
c     thrsh is the argument for which anorm = 0.75.
c
c ------------------------------------------------------------------
c     Next two lines represent limits on double precision arithmetic
c     for this problem.  One may want to test on one's own machine
c     and change -- the chips hold data on a register with more
c     precision than double, and this can effect results obtained.
c
c     LNNORM_MAX_X
c     LNNORM_MIN_X
c 	defined in norm.inc

c------------------------------------------------------------------
c Coefficients for approximation in first interval
c------------------------------------------------------------------
      real*8, parameter, dimension(0:4) :: a = (/
     &     2.2352520354606839287d+0,
     &     1.6102823106855587881d+2,
     &     1.0676894854603709582d+3,
     &     1.8154981253343561249d+4,
     &     6.5682337918207449113d-2 /)

      real*8, parameter, dimension(0:3) :: b = (/
     &     4.7202581904688241870d+1,
     &     9.7609855173777669322d+2,
     &     1.0260932208618978205d+4,
     &     4.5507789335026729956d+4  /)

c------------------------------------------------------------------
c Coefficients for approximation in second interval
c------------------------------------------------------------------
      real*8, parameter, dimension(0:8) :: c = (/
     &     3.9894151208813466764d-1,
     &     8.8831497943883759412d+0,
     &     9.3506656132177855979d+1,
     &     5.9727027639480026226d+2,
     &     2.4945375852903726711d+3,
     &     6.8481904505362823326d+3,
     &     1.1602651437647350124d+4,
     &     9.8427148383839780218d+3,
     &     1.0765576773720192317d-8 /)

      real*8, parameter, dimension(0:7) :: d = (/
     &     2.2266688044328115691d+1,
     &     2.3538790178262499861d+2,
     &     1.5193775994075548050d+3,
     &     6.4855582982667607550d+3,
     &     1.8615571640885098091d+4,
     &     3.4900952721145977266d+4,
     &     3.8912003286093271411d+4,
     &     1.9685429676859990727d+4 /)

c------------------------------------------------------------------
c Coefficients for approximation in third interval
c------------------------------------------------------------------
      real*8, parameter, dimension(0:5) :: p = (/
     &     2.158985340579569900d-1,
     &     1.274011611602473639d-1,
     &     2.223527787064980700d-2,
     &     1.421619193227893466d-3,
     &     2.911287495116879200d-5,
     &     2.307344176494017303d-2 /)

      real*8, parameter, dimension(0:4) :: q = (/
     &     1.28426009614491121d+0,
     &     4.68238212480865118d-1,
     &     6.59881378689285515d-2,
     &     3.78239633202758244d-3,
     &     7.29751555083966205d-5  /)

c------------------------------------------------------------------
      if (arg > LNNORM_MAX_X) then
         lnanorm = 0.d0
         return
      endif
c     Call through to lnnorm to get correct behavior in the lower tail
      if (arg < LNANORM_MIN_X) then
         lnanorm = lnnorm(arg)
         return
      endif
c*------------------------------------------------------------------
	x = arg
	y = abs(x)
        if (y <= thrsh) then
c------------------------------------------------------------------
c     Evaluate  lnanorm  for  |x| <= 0.67448975
c------------------------------------------------------------------
           xsq = zero
           if (y > eps) xsq = x**2
           xnum = a(4)*xsq
           xden = xsq
           do i=0,2
              xnum = (xnum + a(i)) * xsq
              xden = (xden + b(i)) * xsq
           enddo
           result = x * (xnum + a(3)) / (xden + b(3))
           result = log(half + result)
c------------------------------------------------------------------
c     Evaluate  lnanorm  for 0.67448975 <= |x| <= sqrt(32)
c------------------------------------------------------------------
	else if (y <= root32) then
           xnum = c(8)*y
           xden = y
           do i=0,6
              xnum = (xnum + c(i)) * y
              xden = (xden + d(i)) * y
           enddo
           result = (xnum + c(7)) / (xden + d(7))
           xsq = aint(y*sixten) / sixten
           del = (y-xsq)*(y+xsq)
           result = exp(-xsq*xsq*half)*exp(-del*half)*result
           if (x > zero) then
              if (result < 0.1d0) then ! close to 0
                 result = LOG1P(-result)
              else ! far from 0
                 result = log(one-result)
              endif
           else
              result = log(result)
           endif
c------------------------------------------------------------------
c     Evaluate  lnanorm  for |x| > sqrt(32)
c------------------------------------------------------------------
	else
           result = zero
           xsq = one / x ** 2
           xnum = p(5)*xsq
           xden = xsq
           do i=0,3
              xnum = (xnum + p(i)) * xsq
              xden = (xden + q(i)) * xsq
           enddo
           result = xsq *(xnum + p(4)) / (xden + q(4))
           result = (sqrpi -  result) / y
           xsq = aint(x * sixten) / sixten
           del = (x - xsq) * (x + xsq)
           result = exp(-xsq*xsq*half)*exp(-del*half)*result
           if (x > zero) then
              result = LOG1P(-result)
           else
              result = log(result)
           endif
        endif
        lnanorm = result
        return
      end

c#######################################################################

c     C code for the lnnorm(x) function to return the logarithm of the normal
c     distribution.
c
c     It is based on the attached code for the normal distribution, originally
c     written at StataCorp LP, and based on Algorithm 304 and the remarks by
c     Adams and Holmgren.
c
c     This code is written by Jean Marie Linhart
c                             StataCorp LP
c                             jlinhart@stata.com
c
c     Last modified January 4, 2008
c
c     Ported to FORTRAN 20140323 by Alexander Heger
c
c     The following is a modification of Algorithm 304 to give the log of the
c     normal distribution function at x.
c
c     It takes a double precision argument and returns double precision
c

      function lnnorm(zi)

      implicit None

      real*8 :: lnnorm
      real*8, intent(IN) :: zi

      logical :: lower
      real*8 :: z, z2, y, s, p1, q1, p2, q2, t, a1, a2
      real*8 :: n, m

      real*8, external :: log1p

      include 'norm.inc'

      z = zi

      if (z == 0.d0) then
         lnnorm = log(0.5d0)
         return
      endif

c     Next two lines represent limits on double precision arithmetic for
c     this problem.  One may want to test on one's own machine and
c     change -- the chips hold data on a register with more precision
c     than double, and this can effect results obtained.
c
c
c     LNNORM_MAX_X
c     LNNORM_MIN_X
c         defined in norm.inc
c
      if (z > LNNORM_MAX_X) then
         lnnorm = 0.d0
         return
      endif
      if (z <= LNNORM_MIN_X) then
         lnnorm = -0.5d0 * z**2
         return
      endif


c        In the original algorithm, the logical variable upper served a
c        dual purpose.  On input, it indicated whether the user wanted
c        F(z) or 1-F(z).  Inside the routine, it combined this with
c        information on the positivity of z.  In this version, the
c        algorithm always returns F(z)

      if (z < 0.d0) then
         z = -z
         lower = .True.
      else
         lower = .False.
      endif

      z2 = z**2

c     y is the standard normal density function evaluated at z.

      y = exp(-0.5d0 * z2) / sqr2pi
      n = y / z

c     The original Algorithm 304 had z < 2.32 or z < 3.5 depending on
c     whether you were doing upper or lower.  This was changed for
c     (slightly) greater accuracy
      if ( z <= 2.d0 ) then
c     This is the series representation for F(z)
c     where z lies near center of the distribution
         z = z * y
         s = z
         t = 0.d0
c     Comparing s!=t is done here to insure that this calculation
c     continues until the accuracy of the machine is reached.  At some
c     point, the value is not being updated in successive iterations,
c     that is time to quit.
         n = 3.d0
         do while (s /= t)
            t = s
            z = z * (z2 / n)
            s = s + z
            n = n + 2.d0
         enddo
         if (lower) then
            lnnorm = log(0.5d0 - s)
         else
            lnnorm = log(0.5d0 + s)
         endif
         return
      endif

c     This is the continued fraction representation for Q(z) which
c     is used to evaluate F(z) when z lies in one of the tails
c     of this distribution.  This section contains the modifications
c     suggested by Homgren and Adams.


        a1 = 2.d0
        a2 = 0.d0
        n = z2 + 3.d0
        p1 = 1.d0
        q1 = z
        p2 = n - 1.d0
        q2 = n * z
        m = p1 / q1
        t = p2 / q2
        s = m  ! dummy assignment to not stop while

c Comparing m!=t and s!=t is done here to insure that this calculation
c continues until the accuracy of the machine is reached.  At some
c point, the value is not being updated in successive iterations, that
c is time to quit.

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
           if (q2 > 1.d30) then
              p1 = p1 / 1.d30
              p2 = p2 / 1.d30
              q1 = q1 / 1.d30
              q2 = q2 / 1.d30
           endif
           t = p2 / q2
        enddo
        if (lower) then
           lnnorm = log(t) - 0.5d0 * z2 - log(sqr2pi)
        else
           lnnorm =  log1p(-y*t)
        endif
        return
        end

c#######################################################################

c     This algorithm for the standard normal distribution is written by
c     Jean Marie Linhart
c     StataCorp LP
c
c     Last modified: October 16, 2007
c     converted to Fortran 20140324 by Alexander Heger

      function enorm(xi)

      implicit none

      real*8 :: enorm
      real*8, intent (IN) :: xi

      real*8 :: ex, ans, x
      logical :: upper

      real*8, parameter :: sqrt2i = sqrt(0.5d0)
      real*8, parameter :: sqrt8i = sqrt(0.125d0)

      include 'norm.inc'

      enorm = 0.d0

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
            enorm = 0.5d0 + ex
         else
            enorm = 0.5d0 - ex
         endif
      else
         ex = 0.5d0 * erfc(x)
         if (upper) then
            enorm = 1.d0 - ex
         else
            enorm = ex
         endif
      endif
      return
      end
