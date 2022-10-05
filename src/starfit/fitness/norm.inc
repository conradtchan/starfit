c     these are constants relating to double precision arithmetic */

      real*8, parameter :: LNNORM_MAX_X = 38.d0
      real*8, parameter :: LNNORM_MIN_X = -1.d9
      real*8, parameter :: LNANORM_MIN_X = -37.519d0

c     Some more constants.  If upping precision from double, these next
c     three will need to be made more precise.

c$$$      real*8, parameter :: sqrpi = 3.989422804014326779399461d-1
c$$$      real*8, parameter :: sqr2pi = 2.506628274631000502415765d0
c$$$      real*8, parameter :: sqrt2 = 1.4142135623730950488d0

      real*8, parameter :: pi = 4.d0*atan(1.d0)
      real*8, parameter :: sqr2pi = sqrt(2.d0*pi)
      real*8, parameter :: sqrpi = 1.d0 / sqr2pi
      real*8, parameter :: sqrt2 = sqrt(2.d0)

c     these next two are used for comparison, not calculation, they do not
c     need to be as precise as the above

      real*8, parameter :: thrsh = 0.67448975d0
c$$$      real*8, parameter :: root32 = 5.65685424949238058d0
      real*8, parameter :: root32 = sqrt(32.d0)
