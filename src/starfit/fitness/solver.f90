! interface module for f2py

! include and wrap only what is actually needed

subroutine fitness_(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

  use fitting, only: &
       fitness

  use typedef, only: &
       real64, int64

  implicit none

  integer(kind=int64), intent(in) :: &
       nstar, nel, nsol, ncov

  real(kind=real64), dimension(nel), intent(in) :: &
       obs, err, det
  real(kind=real64), dimension(nel, ncov), intent(in) :: &
       cov
  real(kind=real64), dimension(nsol, nstar, nel), intent(in) :: &
       abu
  integer(kind=int64), intent(in) :: &
       ls
  integer(kind=int64), intent(in) :: &
       icdf, flags

  !f2py real(kind=real64), intent(out), dimension(nsol) :: f
  real(kind=real64), dimension(nsol), intent(out) :: &
       f
  !f2py real(kind=real64), intent(in,out), dimension(nsol, nstar) :: c
  real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
       c

  call fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

end subroutine fitness_
