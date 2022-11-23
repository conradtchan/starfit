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


subroutine fitness_m_(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

  use fitting, only: &
       fitness_m

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

  !f2py real(kind=real64), intent(out), dimension(nsol, nel, nel) :: f
  real(kind=real64), dimension(nsol, nel, nel), intent(out) :: &
       f
  !f2py real(kind=real64), intent(in,out), dimension(nsol, nstar) :: c
  real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
       c

  call fitness_m(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

end subroutine fitness_m_


subroutine get_complete_matrix_(m, obs, err, det, cov, nel, ncov, icdf)

  use star_data, only: &
       set_star_data, get_complete_matrix

  use typedef, only: &
       real64, int64

  implicit none

  integer(kind=int64), intent(in) :: &
       nel, ncov

  real(kind=real64), dimension(nel), intent(in) :: &
       obs, err, det
  real(kind=real64), dimension(nel, ncov), intent(in) :: &
       cov
  integer(kind=int64), intent(in) :: &
       icdf

  !f2py real(kind=real64), intent(out), dimension(nel, nel) :: f
  real(kind=real64), dimension(nel,nel), intent(out) :: &
       m

  call set_star_data(obs, err, det, cov, nel, ncov, icdf)
  call get_complete_matrix(m)

end subroutine get_complete_matrix_


subroutine get_complete_inverse_(m1, obs, err, det, cov, nel, ncov, icdf)

  use star_data, only: &
       set_star_data, get_complete_inverse

  use typedef, only: &
       real64, int64

  implicit none

  integer(kind=int64), intent(in) :: &
       nel, ncov

  real(kind=real64), dimension(nel), intent(in) :: &
       obs, err, det
  real(kind=real64), dimension(nel, ncov), intent(in) :: &
       cov
  integer(kind=int64), intent(in) :: &
       icdf

  !f2py real(kind=real64), intent(out), dimension(nel, nel) :: f
  real(kind=real64), dimension(nel,nel), intent(out) :: &
       m1

  call set_star_data(obs, err, det, cov, nel, ncov, icdf)
  call get_complete_inverse(m1)

end subroutine get_complete_inverse_
