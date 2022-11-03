program test

  use type_def, only: &
       real64, int64

  use fitting, only: &
       chi2

  implicit none

  ! integer(kind=int64), parameter :: &
  !      nstar = 3, &
  !      nel = 9
  ! integer(kind=int64), parameter :: &
  !      nstar = 1, &
  !      nel = 1, &
  !      ncov = 0
  ! real(kind=real64), dimension(nel) :: &
  !      obs, err
  ! real(kind=real64), dimension(nel, ncov) :: &
  !      cov
  ! real(kind=real64), dimension(nstar, nel) :: &
  !      abu
  ! real(kind=real64), dimension(nstar) :: &
  !      c, f1, f2, f1up, f1do, f2up, f2do, cup, cdo
  ! real(kind=real64) :: &
  !      f, fup, fdo, h
  ! integer(kind=int64) :: &
  !      i

    ! obs = (/ -5.86418771d0, -6.19418771d0, -6.23418771d0, &
    ! -9.26418771d0, -8.39418771d0, -9.48418771d0, -10.60418771d0, &
    ! -12.36418771d0, -10.20418771d0 /)
    ! err = (/0.24d0, 0.3d0, 0.24d0, 0.08d0, 0.07d0, &
    ! 0.09d0, 0.2d0, 0.2d0, 0.15d0 /)
    !
    ! abu(1,:) = (/ 9.23837878d-04, 2.86915484d-08, 9.86329569d-04, &
    ! 7.22169530d-07, 1.44506086d-05, 2.98930995d-07, 3.79484481d-13, &
    ! 3.49929099d-15, 9.71069640d-18 /)
    ! abu(2,:) = (/ 1.29517519d-03, 3.96547937d-07, 7.74776042d-03, &
    ! 5.80834667d-06, 1.79736058d-04, 2.76109507d-06, 1.32311697d-05, &
    ! 8.84902251d-08, 1.15314210d-05 /)
    ! abu(3,:) = (/ 1.06639960d-06, 1.63124275d-05, 3.58418401d-07, &
    ! 8.92132488d-12, 1.22640917d-11, 2.84983024d-13, 8.85979959d-13, &
    ! 5.86733432d-16, 1.37888246d-17 /)
    !
    ! c(:) = 1.0d-4
    !
    ! obs = (/ -5.0d0 /)
    ! err = (/ -0.5d0 /)
    ! abu(1,:) = (/ 1.0d-05 /)
    !
    ! c = (/ 2.0d0/)
    !
    ! h = 1.0d-7
    ! cdo = c
    ! cup = c
    ! cdo(1) = cdo(1) - h
    ! cup(1) = cup(1) + h

    ! call chi2(f, f1, f2, c, obs, err, abu, nstar, nel, 1)
    ! call chi2(fdo, f1do, f2do, cdo, obs, err, abu, nstar, nel, 1)
    ! call chi2(fup, f1up, f2up, cup, obs, err, abu, nstar, nel, 1)

    ! print*,'f:'
    ! print*,f
    ! print*,fdo
    ! print*,fup

    ! print*,f1(1)
    ! print*,(fup - fdo)/(2*h)
    !
    ! print*,f2(1)
    ! print*,(fup - 2*f + fdo)/h**2
    ! print*,(f1up - f1do)/(2*h)
    ! print*,f1
    ! print*,f2

    ! call newton(c, obs, err, abu, nstar, nel, 1)

  print *, 'Hello world!'

end program test
