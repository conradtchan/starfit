module callfun_data

  use type_def, only: &
       real64, int64

  implicit none

  save

  integer(kind=int64) :: &
       nel_, nstar_
  real(kind=real64), dimension(:), allocatable :: &
       obs_, err_
  real(kind=real64), dimension(:, :), allocatable :: &
       abu_
  integer(kind=int64) :: &
       icdf_

contains

  subroutine set_data(nel, nstar, obs, err, abu, icdf)

    implicit none

    integer(kind=int64), intent(in) :: &
         nel, nstar, icdf
    real(kind=real64), dimension(:), intent(in) :: &
         obs, err
    real(kind=real64), dimension(:,:), intent(in) :: &
         abu

    if (allocated(abu_)) then
       deallocate(abu_, obs_, err_)
    endif

    if (.not. (size(obs, 1) == nel)) error stop 'obs dimeension mismatch'
    if (.not. (size(err, 1) == nel)) error stop 'err dimeension mismatch'
    if (.not. (size(abu, 2) == nel)) error stop 'abu dimeension 2 mismatch'
    if (.not. (size(abu, 1) == nstar)) error stop 'abu dimeension 1 mismatch'

    nel_ = nel
    nstar_ = nstar
    icdf_ = icdf

    ! allocate implicitly

    obs_ = obs
    err_ = err
    abu_ = abu

  end subroutine set_data

end module callfun_data



module solver

  use type_def, only: &
       real64, int64

  implicit none

  private

  public :: &
       chisq, newton, fitness, &
       prime, analyticsolve, singlesolve, psolve

contains

  subroutine chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf, ider)

    use norm, only: &
         logcdf, logcdfp

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10, &
         ln10i2 = 2.d0 / ln10, &
         ln10i2m = -ln10i2

    integer(kind=int64), intent(in) :: &
         nstar, nel, icdf
    real(kind=real64), intent(in), dimension(nstar)  :: &
         c
    real(kind=real64), intent(in), dimension(nel) :: &
         obs, err
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu
    integer(kind=int64), optional :: &
         ider
    !f2py integer(kind=int64),optional,intent(in) :: ider=0

    real(kind=real64), intent(out) :: &
         f
    real(kind=real64), intent(out), dimension(nstar) :: &
         f1, f2

    real(kind=real64), dimension(nel) :: &
         summed, diff, diffe, diffe2, errsum, erri, errsumi
    real(kind=real64) :: &
         fd1, fd2
    integer(kind=int64) :: &
         i, j
    logical :: &
         lder

    summed(:) = 0.d0
    do i = 1, nstar
        do j = 1, nel
            summed(j) = summed(j) + c(i) * abu(i,j)
        enddo
    enddo

    diff(:) = obs(:) - log10(summed(:))
    erri(:) = 1.d0 / abs(err(:))
    diffe(:) = diff(:) * erri(:)
    diffe2(:) = diffe(:)**2

    f = 0.d0
    do i = 1, nel
        if (err(i) > 0.d0) then
            f = f + diffe2(i)
        else
            if (icdf == 1) then
                f = f - 2.d0*logcdf(diffe(i))
            else
                if (diff(i) < 0.d0) then
                    f = f + diffe2(i)
                endif
            endif
        endif
    enddo

    ! Enable derivatives for using the NR solver

    if (present(ider)) then
       lder = ider == 1
    else
       lder = .False.
    endif

    if (lder) then
       f1(:) = 0.d0
       f2(:) = 0.d0

       ! This needs to be checked, the sums seem off

       errsum(:) = err(:) * summed(:)
       errsumi(:) = 1.d0 / errsum(:)

       do i = 1, nstar
          do j = 1, nel
             if (err(j) > 0.d0) then
                f1(i) = f1(i) + (diffe(j) * abu(i, j)) * errsum(j)
                f2(i) = f2(i) + (abu(i, j) * errsumi(j))**2 * (ln10i + diff(j))
             else
                if (icdf == 1) then
                   fd1 = logcdfp(diffe(j))
                   fd2 = -fd1 * (fd1 + diffe(j))
                   f1(i) = f1(i) - abu(i, j) * fd1 * errsumi(j)
                   f2(i) = f2(i) - erri(j) * (abu(i, j) / summed(j))**2 * (fd2/(ln10 * abs(err(j))) + fd1)
                else
                   if (diff(i) < 0) then
                      f1(i) = f1(i) + (diffe(j) * abu(i, j)) * errsumi(j)
                      f2(i) = f2(i) + (abu(i, j) * errsumi(j))**2 * (ln10i + diff(j))
                   endif
                endif
             endif
          enddo
       enddo

       f1(:) = f1(:) * ln10i2m
       f2(:) = f2(:) * ln10i2
    endif

  end subroutine chisq


  subroutine newton(c, obs, err, abu, nstar, nel, icdf)

    implicit none

    integer(kind=int64), parameter :: &
         trials = 20

    integer(kind=int64), intent(in) :: &
         nstar, nel
    real(kind=real64), dimension(nel), intent(in) :: &
         obs, err
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu
    integer(kind=int64), intent(in) :: &
         icdf

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c
    real(kind=real64), dimension(nstar) :: &
         cold
    real(kind=real64), dimension(trials, nstar) :: &
         starts, recordc
    real(kind=real64), dimension(trials) :: &
         recordf
    real(kind=real64), dimension(nstar) :: &
         a

    real(kind=real64) :: &
         f, fold, &
         maxstep, steplimit, deltalen
    real(kind=real64), dimension(nstar) :: &
         delta, deltaold, f1, f2, randoms

    integer(kind=int64) :: &
         i, j, jstop, ider, inoder

    !Initial N-R values

    delta(:) = 0
    maxstep = 1.0d0
    fold = 1.d6
    f = 1.d3
    ider = 1
    inoder = 0

    starts(1,:) = c(:)
    do j = 1, trials

       ! Find random starting points

       if (j > 1) then
          call random_number(randoms)
          starts(j,:) = 10.0d0**(-5.0d0*randoms + 1.0d0)
       endif

       c(:) = starts(j,:)
       cold(:) = c(:)
       a(:) = log(c(:))
       do i = 1, 100
          fold = f

          ! Evaluate function and derivatives

          call chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf, ider)

          ! If any NaNs are present, redo step with half the step size

          if (any(isnan(delta))) then
             delta(:) = 0.50d0*deltaold(:)
             deltaold(:) = delta(:)
             c(:) = cold(:)
             a(:) = log(c(:))
          else

             ! Calculate the step

             deltaold(:) = delta(:)
             delta(:) = (f1(:) / (f1(:) + f2(:)*c(:)))
          endif

          ! Calculate the length of the step

          deltalen = sqrt(dot_product(delta,delta))

          ! Calculate how much the step needs to be scaled back

          steplimit = min(0.5d0, maxstep/deltalen)

          ! Apply step
          ! print*,i, '|', f, '|', a, '|', c, '|', f1, '|', f2, '|', delta, '|', deltalen, maxstep, steplimit
          a(:) = a(:) - steplimit*delta(:)
          c(:) = exp(a(:))

          ! Convergence detection

          if ((fold - f) < 1.d-6) then
             exit
          end if
       end do

       call chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf, inoder)

       recordc(j,:) = c
       recordf(j) = f

       if (j > 1) then
          jstop = j
          if (recordf(j) - recordf(j-1) < 1.d-6) then
             exit
          endif
       endif
    enddo

    c = recordc(minloc(recordf(1:jstop),1),:)

  end subroutine newton


  subroutine psolve(c, obs, err, abu, nstar, nel, icdf)

    use callfun_data, only: &
         set_data

    use powell, only: &
         uobyqa

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar, nel
    real(kind=real64), intent(in), dimension(nel) :: &
         obs, err
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu
    integer(kind=int64), intent(in) :: &
         icdf

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         x
    real(kind=real64) :: rhobeg, rhoend
    integer(kind=int64) :: calls, iprint

    ! save module data for calfun

    ! an pointer to an allocated data structure should be used instead.

    call set_data(nel, nstar, obs, err, abu, icdf)

    ! Options for the uobyqa solver

    rhobeg = 1.d0
    rhoend = 1.d-5
    iprint = 0
    calls = 50 * nstar

    if (any(c >= 1.d0)) then
       print*, '[psolve] DEBUG IN: c = ', c
       error stop 'c >= 1'
    endif

    !Convert offsets to solver space
    x = atanh(c * 2.d0 - 1.d0)

    !Call solver
    call uobyqa(nstar, x, rhobeg, rhoend, iprint, calls)
    !Convert solver space to offsets
    c = 0.5d0 * (1.d0 + tanh(x))

  end subroutine psolve


  subroutine singlesolve(c, obs, err, abu, nel, icdf)

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10

    integer(kind=int64), intent(in) :: &
         nel
    real(kind=real64), intent(in), dimension(nel) :: &
         obs, err, abu
    integer(kind=int64), intent(in) :: &
         icdf

    real(kind=real64), intent(inout) :: &
         c

    real(kind=real64) :: &
         x, f, f1, f2, delta

    integer(kind=int64) :: &
         i

    x = log(c) * ln10i
    do i = 1, 20
       call prime(x, f1, f2, obs, err, abu, nel, icdf)
       if (f2 == 0.d0) exit
       delta = f1 / f2
       x = x - delta
       if (abs(delta) < 1.0d-6) then
          c = exp(x * ln10)
          exit
       endif
    enddo

  end subroutine singlesolve


  subroutine analyticsolve(c, obs, err, abu, nel)

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10

    integer(kind=int64), intent(in) :: &
         nel
    real(kind=real64), intent(out) :: &
         c
    real(kind=real64), intent(in), dimension(nel) :: &
         obs, err, abu

    real(kind=real64) :: &
         x, newdiff
    real(kind=real64), dimension(nel)  :: &
         ei2, diff, logabu

    logical, dimension(nel) :: &
         ignore
    logical :: &
         change

    integer(kind=int64) :: &
         i

    logabu(:) = log(abu(:)) * ln10i

    ! First, ignore upper limits

    do i = 1, nel
       ei2(i) = 1.d0 / err(i)**2
       diff(i) = (logabu(i) - obs(i)) * ei2(i)
    enddo
    x = -sum(diff) / sum(ei2)

    ! Check position

    change = .false.
    ignore(:) = .false.
    do i = 1, nel
       if (err(i) < 0.0d0) then
          newdiff = logabu(i) + x - obs(i)
          if (newdiff < 0.0d0) then
             ignore(i) = .true.
             change = .true.
          endif
       endif
    enddo

    ! Make adjustments until okay

    do while (change)

       ! Recalculate with ignored elements

       do i = 1, nel
          if (ignore(i)) then
             ei2(i) = 0.d0
             diff(i) = 0.d0
          else
             ei2(i) = 1.d0 / err(i)**2
             diff(i) = (logabu(i) - obs(i)) * ei2(i)
          endif
       enddo
       x = -sum(diff) / sum(ei2)

       ! Check position

       change = .false.
       do i = 1, nel
          if ((err(i) < 0.d0) .and. (.not. ignore(i))) then
             newdiff = logabu(i) + x - obs(i)
             if (newdiff < 0.d0) then
                ignore(i) = .true.
                change = .true.
             endif
          endif
       enddo
    enddo

    c = exp(x * ln10)

  end subroutine analyticsolve

  subroutine prime(x, f1, f2, obs, error, abu, nel, icdf)

    use norm, only: &
         logcdfp

    implicit none

    ! Returns the derivatives with respect to log(s)
    ! For single star fits only!

    integer(kind=int64), intent(in) :: &
         nel
    real(kind=real64), intent(in), dimension(nel) :: &
         abu, obs, error
    real(kind=real64), intent(in) :: &
         x
    integer(kind=int64), intent(in) :: &
         icdf

    real(kind=real64), intent(out) :: f1, f2

    real(kind=real64), dimension(nel) :: &
         ei, ei2
    real(kind=real64) :: &
         diff, de, df, e2

    integer(kind=int64) :: &
         i

    f1 = 0.d0
    f2 = 0.d0

    if (icdf == 0) then
       ei2 = 1.d0 / error**2
       do i = 1, nel
          diff = log10(abu(i)) + x - obs(i)

          ! Upper limits

          if ( (diff > 0.d0) .or. (error(i) > 0.d0) ) then
             f1 = f1 + diff * ei2(i)
             f2 = f2 + ei2(i)
          endif
       enddo
    else
       ei = 1.d0 / error
       ei2 = ei**2
       do i = 1, nel
          diff = (log10(abu(i)) + x - obs(i)) * ei(i)

          ! Upper limits if error < 0 (NOTE: changes sign of diff)

          if (error(i) > 0.d0)  then
             f1  = f1 + diff * ei(i)
             f2 = f2 + ei2(i)
          else
             df = - logcdfp(diff) * ei(i) !this has correct sign now
             de = - diff * ei(i)
             f1 = f1 + df
             f2 = f2 + df * (df + de)
          endif
       enddo
    endif
  end subroutine prime


  subroutine fitness(f, c, obs, err, corr, abu, nstar, nel, nsol, ncorr, ls, icdf)

    use type_def, only: &
         real64, int64

    ! use solver, only: &
    !      analyticsolve, singlesolve, psolve, chisq

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar, nel, nsol, ncorr

    real(kind=real64), dimension(nel), intent(in) :: &
         obs, err
    real(kind=real64), dimension(nel, ncorr), intent(in) :: &
         corr
    real(kind=real64), dimension(nsol, nstar, nel), intent(in) :: &
         abu
    logical, intent(in) :: &
         ls
    integer(kind=int64), intent(in) :: &
         icdf

    real(kind=real64), dimension(nstar) :: &
         f1, f2
    real(kind=real64) :: &
         scale
    real(kind=real64), dimension(nel) :: &
         summed

    integer(kind=int64) :: &
         i, j, k

    !f2py real(kind=real64), intent(out), dimension(nsol) :: f
    real(kind=real64), dimension(nsol), intent(out) :: &
         f
    !f2py real(kind=real64), intent(in,out), dimension(nsol, nstar) :: c
    real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
         c

    ! If localsearch is enabled, modify the offsets first

    if (ls .eqv. .true.) then
       do i = 1, nsol
          if (nstar == 1) then
             if (icdf == 0) then
                call analyticsolve(c(i,1), obs, err, abu(i,1,:), nel)
             else
                call singlesolve(c(i,1), obs, err, abu(i,1,:), nel, icdf)
             endif
          else
             !! Dodgy NR solver
             ! call newton(c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
             ! Slower, but more reliable UOBYQA solver
             !return
             call psolve(c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
          endif
          call chisq(f(i), f1, f2, c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
       end do
    else
       return
       do i = 1, nsol
          scale = 1.d0
          summed(:) = 0.d0
          do j = 1, nstar
             do k = 1, nel
                summed(k) = summed(k) + c(i,j) * abu(i,j,k)
             enddo
          enddo
          if (icdf == 0) then
             call analyticsolve(scale, obs, err, summed, nel)
          else
             call singlesolve(scale, obs, err, summed, nel, icdf)
          endif
          c(i,:) = c(i,:) * scale
          call chisq(f(i), f1, f2, c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
       enddo
    endif

  end subroutine fitness

end module solver

! =======================================================================


subroutine calfun(n, x, f)

  use type_def, only: &
       int64, real64

  use callfun_data, only: &
       obs_, err_, abu_, nel_, icdf_

  use solver, only: &
       chisq

  implicit none

  integer(kind=int64), intent(in) :: &
       n
  real(kind=real64), intent(in), dimension(n) :: &
       x

  real(kind=real64), intent(out) :: &
       f

  real(kind=real64), dimension(n) :: &
       tanhx, f1, f2
  integer(kind=int64) :: &
       i

  ! Transform from -inf/inf to 0/1
  ! x is used for the solver, tanhx is physically meaningful

  tanhx = 0.5d0 * (1.0d0 + tanh(x))

  ! Calculate the chisq

  call chisq(f, f1, f2, tanhx, obs_, err_, abu_(1:n,1:nel_), n, nel_, icdf_)

  ! Build a wall at zero

  do i = 1, n
     if (abs(x(i)) > 12.d0) then
        f = f * exp((abs(x(i)) - 12.d0)**2)
     endif
  enddo

end subroutine calfun


! =======================================================================

program test

  use type_def, only: &
       real64, int64

  use solver, only: &
       chisq, newton

  implicit none

  ! integer(kind=int64), parameter :: &
  !      nstar = 3, &
  !      nel = 9
  integer(kind=int64), parameter :: &
       nstar = 1, &
       nel = 1
  real(kind=real64), dimension(nel) :: &
       obs, err
  real(kind=real64), dimension(nstar, nel) :: &
       abu
  real(kind=real64), dimension(nstar) :: &
       c, f1, f2, f1up, f1do, f2up, f2do, cup, cdo
  real(kind=real64) :: &
       f, fup, fdo, h
  integer(kind=int64) :: &
       i

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

    ! call chisq(f, f1, f2, c, obs, err, abu, nstar, nel, 1)
    ! call chisq(fdo, f1do, f2do, cdo, obs, err, abu, nstar, nel, 1)
    ! call chisq(fup, f1up, f2up, cup, obs, err, abu, nstar, nel, 1)

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

end program test
