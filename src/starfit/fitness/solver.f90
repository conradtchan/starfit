module solver

  use type_def, only: &
       real64, int64

  implicit none

  private

  public :: &
       chi2, newton, fitness, &
       prime, analyticsolve, singlesolve, psolve

contains

  ! Due to poor f2py / numpy.distutils code we need to carry nel
  ! rather than taking it from star_data

  subroutine chi2(f, c, abu, nel, nstar)

    use star_data, only: &
         icdf, obs, err, det, measu

    ! use star_data, only: &
    !      diff_covariance, iuncor, iupper, nupper

    use norm, only: &
         logcdf, logcdfp

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10, &
         ln10i2 = 2.d0 / ln10, &
         ln10i2m = -ln10i2

    integer(kind=int64), intent(in) :: &
         nel, nstar
    real(kind=real64), intent(in), dimension(nstar)  :: &
         c
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(out) :: &
         f

    real(kind=real64), dimension(nel) :: &
         logy, diff, diffe, diffe2, erri, &
         ddif, ddife, ddife2
    integer(kind=int64) :: &
         i, j

    logy(:) = 0.d0
    do i = 1, nstar
       do j = 1, nel
          logy(j) = logy(j) + c(i) * abu(i,j)
       enddo
    enddo
    logy(:) = log(logy(:)) * ln10i

    diff(:) = logy(:) - obs(:)

    erri(:) = 1.d0 / err(:)
    diffe(:)  = diff(:) * erri(:)
    diffe2(:) = diffe(:)**2

    ddif(:) = logy(:) - det(:)
    ddife(:)  = ddif(:) * erri(:)
    ddife2(:) = ddife(:)**2

    f = 0.d0
    if (icdf == 0) then
       do i = 1, nel
          if (measu(i)) then
             f = f + diffe2(i)
             if (ddif(i) < 0.d0) then
                f = f - ddife2(i)
             endif
          else
             if (diff(i) > 0.d0) then
                f = f + diffe2(i)
             endif
          endif
       enddo
    else
       do i = 1, nel
          if (measu(i)) then
             f = f + diffe2(i) + 2.d0*logcdf(ddife(i))
          else
             f = f - 2.d0*logcdf(diffe(i))
          endif
       enddo
    endif
    ! do i = 1, nel
    !    if (measu(i)) then
    !       f = f + diffe2(i)
    !       if (icdf == 1) then
    !          f = f + 2.d0*logcdf(ddife(i))
    !       else
    !          if (ddif(i) < 0.d0) then
    !             f = f - ddife2(i)
    !          endif
    !       endif
    !    else
    !       if (icdf == 1) then
    !          f = f - 2.d0*logcdf(diffe(i))
    !       else
    !          if (diff(i) > 0.d0) then
    !             f = f + diffe2(i)
    !          endif
    !       endif
    !    endif
    ! enddo

    ! erri(:) = 1.d0 / err(:)
    ! f = diff_covariance(diff)
    ! f = f + sum((diff(iuncor) * erri(iuncor))**2)
    ! if (icdf == 1) then
    !    f = f - 2.d0* sum(logcdf(diff(iupper) * erri(iupper)))
    !    ! do i=1, nupper
    !    !    j = iupper(i)
    !    !    f = f - 2.d0*logcdf(diff(j) * erri(j))
    !    ! enddo
    !    f = f + 2.d0* sum(logcdf(ddif(imeasu) * erri(imeasu)))
    !    ! do i=1, nmeasu
    !    !    j = imeasu(i)
    !    !    f = f + 2.d0*logcdf(ddif(j) * erri(j))
    !    ! enddo
    ! else
    !    do i=1, nupper
    !       j = iupper(i)
    !       if (diff(j) > 0.d0) then
    !          f = f + (diff(j) * erri(j))**2
    !       endif
    !    enddo
    !    ! f = f + sum((min(diff(iupper), 0.d0) * erri(iupper))**2)
    !    do i=1, nmeasu
    !       j = imeasu(i)
    !       if (ddif(j) < 0.d0) then
    !          f = f - (ddif(j) * erri(j))**2
    !       endif
    !    enddo
    !    ! f = f - sum((max(ddff(imeasu), 0.d0) * erri(imeasu))**2)
    ! endif

  end subroutine chi2


  subroutine chi2_prime(f1, f2, c, abu, nel, nstar)

    ! return first and second derivative of chi2

    use star_data, only: &
         icdf, obs, err, det, measu

    use norm, only: &
         logcdf, logcdfp

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10, &
         ln10i2 = 2.d0 / ln10, &
         ln10i2m = -ln10i2

    integer(kind=int64), intent(in) :: &
         nel, nstar
    real(kind=real64), intent(in), dimension(nstar)  :: &
         c
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(out), dimension(nstar) :: &
         f1
    real(kind=real64), intent(out), dimension(nstar, nstar) :: &
         f2

    real(kind=real64), dimension(nel) :: &
         y, logy, yi, yi2, diff, diffe, erri, erri2, &
         ddif, ddife
    real(kind=real64) :: &
         ! fd1, fd2, &
         fa, fi1, fi2
    integer(kind=int64) :: &
         i, j, k

    y(:) = 0.d0
    do i = 1, nel
       do j = 1, nstar
          y(i) = y(i) + c(j) * abu(j,i)
       enddo
    enddo
    logy(:) = log(y(:)) * ln10i

    erri(:) = 1.d0 / err(:)
    erri2(:) = erri(:)**2

    diff(:) = logy(:) - obs(:)
    diffe(:)  = diff(:) * erri(:)

    ddif(:) = logy(:) - det(:)
    ddife(:)  = ddif(:) * erri(:)

    ! This needs to be checked, the sums may be off

    yi(:) = 1.d0 / yi(i)
    yi2(:) = yi(:)**2

    f1(:)   = 0.d0
    f2(:,:) = 0.d0

    if (icdf == 0) then
       do i = 1, nel
          fi1 = erri(i) * yi(i)
          fi2 = yi2(i) * (erri2(i) - erri(i) * diff(i))
          if (measu(i)) then
             do j=1, nel
                fa = fi1 * abu(j, i)
                f1(j) = f1(j) + diffe(i) * fa
                if (ddif(i) < 0.d0) then
                   f1(j) = f1(j) - ddife(i) * fa
                endif
                do k=1, nel
                   ! this is where NR may fail, discontinuous change in derivative
                   if (ddif(i) > 0.d0) then
                      f2(j,k) = f2(j,k) + abu(j,i) * abu(k,i) * fi2
                   endif
                enddo
             enddo
          else
             ! if (diff(i) > 0.d0) then
             !    f = f + diffe2(i)
             ! endif
          endif
       enddo
    else
       ! do i = 1, nel
       !    if (measu(i)) then
       !       f = f + diffe2(i)
       !       f = f + 2.d0*logcdf(ddife(i))
       !    else
       !       f = f - 2.d0*logcdf(diffe(i))
       !    endif
       ! enddo
    endif


    ! do i = 1, nstar
    !    do j = 1, nel
    !       if (measu(j)) then
    !          f1(i) = f1(i) + diffe(j) * erri(j) * abu(i, j) * yi(j)
    !          ! here we actully need to comoute the Jacobian!
    !          do k=1, nstar
    !          f2(i) = f2(i) + (abu(i, j) * errsumi(j))**2 * (ln10i - diff(j))
    !          if (icdf == 1) then
    !             fd1 = -logcdfp(ddife(j))
    !             fd2 = fd1 * (fd1 + ddife(j))
    !             f1(i) = f1(i) - abu(i, j) * fd1 * errsumi(j)
    !             f2(i) = f2(i) + erri(j) * (abu(i, j) / summed(j))**2 * (fd2/(-ln10 * err(j)) + fd1)
    !          else
    !             if (ddif(i) > 0.d0) then
    !                f1(i) = f1(i) + (ddife(j) * abu(i, j)) * errsumi(j)
    !                f2(i) = f2(i) - (abu(i, j) * errsumi(j))**2 * (ln10i - ddif(j))
    !             endif
    !          endif
    !       else
    !          if (icdf == 1) then
    !             fd1 = logcdfp(diffe(j))
    !             fd2 = -fd1 * (fd1 + diffe(j))
    !             f1(i) = f1(i) - abu(i, j) * fd1 * errsumi(j)
    !             f2(i) = f2(i) + erri(j) * (abu(i, j) / summed(j))**2 * (fd2/(-ln10 * err(j)) + fd1)
    !          else
    !             if (diff(i) < 0.d0) then
    !                f1(i) = f1(i) - (diffe(j) * abu(i, j)) * errsumi(j)
    !                f2(i) = f2(i) + (abu(i, j) * errsumi(j))**2 * (ln10i - diff(j))
    !             endif
    !          endif
    !       endif
    !    enddo
    ! enddo

    f1(:) = f1(:) * ln10i2m
    f2(:,:) = f2(:,:) * ln10i2

  end subroutine chi2_prime


  subroutine newton(c, abu, nel, nstar)

    implicit none

    integer(kind=int64), parameter :: &
         trials = 20

    integer(kind=int64), intent(in) :: &
         nel, nstar
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu

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
         delta, deltaold, f1, randoms
    real(kind=real64), dimension(nstar, nstar) :: &
         f2

    integer(kind=int64) :: &
         i, j, jstop

    !Initial N-R values

    delta(:) = 0
    maxstep = 1.0d0
    fold = 1.d6
    f = 1.d3

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
       call chi2(f, c, abu, nel, nstar)
       do i = 1, 100

          ! Evaluate function and derivatives

          call chi2_prime(f1, f2, c, abu, nel, nstar)

          ! If any NaNs are present, redo step with half the step size

          if (any(isnan(delta))) then
             delta(:) = 0.50d0*deltaold(:)
             deltaold(:) = delta(:)
             c(:) = cold(:)
             a(:) = log(c(:))
          else

             ! Calculate the step

             deltaold(:) = delta(:)

             ! BUG - replace by LEQS
             ! delta(:) = (f1(:) / (f1(:) + f2(:)*c(:)))
          endif

          ! Calculate the length of the step

          deltalen = sqrt(dot_product(delta,delta))

          ! Calculate how much the step needs to be scaled back

          steplimit = min(0.5d0, maxstep/deltalen)

          ! Apply step
          a(:) = a(:) - steplimit*delta(:)
          c(:) = exp(a(:))

          ! Convergence detection

          fold = f
          call chi2(f, c, abu, nel, nstar)

          ! print*,i, '|', f, '|', a, '|', c, '|', f1, '|', f2, '|', delta, '|', deltalen, maxstep, steplimit

          if ((fold - f) < 1.d-6) then
             exit
          end if
       end do

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


  subroutine psolve(c, abu, nel, nstar)

    use abu_data, only: &
         set_abu_data

    use powell, only: &
         uobyqa

    implicit none

    integer(kind=int64), intent(in) :: &
         nel, nstar
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         x
    real(kind=real64) :: rhobeg, rhoend
    integer(kind=int64) :: calls, iprint

    ! save module data for calfun

    ! an pointer to an allocated data structure should be used instead.

    call set_abu_data(abu, nel, nstar)

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


  subroutine prime(x, f1, f2, abu, nel)

    ! Returns the first and second derivatives of (1/2)*chi^2 with
    ! respect to x for single star fits

    use star_data, only: &
         icdf, obs, err, det, measu

    use norm, only: &
         logcdfp

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10

    real(kind=real64), intent(in) :: &
         x
    real(kind=real64), intent(in), dimension(nel) :: &
         abu
    integer(kind=int64), intent(in) :: &
         nel

    real(kind=real64), intent(out) :: &
         f1, f2

    real(kind=real64), dimension(nel) :: &
         ei, ei2
    real(kind=real64) :: &
         diff, de, df, logabu

    integer(kind=int64) :: &
         i

    f1 = 0.d0
    f2 = 0.d0

    if (icdf == 0) then
       ei2 = 1.d0 / err**2
       do i = 1, nel
          logabu = log(abu(i)) * ln10i

          ! Upper limits

          diff = logabu - obs(i) + x
          if (measu(i) .or. (diff > 0.d0)) then
             f1 = f1 + diff * ei2(i)
             f2 = f2 + ei2(i)
          endif

          ! Detection thresholds

          if (measu(i)) then
             diff = logabu - det(i) + x
             if (diff < 0.d0) then
                f1 = f1 - diff * ei2(i)
                f2 = f2 - ei2(i)
             endif
          endif
       enddo
    else
       ei = 1.d0 / err
       ei2 = ei**2
       do i = 1, nel
          logabu = log(abu(i)) * ln10i

          ! Upper limits if error < 0 (NOTE: changes sign of diff)

          diff = (logabu + x - obs(i)) * ei(i)
          if (measu(i))  then
             f1 = f1 + diff * ei(i)
             f2 = f2 + ei2(i)
          else
             df = - logcdfp(diff) * ei(i)
             de = - diff * ei(i)
             f1 = f1 + df
             f2 = f2 + df * (df + de)
          endif

          ! Detection thresholds

          if (measu(i)) then
             diff = (logabu + x - det(i)) * ei(i)
             df = + logcdfp(diff) * ei(i)
             de = - diff * ei(i)
             f1 = f1 + df
             f2 = f2 + df * (df + de)
          endif
       enddo
    endif
  end subroutine prime


  subroutine singlesolve(c, abu, nel)

    ! single NR solver

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10

    real(kind=real64), intent(in), dimension(nel) :: &
         abu
    integer(kind=int64), intent(in) :: &
         nel

    real(kind=real64), intent(inout) :: &
         c

    real(kind=real64) :: &
         x, f1, f2, delta

    integer(kind=int64) :: &
         i

    x = log(c) * ln10i
    do i = 1, 20
       call prime(x, f1, f2, abu, nel)
       if (f2 == 0.d0) exit
       delta = f1 / f2
       x = x - delta
       if (abs(delta) < 1.0d-6) then
          c = exp(x * ln10)
          exit
       endif
    enddo

  end subroutine singlesolve


  subroutine analyticsolve(c, abu, nel)

    use star_data, only: &
         obs, err, det, measu, upper

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10

    real(kind=real64), intent(in), dimension(nel) :: &
         abu
    real(kind=real64), intent(out) :: &
         c
    integer(kind=int64), intent(in) :: &
         nel

    real(kind=real64) :: &
         x, diff, ei2s
    real(kind=real64), dimension(nel)  :: &
         ei2, logabu, diff_obs, diff_det

    logical, dimension(nel) :: &
         include_obs, include_det
    logical :: &
         change

    integer(kind=int64) :: &
         i

    logabu(:) = log(abu(:)) * ln10i

    diff_obs(:) = logabu - obs(:)
    diff_det(:) = logabu - det(:)

    ! First, include all upper limits as if they were detections, and
    ! exclude all detection thresholds as it you were above

    ei2(:) = 1.d0 / err(:)**2

    change = .true.

    include_obs(:) = measu(:)
    include_det(:) = .false.

    do while (change)

       ! (re)calculate with ignored elements

       ei2s = 0.d0
       diff = 0.d0
       do i = 1, nel
          if (include_obs(i)) then
             ei2s = ei2s + ei2(i)
             diff = diff + diff_obs(i) * ei2(i)
          endif
          if (include_det(i)) then
             diff = diff - diff_det(i) * ei2(i)
          endif
       enddo

       x = -diff / ei2s

       change = .false.
       do i = 1, nel
          if (upper(i)) then
             diff = diff_obs(i) + x
             if (diff < 0.0d0) then
                include_obs(i) = .false.
                change = .true.
             endif
          else
             diff = diff_det(i) + x
             if (diff < 0.0d0) then
                include_det(i) = .true.
                change = .true.
             endif
          endif
       enddo
    enddo

    c = exp(x * ln10)

  end subroutine analyticsolve


  subroutine fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf)

    use star_data, only: &
         set_star_data, abu_covariance

    use type_def, only: &
         real64, int64

    implicit none

    integer(kind=int64), parameter :: &
         isolve = 2

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
         icdf

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
    ! otherwise keep relative weights fixed

    call set_star_data(obs, err, det, cov, nel, ncov, icdf)

    if (ls == 1) then
       do i = 1, nsol
          if (nstar == 1) then
             if (icdf == 0) then
                call analyticsolve(c(i,1), abu(i,1,:), nel)
             else
                call singlesolve(c(i,1), abu(i,1,:), nel)
             endif
          else
             if (isolve == 1) then
                ! Dodgy NR solver
                call newton(c(i,:), abu(i,:,:), nel, nstar)
             else
                ! Slower, but more reliable UOBYQA solver
                call psolve(c(i,:), abu(i,:,:), nel, nstar)
             endif
          endif
          scale = sum(c(i,:))
          if (scale > 1.d0) then
             c(i,:) = c(i,:) * (1.d0 / scale)
          endif
       end do
    else if (ls == 0) then
       do i = 1, nsol
          scale = 1.d0
          summed(:) = 0.d0
          do j = 1, nstar
             do k = 1, nel
                summed(k) = summed(k) + c(i,j) * abu(i,j,k)
             enddo
          enddo
          if (icdf == 0) then
             call analyticsolve(scale, summed, nel)
          else
             call singlesolve(scale, summed, nel)
          endif
          c(i,:) = c(i,:) * min(scale, 1.d0 / sum(c(i,:)))
       enddo
    else if (ls /= 2) then
       error stop 'invalid local search option'
    else
    endif

    ! get final fit value for possibly adjusted scales

    do i = 1, nsol
       call chi2(f(i), c(i,:), abu(i,:,:), nel, nstar)
    enddo

  end subroutine fitness

end module solver

! =======================================================================


subroutine calfun(nstar, x, f)

  use type_def, only: &
       int64, real64

  use abu_data, only: &
       abu

  use star_data, only: &
       nel

  use solver, only: &
       chi2

  implicit none

  integer(kind=int64), intent(in) :: &
       nstar
  real(kind=real64), intent(in), dimension(nstar) :: &
       x

  real(kind=real64), intent(out) :: &
       f

  real(kind=real64), dimension(nstar) :: &
       tanhx
  integer(kind=int64) :: &
       i

  ! Transform from -inf/inf to 0/1
  ! x is used for the solver, tanhx is physically meaningful

  tanhx = 0.5d0 * (1.0d0 + tanh(x))

  ! Calculate the chi2

  call chi2(f, tanhx, abu(1:nstar,1:nel), nel, nstar)

  ! Build a wall at zero

  do i = 1, nstar
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
       chi2, newton

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
