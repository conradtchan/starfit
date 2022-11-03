module fitting

  use type_def, only: &
       real64, int64

  implicit none

  private

  public :: &
       fitness, chi2

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

    ! ! A potntiall faster version (to be timed)
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
         fa, fb, fc, fi1, fi2, fi3
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

    yi(:) = 1.d0 / yi(i)
    yi2(:) = yi(:)**2

    f1(:)   = 0.d0
    f2(:,:) = 0.d0

    ! This needs to be checked

    if (icdf == 0) then
       do i = 1, nel
          if (measu(i)) then
             fi1 = erri(i) * yi(i)
             fi2 = yi2(i) * (erri2(i) - erri(i) * diff(i))
             if (ddif(i) > 0.d0) then
                fi3 = yi2(i) * (erri2(i) - erri(i) * ddif(i))
             endif
             do j=1, nel
                fa = fi1 * abu(j, i)
                f1(j) = f1(j) + diffe(i) * fa
                if (ddif(i) < 0.d0) then
                   f1(j) = f1(j) - ddife(i) * fa
                endif
                do k=1, nel
                   fb = abu(j,i) * abu(k,i)
                   f2(j,k) = f2(j,k) + fb * fi2
                   if (ddif(i) > 0.d0) then
                      f2(j,k) = f2(j,k) - fb * fi3
                   endif
                enddo
             enddo
          else if (diff(i) > 0.d0) then
             fi1 = erri(i) * yi(i)
             fi2 = yi2(i) * (erri2(i) - erri(i) * diff(i))
             do j=1, nel
                fa = fi1 * abu(j, i)
                f1(j) = f1(j) - diffe(i) * fa
                do k=1, nel
                   fb = abu(j,i) * abu(k,i)
                   f2(j,k) = f2(j,k) + fb * fi2
                enddo
             enddo
          endif
       enddo
    else
       do i = 1, nel
          fi1 = erri(i) * yi(i)
          fc = logcdfp(ddife(i))
          if (measu(i)) then
             fi2 = yi2(i) * (erri2(i) - erri(i) * diff(i))
             do j=1, nel
                fa = fi1 * abu(j, i)
                f1(j) = f1(j) + diffe(i) * fa
                f1(j) = f1(j) + fc * abu(j,i) * yi(i)
                do k=1, nel
                   fb = abu(j,i) * abu(k,i)
                   f2(j,k) = f2(j,k) + fb * fi2
                   f2(j,k) = f2(j,k) - fc * fb * yi(i) * &
                        (1 + (ddife(i) + fc) * yi(i))
                enddo
             enddo
          else
             fi1 = erri(i) * yi(i)
             fc = logcdfp(diffe(i))
             do j=1, nel
                fa = fi1 * abu(j, i)
                f1(j) = f1(j) + fc * abu(j,i) * yi(i)
                do k=1, nel
                   fb = abu(j,i) * abu(k,i)
                   f2(j,k) = f2(j,k) + fc * fb * yi(i) * &
                        (1 + (diffe(i) + fc) * yi(i))
                enddo
             enddo
          endif
        enddo
    endif

    f1(:) = f1(:) * ln10i2
    f2(:,:) = f2(:,:) * ln10i2

  end subroutine chi2_prime


  subroutine newton(c, abu, nel, nstar)

    use mleqs, only: &
         leqs

    implicit none

    integer(kind=int64), parameter :: &
         max_steps = 100

    integer(kind=int64), intent(in) :: &
         nel, nstar
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         dc, f1
    real(kind=real64), dimension(nstar, nstar) :: &
         f2
    real(kind=real64) :: &
         dcr

    integer(kind=int64) :: &
         i

    !Initial N-R values

    c(:) = max(min(c(:), 1.d0), 1e-8)

    do i = 1, max_steps
       call chi2_prime(f1, f2, c, abu, nel, nstar)
       dc(:) = leqs(f2, f1, nstar)
       dcr = maxval(abs(dc) / dc)
       if (dcr > 0.5) then
          dc(:) = dc(:) * (0.5d0 / dcr)
       endif
       if (dcr < 1.0d-6) then
          exit
       endif
    enddo

    if (i >= max_steps) then
       error stop '[newton] did not converge.'
    endif

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

    ! should not actually be used for icdf == 0

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

    ! should not be used for icdf == 0

    implicit none

    real(kind=real64), parameter :: &
         ln10 = log(10.0d0), &
         ln10i = 1.d0 / ln10
    integer(kind=int64), parameter :: &
         max_steps = 20


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
    do i = 1, max_steps
       call prime(x, f1, f2, abu, nel)
       if (f2 == 0.d0) exit
       delta = f1 / f2
       x = x - delta
       if (abs(delta) < 1.0d-6) then
          c = exp(x * ln10)
          exit
       endif
    enddo

    if (i >= max_steps) then
       error stop '[singlesolve] did not converge.'
    endif

  end subroutine singlesolve


  subroutine analyticsolve(c, abu, nel)

    use star_data, only: &
         obs, err, det, upper

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

    include_obs(:) = .true.
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
             if (include_obs(i)) then
                diff = diff_obs(i) + x
                if (diff < 0.0d0) then
                   include_obs(i) = .false.
                   change = .true.
                endif
             endif
          else
             if (.not.include_det(i)) then
                diff = diff_det(i) + x
                if (diff < 0.0d0) then
                   include_det(i) = .true.
                   change = .true.
                endif
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
         y

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

             ! This is what it should be when chi2_prime has been tested
             ! if (icdf == 0) then
             !    call psolve(c(i,:), abu(i,:,:), nel, nstar)
             ! else
             !    call newton(c(i,:), abu(i,:,:), nel, nstar)
             ! endif

             if (isolve == 1) then
                ! Dodgy NR solver
                ! should work for cdf == 1 but may fail for cdf==0

                call newton(c(i,:), abu(i,:,:), nel, nstar)
             else
                ! Slower UOBYQA solver that does not depend on C(2)
                ! function for chi2

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
          y(:) = 0.d0
          do j = 1, nstar
             do k = 1, nel
                y(k) = y(k) + c(i,j) * abu(i,j,k)
             enddo
          enddo
          if (icdf == 0) then
             call analyticsolve(scale, y, nel)
          else
             call singlesolve(scale, y, nel)
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

end module fitting

! =======================================================================

subroutine calfun(nstar, x, f)

  use type_def, only: &
       int64, real64

  use abu_data, only: &
       abu, nel

  use fitting, only: &
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

  tanhx(:) = 0.5d0 * (1.0d0 + tanh(x(:)))

  ! Calculate the chi2

  call chi2(f, tanhx(1:nstar), abu(1:nstar,1:nel), nel, nstar)

  ! Build a wall at zero

  do i = 1, nstar
     if (abs(x(i)) > 12.d0) then
        f = f * exp((abs(x(i)) - 12.d0)**2)
     endif
  enddo

end subroutine calfun
