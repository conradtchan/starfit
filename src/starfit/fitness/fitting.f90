module fitting

  use typedef, only: &
       real64, int64

  implicit none

  logical :: &
       stop_on_nonconvergence = .false., &
       stop_on_large_offset = .false., &
       stop_on_zero_offset = .false.

  integer(kind=int64), parameter :: &
       FLAGS_LIMIT_SOLUTION_BIT = 0, &
       FLAGS_LIMITED_SOLVER_BIT = 1, &
       FLAGS_NO_CHI2 = 2


  private

  real(kind=real64), parameter :: &
       ln10 = log(10.0d0), &
       ln10i = 1.d0 / ln10, &
       ln10i2 = 2.d0 * ln10i, &
       ln10i2p2 = ln10i2**2, &
       ln10ip22 = 2.d0 * ln10i**2, &
       ln10i2m = -ln10i2, &
       ALMOST_ONE = 1.d0 - 1.d-14

  logical :: &
       wall_chi2_prime = .false.

  public :: &
       fitness, fitness_m

contains

  subroutine fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

    ! If localsearch is enabled, adjust all offsets first
    ! otherwise keep relative weights fixed
    ! If ls < 0, only return chi2, no optimisation
    ! If ls == 2, force use of psolve
    ! if ls == 3, use classic NR solver (not recommended)
    !
    ! flags:
    !    bit 0, "lsolve":
    !       set to limit sum of weight of all solutions to 1
    !    bit 1:
    !       use tanh solver that limits each solution to 1
    !    bit 2:
    !       do not compute final chi**2, only have adjusted offsets (c)

    use utils, only: &
         signan

    use star_data, only: &
         set_star_data, abu_covariance, &
         init_ei2, &
         init_covaricance_const

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar, nel, nsol, ncov, flags

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

    real(kind=real64), dimension(nsol), intent(out) :: &
         f
    real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
         c

    real(kind=real64) :: &
         scale
    real(kind=real64), dimension(nel) :: &
         y

    integer(kind=int64) :: &
         i, k, ierr

    call set_star_data(obs, err, det, cov, nel, ncov, icdf)

    if ((ls > 3).or.(ls < -1)) then
       print*, '[fitness] ls=', ls
       error stop '[fitness] invalid local search option'
    endif
    if (ls == -1) then
       goto 1000
    endif

    if ((ls == 0).and.(nstar == 1)) then
       call init_ei2()
       call init_covaricance_const()
       if (icdf == 0) then
          do k = 1, nsol
             call analytic_solve(c(k,1), abu(k,1,:))
          enddo
       else
          do k = 1, nsol
             call single_solve(c(k,1), abu(k,1,:))
          enddo
       endif
    else if ((ls == 0).and.(nstar > 1)) then
       do k = 1, nsol
          scale = 0.1d0
          do i = 1, nel
             y(i) = sum(c(k,:) * abu(k,:,i))
          enddo
          call init_ei2()
          call init_covaricance_const()
          if (icdf == 0) then
             call analytic_solve(scale, y)
          else
             call single_solve(scale, y)
          endif
          c(k,:) = c(k,:) * min(scale, 1.d0 / sum(c(k,:)))
       enddo
    else if ((ls == 1).and.(icdf == 1)) then

       ! NR solver with modified x-axis

       call init_ei2()
       call init_covaricance_const()
       if (btest(flags, FLAGS_LIMITED_SOLVER_BIT)) then
          do k = 1, nsol
             call newton(c(k,:), abu(k,:,:), nstar, ierr)
             if (ierr == 1) then
                call psolve(c(k,:), abu(k,:,:), nstar)
                call newton(c(k,:), abu(k,:,:), nstar, ierr)
             endif
          enddo
       else
          do k = 1, nsol
             call newton2(c(k,:), abu(k,:,:), nstar, ierr)
             if (ierr == 1) then
                call psolve2(c(k,:), abu(k,:,:), nstar)
                call newton2(c(k,:), abu(k,:,:), nstar, ierr)
             endif
          enddo
       end if
    else if ((ls == 3).and.(icdf == 1)) then

       ! classical NR solver that converges poorly due to stiffness of log/exp

       call init_ei2()
       call init_covaricance_const()
       do k = 1, nsol
          call newton_classic(c(k,:), abu(k,:,:), nstar, ierr)
          if (ierr == 1) then
             call psolve(c(k,:), abu(k,:,:), nstar)
             call newton_classic(c(k,:), abu(k,:,:), nstar, ierr)
          endif
       enddo
    else

       ! Slower UOBYQA solver that does not depend on C(2)
       ! function for chi2

       if (btest(flags, FLAGS_LIMITED_SOLVER_BIT)) then
          do k = 1, nsol
             call psolve(c(k,:), abu(k,:,:), nstar)
          enddo
       else
          do k = 1, nsol
             call psolve2(c(k,:), abu(k,:,:), nstar)
          enddo
       endif
    endif

    if (btest(flags, FLAGS_LIMIT_SOLUTION_BIT)) then
       do k = 1, nsol
          scale = sum(c(k,:))
          if (scale > ALMOST_ONE) then
             scale = ALMOST_ONE / scale
             c(k,:) = c(k,:) * scale
          endif
       end do
    endif

    ! get final fit value for possibly adjusted scales

1000 continue

    if (btest(flags, FLAGS_NO_CHI2)) then
       goto 2000
    endif
    do k = 1, nsol
       call chi2(f(k), c(k,:), abu(k,:,:), nstar)
    enddo

2000 continue

  end subroutine fitness


  subroutine fitness_m(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)

    use utils, only: &
         signan

    use star_data, only: &
         set_star_data, abu_covariance, &
         init_ei2, &
         init_covaricance_const

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar, nel, nsol, ncov, flags

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

    real(kind=real64), dimension(nsol, nel, nel), intent(out) :: &
         f
    real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
         c

    integer(kind=int64) :: &
         k, xflags

    xflags = ibset(flags, FLAGS_NO_CHI2)

    call fitness(f(:,1,1), c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, xflags)

    do k = 1, nsol
       call chi2m(f(k,:,:), c(k,:), abu(k,:,:), nstar)
    enddo

  end subroutine fitness_m


  subroutine chi2(f, c, abu, nstar)

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         eri

    use star_data, only: &
         diff_covariance, &
         iupper, idetec, ierinv, iuncor, &
         nupper, ndetec

    use norm, only: &
         logcdf, logcdfp

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), intent(in), dimension(nstar)  :: &
         c
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(out) :: &
         f

    real(kind=real64), dimension(nel) :: &
         logy, diff_obs, diff_det
    integer(kind=int64) :: &
         i, i1

    do i = 1, nel
       logy(i) = log(sum(c(:) * abu(:,i))) * ln10i
    enddo

    diff_obs(:) = logy(:) - obs(:)
    diff_det(ierinv) = logy(ierinv) - det(ierinv)

    f = diff_covariance(diff_obs)
    f = f + sum((diff_obs(iuncor) * eri(iuncor))**2)
    if (icdf == 0) then
       do i1=1, nupper
          i = iupper(i1)
          if (diff_obs(i) > 0.d0) then
             f = f + (diff_obs(i) * eri(i))**2
          endif
       enddo
       do i1=1, ndetec
          i = idetec(i1)
          if (diff_det(i) < 0.d0) then
             f = f - (diff_det(i) * eri(i))**2
          endif
       enddo
    else
       f = f - 2.d0* sum(logcdf(diff_obs(iupper) * eri(iupper)))
       f = f + 2.d0* sum(logcdf(diff_det(idetec) * eri(idetec)))
    endif

  end subroutine chi2


  subroutine chi2m(f, c, abu, nstar)

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         eri

    use star_data, only: &
         diff_covariance_m, &
         iupper, idetec, ierinv, iuncor, &
         nupper, ndetec, icovar, nuncor

    use norm, only: &
         logcdf, logcdfp

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), intent(in), dimension(nstar)  :: &
         c
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(out), dimension(nel, nel) :: &
         f

    real(kind=real64), dimension(nel) :: &
         logy, diff_obs, diff_det
    integer(kind=int64) :: &
         i, i1

    do i = 1, nel
       logy(i) = log(sum(c(:) * abu(:,i))) * ln10i
    enddo

    diff_obs(:) = logy(:) - obs(:)
    diff_det(ierinv) = logy(ierinv) - det(ierinv)

    f(:,:) = 0.d0
    f(icovar,icovar) = diff_covariance_m(diff_obs)

    do i1=1,nuncor
       i = iuncor(i1)
       f(i,i) = f(i,i) + (diff_obs(i) * eri(i))**2
    enddo
    if (icdf == 0) then
       do i1=1, nupper
          i = iupper(i1)
          if (diff_obs(i) > 0.d0) then
             f(i,i) = f(i,i) + (diff_obs(i) * eri(i))**2
          endif
       enddo
       do i1=1, ndetec
          i = idetec(i1)
          if (diff_det(i) < 0.d0) then
             f(i,i) = f(i,i) - (diff_det(i) * eri(i))**2
          endif
       enddo
    else
       do i1=1, nupper
          i = iupper(i1)
          f(i,i) = f(i,i) - 2.d0*logcdf(diff_obs(i) * eri(i))
       enddo
       do i1=1, ndetec
          i = idetec(i1)
          f(i,i) = f(i,i) + 2.d0*logcdf(diff_det(i) * eri(i))
       enddo
    endif

  end subroutine chi2m


  subroutine chi2_prime(f1, f2, c)

    ! return first and second derivative of chi2

    use utils, only: &
         signan

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         eri, ei2, &
         icovar, ncovar, &
         nuncor, iuncor, &
         ndetec, idetec, &
         nupper, iupper, &
         ierinv, iernoi, &
         diff_zv, &
         reduced_diff_zv

    use abu_data, only: &
         abu, nstar

    use norm, only: &
         logcdf, logcdfp

    implicit none

    real(kind=real64), parameter :: &
         WALL_LOC = 1.d-15, &
         WALL_POWER = 1.0d0

    real(kind=real64), intent(in), dimension(nstar)  :: &
         c

    real(kind=real64), intent(out), dimension(nstar) :: &
         f1
    real(kind=real64), intent(out), dimension(nstar, nstar) :: &
         f2

    real(kind=real64), dimension(nel) :: &
         y, logy, yi, yi2, &
         diff_obs, diff_obs_eri, diff_obs_ei2, &
         diff_det, diff_det_eri, diff_det_ei2
    real(kind=real64), dimension(ncovar) :: &
         zv, yp, yc, ypp, yic, yc2, yp2, zvp
    real(kind=real64) :: &
         fi1, fi2, fa, fc
    integer(kind=int64) :: &
         i, i1, j, k

    do i = 1, nel
       y(i) = sum(c(:) * abu(:,i))
    enddo
    logy(:) = log(y(:)) * ln10i

    yi(:) = 1.d0 / y(:)
    yi2(:) = yi(:)**2

    diff_obs(:) = logy(:) - obs(:)
    diff_obs_ei2(ierinv) = diff_obs(ierinv) * ei2(ierinv)
    diff_obs_ei2(iernoi) = signan()

    diff_det(:) = logy(:) - det(:)
    diff_det_ei2(ierinv) = diff_det(ierinv) * ei2(ierinv)
    diff_det_ei2(iernoi) = signan()

    f1(:)   = 0.d0
    f2(:,:) = 0.d0

    ! only compute lower triangle of the Jacobian matrix as it is symmetric

    ! covariant errors

    zv(:) = diff_zv(diff_obs)

    do j = 1, nstar
       yc(:) = abu(j, icovar)
       yic(:) = yi(icovar)
       yp(:) = yc(:) * yic(:)
       f1(j) = f1(j) + sum(zv(:) * yp(:))
       zvp(:) = reduced_diff_zv(yp)
       do k = 1, j
          yc2(:) = abu(k, icovar)
          yp2(:) = yc2(:) * yic(:)
          ypp(:) = yp(:) * yic(:) * yc2(:)
          f2(j,k) = f2(j,k) - sum(zv(:) * ypp(:)) + sum(zvp(:) * yp2(:)) * ln10i
       enddo
    enddo

    ! uncorrelated errors

    do i1 = 1, nuncor
       i = iuncor(i1)
       fi1 = diff_obs_ei2(i) * yi(i)
       fi2 = (ei2(i) * ln10i - diff_obs_ei2(i)) * yi2(i)
       do j = 1, nstar
          f1(j) = f1(j) + fi1 * abu(j,i)
          fa = fi2 * abu(j,i)
          do k=1,j
             f2(j,k) = f2(j,k) + fa * abu(k,i)
          enddo
       enddo
    enddo

    if (icdf == 0) then

       ! upper limits

       do i1 = 1, nupper
          i = iupper(i1)
          if (diff_obs(i) <= 0.d0) &
               cycle
          fi1 = diff_obs_ei2(i) * yi(i)
          fi2 = (ei2(i) * ln10i - diff_obs_ei2(i)) * yi2(i)
          do j = 1, nstar
             f1(j) = f1(j) + fi1 * abu(j,i)
             fa = fi2 * abu(j,i)
             do k=1,j
                f2(j,k) = f2(j,k) + fa * abu(k,i)
             enddo
          enddo
       enddo

       ! detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          if (diff_det(i) >= 0.d0) &
               cycle
          fi1 = diff_det_ei2(i) * yi(i)
          fi2 = (ei2(i) * ln10i - diff_det_ei2(i)) * yi2(i)
          do j = 1, nstar
             f1(j) = f1(j) - fi1 * abu(j,i)
             fa = fi2 * abu(j,i)
             do k=1,j
                f2(j,k) = f2(j,k) - fa * abu(k,i)
             enddo
          enddo
       enddo

    else

       diff_obs_eri(ierinv) = diff_obs(ierinv) * eri(ierinv)
       diff_det_eri(ierinv) = diff_det(ierinv) * eri(ierinv)
       diff_obs_eri(iernoi) = signan()
       diff_det_eri(iernoi) = signan()

       ! upper limits

       do i1 = 1, nupper
          i = iupper(i1)
          fc = logcdfp(diff_obs(i) * eri(i))
          fi1 = eri(i) * yi(i) * fc
          fi2 = yi(i) * (1.d0 + ln10i * eri(i) * (fc + diff_obs(i) * eri(i)))
          do j=1, nstar
             fa = fi1 * abu(j, i)
             f1(j) = f1(j) - fa
             do k=1, j
                f2(j,k) = f2(j,k) + fi2 * fa * abu(k,i)
             enddo
          enddo
       enddo

       ! detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          fc = logcdfp(diff_det(i) * eri(i))
          fi1 = eri(i) * yi(i) * fc
          fi2 = yi(i) * (1.d0 + ln10i * eri(i) * (fc + diff_det(i) * eri(i)))
          do j=1, nstar
             fa = fi1 * abu(j, i)
             f1(j) = f1(j) + fa
             do k=1, j
                f2(j,k) = f2(j,k) - fi2 * fa * abu(k,i)
             enddo
          enddo
       enddo

    endif

    ! these identical factors make no difference in solver
    ! Just for correctness while debugging

    do j=1, nstar
       f1(j) = f1(j) * ln10i2
       do k=1, j
          f2(j,k) = f2(j,k) * ln10i2
       enddo
    enddo

    ! fill in symmetric rest of jacobian

    do j = 1, nstar-1
       do k = j+1, nstar
          f2(j,k) = f2(k,j)
       enddo
    enddo

    ! penalty for low c - for newton_classic

    if (wall_chi2_prime) then

       do j = 1, nstar
          fa =  1.d0 / c(j)
          fc = (WALL_LOC * fa) ** WALL_POWER
          fi1 = -fc * fa * WALL_POWER
          fi2 = -fi1 * fa * (WALL_POWER + 1.d0)
          f1(j) = f1(j) + fi1
          f2(j,j) = f2(j,j) + fi2
       enddo
    end if

  end subroutine chi2_prime


  subroutine newton_classic(c, abu, nstar, ierr)

    use mleqs, only: &
         leqs

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data

    implicit none

    real(kind=real64), parameter :: &
         FMIN = 0.5d0

    integer(kind=int64), parameter :: &
         max_steps = 1000

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu
    integer(kind=int64), intent(out) :: &
         ierr

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         cx, dc, f1
    real(kind=real64), dimension(nstar, nstar) :: &
         f2
    real(kind=real64) :: &
         dcr

    integer(kind=int64) :: &
         i

    ! Initial N-R values

    wall_chi2_prime = .true.

    cx(:) = max(min(c(:), 1.d0), 1e-12)

    call set_abu_data(abu, nel, nstar)

    do i = 1, max_steps
       call chi2_prime(f1, f2, cx)
       dc(:) = leqs(f2, f1, nstar)
       dcr = maxval(abs(dc) / c)
       if (dcr > FMIN) then
          dc(:) = dc(:) * (FMIN / dcr)
       endif
       cx(:) = cx(:) - dc(:)
       if (dcr < 1.0d-6) goto 1000
    enddo

    if (stop_on_nonconvergence) then
       error stop '[newton_classic] did not converge.'
    endif
    ierr = 1
    return

1000 continue

    ierr = 0
    c(:) = cx(:)
    return

  end subroutine newton_classic


  subroutine newton_prime(f1, f2, x)

    use abu_data, only: &
         nstar

    implicit none

    real(kind=real64), parameter :: &
         WALL_LOC = 12.d0

    real(kind=real64), intent(in), dimension(nstar) :: &
         x

    real(kind=real64), intent(out), dimension(nstar) :: &
         f1
    real(kind=real64), intent(out), dimension(nstar, nstar) :: &
         f2

    real(kind=real64), dimension(nstar) :: &
         g1
    real(kind=real64), dimension(nstar, nstar) :: &
         g2

    real(kind=real64), dimension(nstar) :: &
         t, c, cp, cpp
    real(kind=real64) :: &
         p, e
    integer(kind=int64) :: &
         j, k
    logical :: &
         wall

    ! Transform from -inf/inf to 0/1
    ! x is used for the solver, tanhx is physically meaningful

    t(:) = tanh(x)
    c(:) = 0.5d0 * (1.0d0 + t(:))

    ! Calculate the chi2

    call chi2_prime(g1, g2, c)

    cp(:) = 0.5d0 * (1.d0 - t(:)**2)
    cpp(:) = -2.d0 * t(:) * cp(:)

    do j = 1, nstar
       do k = 1, nstar
          f2(j,k) = g2(j,k) * cp(j) * cp(k)
       enddo
       f1(j) = g1(j) * cp(j)
       f2(j,j) = f2(j,j) + g1(j) * cpp(j)
    enddo

    ! Build a wall at zero

    do j = 1, nstar
       if (x(j) < -WALL_LOC) then
          p = x(j) + WALL_LOC
          wall = .true.
       else if (x(j) > WALL_LOC) then
          p = x(j) - WALL_LOC
          wall = .true.
       else
          wall=.false.
       endif
       if (wall) then
          e = exp(p**2)
          f1(j) = f1(j) + 2.d0 * p * e
          f2(j,j) = f2(j,j) + (4.d0 * p**2 + 2.d0) * e
       endif
    enddo

  end subroutine newton_prime


  subroutine newton(c, abu, nstar, ierr)

    ! based on apprach used for psolve

    use utils, only: &
         signan

    use mleqs, only: &
         leqs

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data

    implicit none

    real(kind=real64), parameter :: &
         FMIN = 0.5d0

    integer(kind=int64), parameter :: &
         max_steps = 100

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu
    integer(kind=int64), intent(out) :: &
         ierr

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         x, dx, f1
    real(kind=real64), dimension(nstar, nstar) :: &
         f2
    real(kind=real64) :: &
         dxr, c1, c2

    integer(kind=int64) :: &
         iter

    ! Initial N-R values

    wall_chi2_prime = .false.

    if (any(c >= 1.d0)) then
       if (stop_on_large_offset) then
          print*, '[newton] DEBUG IN: c = ', c
          error stop '[newton] c >= 1'
       endif
       c = signan()
       return
    endif

    call set_abu_data(abu, nel, nstar)

    ! Convert offsets to solver space

    x = atanh(max(min(c(:), ALMOST_ONE), 1e-12) * 2.d0 - 1.d0)

    do iter = 1, max_steps
       call newton_prime(f1, f2, x)
       dx(:) = leqs(f2, f1, nstar)
       dxr = maxval(abs(dx))
       if (dxr > FMIN) then
          dx(:) = dx(:) * (FMIN / dxr)
       endif

       x(:) = x(:) - dx(:)
       if (dxr < 1.0d-6) goto 1000
    enddo

    call chi2(c1, 0.5d0 * (1.d0 + tanh(x(:))), abu, nstar)
    call chi2(c2, 0.5d0 * (1.d0 + tanh(x(:) + dx(:))), abu, nstar)
    if (abs(c1 - c2) / (abs(c1 + c2) + 1.d-99) < 1e-8) then
       goto 1000
    endif

    if (stop_on_nonconvergence) then
       error stop '[newton] did not converge.'
    end if
    ierr = 1
    return

1000 continue

    ! Convert solver space to offsets

    ierr = 0
    c(:) = 0.5d0 * (1.d0 + tanh(x(:)))

  end subroutine newton


  subroutine newton2_prime(f1, f2, x)

    use abu_data, only: &
         nstar

    implicit none

    real(kind=real64), parameter :: &
         WALL_LOC = -12.d0 * ln10

    real(kind=real64), intent(in), dimension(nstar) :: &
         x

    real(kind=real64), intent(out), dimension(nstar) :: &
         f1
    real(kind=real64), intent(out), dimension(nstar, nstar) :: &
         f2

    real(kind=real64), dimension(nstar) :: &
         g1
    real(kind=real64), dimension(nstar, nstar) :: &
         g2

    real(kind=real64), dimension(nstar) :: &
         c
    real(kind=real64) :: &
         p, e
    integer(kind=int64) :: &
         j, k

    ! Transform from -inf/inf to 0/1
    ! x is used for the solver, tanhx is physically meaningful

    c(:) = exp(x(:))

    ! Calculate the chi2

    call chi2_prime(g1, g2, c)

    do j = 1, nstar
       do k = 1, nstar
          f2(j,k) = g2(j,k) * c(j) * c(k)
       enddo
       f1(j) = g1(j) * c(j)
       f2(j,j) = f2(j,j) + f1(j)
    enddo

    ! Build a wall at zero

    do j = 1, nstar
       if (x(j) < WALL_LOC) then
          p = (x(j) - WALL_LOC) * ln10i
          e = exp(p**2)
          f1(j) = f1(j) + ln10i2 * p * e
          f2(j,j) = f2(j,j) + (ln10i2p2 * p**2 + ln10ip22) * e
       endif
    enddo

  end subroutine newton2_prime


  subroutine newton2(c, abu, nstar, ierr)

    ! based on apprach used for psolve

    use utils, only: &
         signan

    use mleqs, only: &
         leqs

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data

    implicit none

    real(kind=real64), parameter :: &
         FMIN = 0.5d0

    integer(kind=int64), parameter :: &
         max_steps = 200

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), dimension(nstar, nel), intent(in) :: &
         abu

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c
    integer(kind=int64), intent(out) :: &
         ierr

    real(kind=real64), dimension(nstar) :: &
         x, dx, f1
    real(kind=real64), dimension(nstar, nstar) :: &
         f2
    real(kind=real64) :: &
         dxr, c1, c2

    integer(kind=int64) :: &
         iter

    ! Initial N-R values

    wall_chi2_prime = .false.

    if (any(c <= 0.d0)) then
       if (stop_on_zero_offset) then
          print*, '[newton2] DEBUG IN: c = ', c
          error stop '[newton2] c < 0'
       endif
       c = signan()
       return
    endif

    call set_abu_data(abu, nel, nstar)

    ! Convert offsets to solver space

    x = log(max(c(:), 1e-12))

    do iter = 1, max_steps
       call newton2_prime(f1, f2, x)
       dx(:) = leqs(f2, f1, nstar)
       dxr = maxval(abs(dx))
       if (dxr > FMIN) then
          dx(:) = dx(:) * (FMIN / dxr)
       endif

       x(:) = x(:) - dx(:)
       if (dxr < 1.0d-12) goto 1000
    enddo

    call chi2(c1, exp(x(:)), abu, nstar)
    call chi2(c2, exp(x(:) + dx(:)), abu, nstar)
    if (abs(c1 - c2) / (abs(c1 + c2) + 1.d-99) < 1e-8) then
       goto 1000
    endif

    if (stop_on_nonconvergence) then
       error stop '[newton2] did not converge.'
    end if
    ierr = 1
    return

1000 continue

    ! Convert solver space to offsets

    ierr = 0
    c = exp(x)

  end subroutine newton2


  subroutine single_prime(x, f1, f2)

    ! Returns the first and second derivatives of (1/2)*chi^2 with
    ! respect to x for single star fits

    ! should not actually be used for icdf == 0

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         idetec, iuncor, iupper, idetec, &
         nupper, ndetec, &
         eri, ei2, &
         mp, &
         diff_z

    use abu_data, only: &
         logabu

    use norm, only: &
         logcdfp

    implicit none

    real(kind=real64), intent(in) :: &
         x

    real(kind=real64), intent(out) :: &
         f1, f2

    real(kind=real64), dimension(nel) :: &
         diff_obs, diff_det
    real(kind=real64) :: &
         de, df
    integer(kind=int64) :: &
         i, i1

    diff_obs(:) = logabu(:) - obs(:) + x
    diff_det(idetec) = logabu(idetec) - det(idetec) + x

    ! correlated errors

    f1 = diff_z(diff_obs)
    f2 = mp

    ! Uncorrelated errors

    f1 = f1 + sum(diff_obs(iuncor) * ei2(iuncor))
    f2 = f2 + sum(ei2(iuncor))

    if (icdf == 0) then

       ! This approach is not really good for NR as chi^2 is only C(1)
       ! but we really need at least C(2)

       ! Upper limits if error < 0 (NOTE: changes sign)

       do i1 = 1, nupper
          i = iupper(i1)
          if (diff_obs(i) > 0.d0) then
             f1 = f1 + diff_obs(i) * ei2(i)
             f2 = f2 + ei2(i)
          endif
       enddo

       ! Detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          if (diff_det(i) < 0.d0) then
             f1 = f1 - diff_det(i) * ei2(i)
             f2 = f2 - ei2(i)
          endif
       enddo
    else

       ! Upper limits if error < 0 (NOTE: changes sign)

       do i1 = 1, nupper
          i = iupper(i1)
          df = - logcdfp(diff_obs(i)*eri(i)) * eri(i)
          de = - diff_obs(i) * ei2(i)
          f1 = f1 + df
          f2 = f2 + df * (df + de)
       enddo

       ! Detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          df = + logcdfp(diff_det(i)*eri(i)) * eri(i)
          de = - diff_det(i) * ei2(i)
          f1 = f1 + df
          f2 = f2 + df * (df + de)
       enddo
    endif

    ! these identical factors make no difference in solver
    ! Just for correctness while debugging

    f1 = f1 * 2.d0
    f2 = f2 * 2.d0

  end subroutine single_prime


  subroutine single_solve(c, abu)

    ! single NR solver

    ! should not be used for icdf == 0

    use utils, only: &
         signan

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data, &
         init_logabu

    implicit none

    integer(kind=int64), parameter :: &
         one = 1.d0

    integer(kind=int64), parameter :: &
         max_steps = 20

    real(kind=real64), intent(in), dimension(1, nel) :: &
         abu

    real(kind=real64), intent(inout) :: &
         c

    real(kind=real64) :: &
         x, f1, f2, delta

    integer(kind=int64) :: &
         iter

    call set_abu_data(abu, nel, one)
    call init_logabu()

    x = log(c) * ln10i
    do iter = 1, max_steps
       call single_prime(x, f1, f2)
       if (f2 == 0.d0) exit
       delta = f1 / f2
       x = x - delta
       if (abs(delta) < 1.0d-6) goto 1000
    enddo

    if (stop_on_nonconvergence) then
       error stop '[single_solve] did not converge.'
    end if
    c = signan()
    return

1000 continue
    c = exp(x * ln10)

  end subroutine single_solve


  subroutine analytic_solve(c, abu)

    use utils, only: &
         signan

    use star_data, only: &
         nel, &
         obs, det, upper, &
         ei2, mp, &
         inocov, nnocov, idetec, ndetec, &
         diff_z

    implicit none

    real(kind=real64), intent(in), dimension(nel) :: &
         abu
    real(kind=real64), intent(out) :: &
         c

    real(kind=real64) :: &
         x, diff, ei2s, z
    real(kind=real64), dimension(nel)  :: &
         diff_obs, diff_det, logabu

    logical, dimension(nel) :: &
         include_obs, include_det
    logical :: &
         change

    integer(kind=int64) :: &
         i, i1

    logabu(:) = log(abu(:)) * ln10i

    diff_obs(:) = logabu(:) - obs(:)
    diff_det(:) = logabu(:) - det(:)

    z = diff_z(diff_obs)

    ! First, include all upper limits as if they were detections, and
    ! exclude all detection thresholds as it you were above

    change = .true.

    include_obs(:) = .true.
    include_det(:) = .false.

    do while (change)

       ! (re)calculate with ignored elements and included thresholds

       ei2s = mp
       diff = z
       do i1 = 1, nnocov
          i = inocov(i1)
          if (include_obs(i)) then
             ei2s = ei2s + ei2(i)
             diff = diff + diff_obs(i) * ei2(i)
          endif
       enddo
       do i1 = 1, ndetec
          i = idetec(i1)
          if (include_det(i)) then
             ei2s = ei2s - ei2(i)
             diff = diff - diff_det(i) * ei2(i)
          endif
       enddo

       if (ei2s == 0.d0) then
          if (stop_on_nonconvergence) then
             error stop '[analytic] all values below threshold'
          endif
          c = signan()
          return
       endif

       x = -diff / ei2s

       change = .false.
       do i1 = 1, nnocov
          i = inocov(i1)
          if (upper(i).and.(include_obs(i))) then
             diff = diff_obs(i) + x
             if (diff < 0.0d0) then
                include_obs(i) = .false.
                change = .true.
             endif
          endif
       enddo
       do i1 = 1, ndetec
          i = idetec(i1)
          if (.not.include_det(i)) then
             diff = diff_det(i) + x
             if (diff < 0.0d0) then
                include_det(i) = .true.
                change = .true.
             endif
          endif
       enddo
    enddo

    c = exp(x * ln10)

  end subroutine analytic_solve


  subroutine psolve_chi2(nstar, x, f)

    use abu_data, only: &
         abu

    implicit none

    real(kind=real64), parameter :: &
         WALL_LOC = 12.d0

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), intent(in), dimension(nstar) :: &
         x

    real(kind=real64), intent(out) :: &
         f

    real(kind=real64), dimension(nstar) :: &
         c
    integer(kind=int64) :: &
         j

    ! Transform from -inf/inf to 0/1
    ! x is used for the solver, tanhx is physically meaningful

    c(:) = 0.5d0 * (1.0d0 + tanh(x(:)))

    ! Calculate the chi2

    call chi2(f, c(:), abu(:,:), nstar)

    ! Build a wall at zero

    do j = 1, nstar
       if (abs(x(j)) > WALL_LOC) then
          f = f + (exp((abs(x(j)) - WALL_LOC)**2) - 1.d0)
       endif
    enddo

  end subroutine psolve_chi2


  subroutine psolve(c, abu, nstar)

    use utils, only: &
         signan

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data

    use powell, only: &
         uobyqa

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), intent(in), dimension(nstar, nel) :: &
         abu

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         x
    real(kind=real64) :: &
         rhobeg, rhoend
    integer(kind=int64) :: &
         calls, iprint

    ! save module data for calfun

    ! an pointer to an allocated data structure should be used instead.

    call set_abu_data(abu, nel, nstar)

    ! Options for the uobyqa solver

    rhobeg = 1.d0
    rhoend = 1.d-6
    iprint = 0

    calls = nstar**4 + 8 * nstar**3 + 23 * nstar**2 + 42 * nstar + &
         max(2 * nstar**2 + 4, 18 * nstar)

    if (any(c >= 1.d0)) then
       if (stop_on_large_offset) then
          print*, '[psolve] DEBUG IN: c = ', c
          error stop '[psolve] c >= 1'
       endif
       c = signan()
       return
    endif

    ! Convert offsets to solver space

    x = atanh(c * 2.d0 - 1.d0)

    ! Call solver

    call uobyqa(psolve_chi2, nstar, x, rhobeg, rhoend, iprint, calls)

    ! Convert solver space to offsets

    c = 0.5d0 * (1.d0 + tanh(x))

  end subroutine psolve


  subroutine psolve2_chi2(nstar, x, f)

    use abu_data, only: &
         abu

    implicit none

    real(kind=real64), parameter :: &
         WALL_LOC = -12.d0 * ln10

    integer(kind=int64), intent(in) :: &
         nstar
    real(kind=real64), intent(in), dimension(nstar) :: &
         x

    real(kind=real64), intent(out) :: &
         f

    real(kind=real64), dimension(nstar) :: &
         c
    integer(kind=int64) :: &
         j

    ! Transform log to linear
    ! x is used for the solver, log(x) is physically meaningful

    c = exp(x)

    ! Calculate the chi2

    call chi2(f, c(:), abu(:,:), nstar)

    ! Build a wall at zero

    do j = 1, nstar
       if (x(j) < WALL_LOC) then
          f = f + (exp((x(j) * ln10i - WALL_LOC)**2) - 1.d0)
       endif
    enddo

  end subroutine psolve2_chi2


  subroutine psolve2(c, abu, nstar)

    use utils, only: &
         signan

    use star_data, only: &
         nel

    use abu_data, only: &
         set_abu_data

    use powell, only: &
         uobyqa

    implicit none

    integer(kind=int64), intent(in) :: &
         nstar
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

    if (any(c <= 0.d0)) then
       if (stop_on_zero_offset) then
          print*, '[psolve2] DEBUG IN: c = ', c
          error stop '[psolve2] c <= 0'
       endif
       c = signan()
       return
    endif

    ! Convert offsets to solver space

    x = log(c)

    ! Call solver

    call uobyqa(psolve2_chi2, nstar, x, rhobeg, rhoend, iprint, calls)

    ! Convert solver space to offsets

    c = exp(x)

  end subroutine psolve2

end module fitting
