module fitting

  use type_def, only: &
       real64, int64

  implicit none

  private

  real(kind=real64), parameter :: &
       ln10 = log(10.0d0), &
       ln10i = 1.d0 / ln10, &
       ln10i2 = 2.d0 / ln10, &
       ln10i2m = -ln10i2

  logical :: &
       wall_chi2_prime = .false.

  public :: &
       fitness, chi2

contains

  ! Due to poor f2py / numpy.distutils code we need to carry nel
  ! rather than taking it from star_data

  ! TODO - potentiall use pointer, e.g., for star object, rather than
  ! soreing data in a module

  subroutine fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf)

    use star_data, only: &
         set_star_data, abu_covariance, &
         init_erri2, &
         init_covaricance_const

    use type_def, only: &
         real64, int64

    implicit none

    integer(kind=int64), parameter :: &
         lsolve = 1 ! set to 1 to limit weight of all solutions to 1

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
         i, k

    real(kind=real64), dimension(nsol), intent(out) :: &
         f
    real(kind=real64), dimension(nsol, nstar), intent(inout) :: &
         c

    ! If localsearch is enabled, modify the offsets first
    ! otherwise keep relative weights fixed

    call set_star_data(obs, err, det, cov, nel, ncov, icdf)

    if ((ls > 3).or.(ls < -1)) then
       error stop 'invalid local search option'
    endif
    if (ls == -1) then
       goto 1000
    endif

    ! if (nstar == 1) then
    if ((ls == 0).and.(nstar == 1)) then
       call init_erri2()
       call init_covaricance_const()
       if (icdf == 0) then

          !print*,'[fitness] analytic'

          do k = 1, nsol
             call analytic_solve(c(k,1), abu(k,1,:))
          enddo
       else

          !print*,'[fitness] single'

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
          call init_erri2()
          call init_covaricance_const()
          if (icdf == 0) then

             !print*,'[fitness] analytic - fixed'

             call analytic_solve(scale, y)
          else

             !print*,'[fitness] single - fixed'

             call single_solve(scale, y)
          endif

          !print*,'[fitness] scale', scale

          c(k,:) = c(k,:) * min(scale, 1.d0 / sum(c(k,:)))
       enddo
    else if ((ls == 1).and.(icdf == 1)) then

       ! faster NR solver

       call init_erri2()
       call init_covaricance_const()

       !print*,'[fittness] newton'
       !print*,'[fittness] nstar',nstar
       !print*,'[fittness] c',c

       do k = 1, nsol
          call newton(c(k,:), abu(k,:,:), nstar)
       enddo

    else if ((ls == 3).and.(icdf == 1)) then

       ! faster NR solver

       call init_erri2()
       call init_covaricance_const()

       !print*,'[fittness] newton_classic'
       !print*,'[fittness] nstar',nstar
       !print*,'[fittness] c',c

       do k = 1, nsol
          call newton_classic(c(k,:), abu(k,:,:), nstar)
       enddo

    else

       ! Slower UOBYQA solver that does not depend on C(2)
       ! function for chi2

       !print*,'[fittness] psolve'

       !print*,'[fittness] nstar',nstar
       !print*,'[fittness] c',c

       do k = 1, nsol
          call psolve(c(k,:), abu(k,:,:), nstar)
       enddo

       !do k = 1, nsol
       !call chi2(f(k), c(k,:), abu(k,:,:), nstar)
       !print*, 'f=', f(k), 'c=', c(k,:)
       !enddo
       !call init_erri2()
       !call init_covaricance_const()

       !print*,'[fittness] newton'
       !print*,'[fittness] nstar',nstar
       !print*,'[fittness] c',c

       !do k = 1, nsol
       !call newton(c(k,:), abu(k,:,:), nstar)
       !enddo

    endif

    if (lsolve == 1) then
       do k = 1, nsol
          scale = sum(c(k,:))
          if (scale > 1.d0) then
             c(k,:) = c(k,:) * (1.d0 / scale)
          endif
       end do
    endif

    ! get final fit value for possibly adjusted scales

1000 continue

    do k = 1, nsol
       call chi2(f(k), c(k,:), abu(k,:,:), nstar)
    enddo

  end subroutine fitness


  subroutine chi2(f, c, abu, nstar)

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         erri

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

    !print*,'[chi2] A', f, nupper, ndetec

    f = f + sum((diff_obs(iuncor) * erri(iuncor))**2)

    !print*,'[chi2] B', f, iuncor
    !print*,'[chi2] B1 ', 'diff', diff_obs(iuncor)
    !print*,'[chi2] B2 ', 'erri', erri(iuncor)

    if (icdf == 0) then
       !print*,'[chi2] C', f, nupper, iupper
       do i1=1, nupper
          i = iupper(i1)
          if (diff_obs(i) > 0.d0) then
             f = f + (diff_obs(i) * erri(i))**2
          endif
       enddo
       ! f = f + sum((min(diff(iupper), 0.d0) * erri(iupper))**2)
       !print*,'[chi2] D', f, ndetec, idetec
       do i1=1, ndetec
          i = idetec(i1)
          if (diff_det(i) < 0.d0) then
             f = f - (diff_det(i) * erri(i))**2
          endif
       enddo
       ! f = f - sum((max(ddff(imeasu), 0.d0) * erri(imeasu))**2)
    else
       !print*,'[chi2] E', f, nupper, iupper
       f = f - 2.d0* sum(logcdf(diff_obs(iupper) * erri(iupper)))
       ! do i1=1, nupper
       !    i = iupper(i1)
       !    f = f - 2.d0*logcdf(diff_obs(j) * erri(j))
       ! enddo
       !print*,'[chi2] F', f, ndetec, idetec
       f = f + 2.d0* sum(logcdf(diff_det(idetec) * erri(idetec)))
       ! do i1=1, ndetec
       !    i = idetec(i1)
       !    f = f + 2.d0*logcdf(dif_det(i) * erri(i))
       ! enddo
    endif

    !print*,'[chi2] G', f

  end subroutine chi2


  subroutine chi2_prime(f1, f2, c)

    ! return first and second derivative of chi2

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         erri, erri2, &
         icovar, ncovar, &
         nuncor, iuncor, &
         ndetec, idetec, &
         nupper, iupper, &
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
         diff_obs, diff_obs_erri, diff_obs_erri2, &
         diff_det, diff_det_erri, diff_det_erri2
    real(kind=real64), dimension(ncovar) :: &
         zv, yp, yc, ypp, yic, yc2, yp2, zvp
    real(kind=real64) :: &
         fi1, fi2, fa, fc
    integer(kind=int64) :: &
         i, i1, j, k

    ! error stop '[chi2_prime] not implemented'

    do i = 1, nel
       y(i) = sum(c(:) * abu(:,i))
    enddo
    logy(:) = log(y(:)) * ln10i

    yi(:) = 1.d0 / y(:)
    yi2(:) = yi(:)**2

    diff_obs(:) = logy(:) - obs(:)
    diff_obs_erri2(:) = diff_obs(:) * erri2(:)

    diff_det(:) = logy(:) - det(:)
    diff_det_erri2(:) = diff_det(:) * erri2(:)

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

    !print*, '[cp] var c', c
    !print*, '[cp] covar'
    !print*, '[cp] f1(:)', f1
    !do j=1, nstar
    !print*, '[cp] f2(j,:)', f2(j,:)
    !enddo

    ! uncorrelated errors

    !print*, '[cp] nuncor, iuncor', nuncor, iuncor

    do i1 = 1, nuncor
       i = iuncor(i1)
       fi1 = diff_obs_erri2(i) * yi(i)
       fi2 = (erri2(i) * ln10i - diff_obs_erri2(i)) * yi2(i)
       do j = 1, nstar
          f1(j) = f1(j) + fi1 * abu(j,i)
          fa = fi2 * abu(j,i)
          do k=1,j
             f2(j,k) = f2(j,k) + fa * abu(k,i)
          enddo
       enddo
    enddo

    !print*, '[cp] uncor'
    !print*, '[cp] var f1(:)', f1
    !do j=1, nstar
    !print*, '[cp] var f2(j,:)', f2(j,:)
    !enddo

    if (icdf == 0) then

       ! upper limits

       do i1 = 1, nupper
          i = iupper(i1)
          if (diff_obs(i) <= 0.d0) &
               cycle
          fi1 = diff_obs_erri2(i) * yi(i)
          fi2 = (erri2(i) * ln10i - diff_obs_erri2(i)) * yi2(i)
          do j = 1, nstar
             f1(j) = f1(j) + fi1 * abu(j,i)
             fa = fi2 * abu(j,i)
             do k=1,j
                f2(j,k) = f2(j,k) + fa * abu(k,i)
             enddo
          enddo
       enddo

       ! detection thresholds

       !print*, 'erri2', erri2
       !print*, 'ndetec, idetec', ndetec, idetec

       do i1 = 1, ndetec
          i = idetec(i1)
          if (diff_det(i) >= 0.d0) &
               cycle
          fi1 = diff_det_erri2(i) * yi(i)
          fi2 = (erri2(i) * ln10i - diff_det_erri2(i)) * yi2(i)
          do j = 1, nstar
             f1(j) = f1(j) - fi1 * abu(j,i)
             fa = fi2 * abu(j,i)
             do k=1,j
                f2(j,k) = f2(j,k) - fa * abu(k,i)
             enddo
          enddo
       enddo

    else

       diff_obs_erri(:)  = diff_obs(:) * erri(:)
       diff_det_erri(:)  = diff_det(:) * erri(:)

       ! upper limits

       do i1 = 1, nupper
          i = iupper(i1)
          fc = logcdfp(diff_obs(i) * erri(i))
          fi1 = erri(i) * yi(i) * fc
          fi2 = yi(i) * (1.d0 + ln10i * erri(i) * (fc + diff_obs(i) * erri(i)))
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
          fc = logcdfp(diff_det(i) * erri(i))
          fi1 = erri(i) * yi(i) * fc
          fi2 = yi(i) * (1.d0 + ln10i * erri(i) * (fc + diff_det(i) * erri(i)))
          do j=1, nstar
             fa = fi1 * abu(j, i)
             f1(j) = f1(j) + fa
             do k=1, j
                f2(j,k) = f2(j,k) - fi2 * fa * abu(k,i)
             enddo
          enddo
       enddo

    endif

    !print*, '[cp] final'
    !print*, '[cp] var f1(:)', f1
    !do j=1, nstar
    !print*, '[cp] var f2(j,:)', f2(j,:)
    !enddo

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

          !print*,'[nr]', j, c(j), fi1, fi2

          f1(j) = f1(j) + fi1
          f2(j,j) = f2(j,j) + fi2
       enddo
    end if

  end subroutine chi2_prime


  subroutine newton_classic(c, abu, nstar)

    use mleqs, only: &
         leqs

    use star_data, only: &
         nel, &
         signan

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
    !integer(kind=int64) :: &
    !     j

    ! Initial N-R values

    wall_chi2_prime = .true.

    c(:) = max(min(c(:), 1.d0), 1e-12)

    call set_abu_data(abu, nel, nstar)

    do i = 1, max_steps
       call chi2_prime(f1, f2, c)
       dc(:) = leqs(f2, f1, nstar)

       !print*,''
       !print*,'[newton] i', i
       !print*,'[newton] c', c
       !print*,'[newton] f1', f1
       !do j=1, nstar
       !print*,'[newton] f2', f2(j, :)
       !enddo
       !print*,'[newton] dc',dc

       dcr = maxval(abs(dc) / c)
       if (dcr > FMIN) then
          dc(:) = dc(:) * (FMIN / dcr)
       endif

       !print*,'[newton] dcr',dcr
       !print*,'[newton] dc',dc

       c(:) = c(:) - dc(:)
       if (dcr < 1.0d-6) return
    enddo

    !print*, '[newton] did not converge after ',i-1,'iterations.'
    c(:) = signan()
    return
    error stop '[newton] did not converge.'

  end subroutine newton_classic


  subroutine newton_prime(f1, f2, x)

    use type_def, only: &
         int64, real64

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
         t, y, yp, ypp
    real(kind=real64) :: &
         p, e
    integer(kind=int64) :: &
         j, k
    logical :: &
         wall

    ! Transform from -inf/inf to 0/1
    ! x is used for the solver, tanhx is physically meaningful

    t(:) = tanh(x)
    y(:) = 0.5d0 * (1.0d0 + t(:))

    ! Calculate the chi2

    call chi2_prime(g1, g2, y)

    yp(:) = 0.5d0 * (1.d0 - t(:)**2)
    ypp(:) = -2.d0 * t(:) * yp(:)

    do j = 1, nstar
       do k = 1, nstar
          f2(j,k) = g2(j,k) * yp(j) * yp(k)
       enddo
       f1(j) = g1(j) * yp(j)
       f2(j,j) = f2(j,j) + g1(j) * ypp(j)
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
          !print*,'[newton_prime] p=', p
          e = exp(p**2)
          f1(j) = f1(j) + 2.d0 * p * e
          f2(j,j) = f2(j,j) + (4.d0 * p**2 + 2.d0) * e
       endif
    enddo

  end subroutine newton_prime


  subroutine newton(c, abu, nstar)

    ! based on apprach used for psolve

    use mleqs, only: &
         leqs

    use star_data, only: &
         nel, &
         signan

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

    real(kind=real64), intent(inout), dimension(nstar) :: &
         c

    real(kind=real64), dimension(nstar) :: &
         x, dx, f1
    real(kind=real64), dimension(nstar, nstar) :: &
         f2
    real(kind=real64) :: &
         dxr

    integer(kind=int64) :: &
         i
    !integer(kind=int64) :: &
    !j

    ! Initial N-R values

    wall_chi2_prime = .false.

    if (any(c >= 1.d0)) then
       print*, '[newton] DEBUG IN: c = ', c
       error stop 'c >= 1'
    endif

    c(:) = max(min(c(:), 1.d0), 1e-12)

    call set_abu_data(abu, nel, nstar)

    !Convert offsets to solver space
    x = atanh(c * 2.d0 - 1.d0)

    do i = 1, max_steps
       call newton_prime(f1, f2, x)
       dx(:) = leqs(f2, f1, nstar)

       !print*,''
       !print*,'[newton] i', i
       !print*,'[newton] x', x
       !print*,'[newton] f1', f1
       !do j=1, nstar
       !print*,'[newton] f2', f2(j, :)
       !enddo
       !print*,'[newton] dx',dx

       dxr = maxval(abs(dx))
       if (dxr > FMIN) then
          dx(:) = dx(:) * (FMIN / dxr)
          !print*,'[newton] dxr',dxr
          !print*,'[newton] dx',dx
       endif

       x(:) = x(:) - dx(:)
       if (dxr < 1.0d-6) goto 1000
    enddo

    !print*, '[newton] did not converge after ',i-1,'iterations.'
    c(:) = signan()
    return
    error stop '[newton] did not converge.'

1000 continue

    ! Convert solver space to offsets

    c = 0.5d0 * (1.d0 + tanh(x))

  end subroutine newton


  subroutine single_prime(x, f1, f2)

    ! Returns the first and second derivatives of (1/2)*chi^2 with
    ! respect to x for single star fits

    ! should not actually be used for icdf == 0

    use star_data, only: &
         nel, &
         icdf, obs, det, &
         idetec, iuncor, iupper, idetec, &
         nupper, ndetec, &
         erri, erri2, &
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

    !print*, '[sp] x', x
    !print*, '[sp] var f1', f1
    !print*, '[sp] var f2', f2

    ! Uncorrelated errors

    f1 = f1 + sum(diff_obs(iuncor) * erri2(iuncor))
    f2 = f2 + sum(erri2(iuncor))

    if (icdf == 0) then

       ! This approach is not really good for NR as chi^2 is only C(1)
       ! but we really need at least C(2)

       ! Upper limits if error < 0 (NOTE: changes sign)

       do i1 = 1, nupper
          i = iupper(i1)
          if (diff_obs(i) > 0.d0) then
             f1 = f1 + diff_obs(i) * erri2(i)
             f2 = f2 + erri2(i)
          endif
       enddo

       ! Detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          if (diff_det(i) < 0.d0) then
             f1 = f1 - diff_det(i) * erri2(i)
             f2 = f2 - erri2(i)
          endif
       enddo
    else

       ! Upper limits if error < 0 (NOTE: changes sign)

       do i1 = 1, nupper
          i = iupper(i1)
          df = - logcdfp(diff_obs(i)*erri(i)) * erri(i)
          de = - diff_obs(i) * erri2(i)
          f1 = f1 + df
          f2 = f2 + df * (df + de)
       enddo

       ! Detection thresholds

       do i1 = 1, ndetec
          i = idetec(i1)
          df = + logcdfp(diff_det(i)*erri(i)) * erri(i)
          de = - diff_det(i) * erri2(i)
          f1 = f1 + df
          f2 = f2 + df * (df + de)
       enddo
    endif

    ! these identical factors make no difference in solver
    ! Just for correctness while debugging

    f1 = f1 * 2.d0
    f2 = f2 * 2.d0

    !print*,'[sp] X', f1, f2

  end subroutine single_prime


  subroutine single_solve(c, abu)

    ! single NR solver

    ! should not be used for icdf == 0

    use star_data, only: &
         nel, &
         signan

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

    ! real(kind=real64), intent(in), dimension(nel), target :: &
    !      abu
    ! real(kind=real64), dimension(:,:), pointer :: &
    !      abu_

    real(kind=real64), intent(inout) :: &
         c

    real(kind=real64) :: &
         x, f1, f2, delta

    integer(kind=int64) :: &
         i

    ! abu_(1:1,1:nel) => abu
    ! call set_abu_data(abu_, nel, one)
    call set_abu_data(abu, nel, one)
    call init_logabu()

    x = log(c) * ln10i
    do i = 1, max_steps
       call single_prime(x, f1, f2)
       if (f2 == 0.d0) exit
       delta = f1 / f2
       x = x - delta
       if (abs(delta) < 1.0d-6) goto 1000
    enddo

    c = signan()
    return
    error stop '[singlesolve] did not converge.'

1000 continue
    c = exp(x * ln10)

  end subroutine single_solve


  subroutine analytic_solve(c, abu)

    use star_data, only: &
         nel, &
         obs, det, upper, &
         erri2, mp, &
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

    !print*, '[analytic] mp, z', mp, z

    ! First, include all upper limits as if they were detections, and
    ! exclude all detection thresholds as it you were above

    change = .true.

    include_obs(:) = .true.
    include_det(:) = .false.

    !print*,'nnocov, inocov', nnocov, inocov

    do while (change)

       ! (re)calculate with ignored elements and included thresholds

       ei2s = mp
       diff = z
       do i1 = 1, nnocov
          i = inocov(i1)
          if (include_obs(i)) then
             ei2s = ei2s + erri2(i)
             diff = diff + diff_obs(i) * erri2(i)
          endif
       enddo
       do i1 = 1, ndetec
          i = idetec(i1)
          if (include_det(i)) then
             ei2s = ei2s - erri2(i)
             diff = diff - diff_det(i) * erri2(i)
          endif
       enddo

       !print*,'[analytic] x', x

       if (ei2s == 0.d0) then
          error stop '[analytic] all values below threshold'
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


  subroutine psolve(c, abu, nstar)

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

    if (any(c >= 1.d0)) then
       print*, '[psolve] DEBUG IN: c = ', c
       error stop 'c >= 1'
    endif

    ! Convert offsets to solver space

    x = atanh(c * 2.d0 - 1.d0)

    ! Call solver
    call uobyqa(nstar, x, rhobeg, rhoend, iprint, calls)

    ! Convert solver space to offsets

    c = 0.5d0 * (1.d0 + tanh(x))

  end subroutine psolve

end module fitting

! =======================================================================

subroutine calfun(nstar, x, f)

  use type_def, only: &
       int64, real64

  use abu_data, only: &
       abu

  use fitting, only: &
       chi2

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
       tanhx
  integer(kind=int64) :: &
       i

  ! Transform from -inf/inf to 0/1
  ! x is used for the solver, tanhx is physically meaningful

  tanhx(:) = 0.5d0 * (1.0d0 + tanh(x(:)))

  ! Calculate the chi2

  call chi2(f, tanhx(:), abu(:,:), nstar)

  ! Build a wall at zero

  do i = 1, nstar
     if (abs(x(i)) > WALL_LOC) then
        f = f + (exp((abs(x(i)) - WALL_LOC)**2) - 1.d0)
     endif
  enddo

end subroutine calfun
