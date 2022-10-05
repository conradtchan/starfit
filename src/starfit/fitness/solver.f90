program test
    implicit none
    integer*8, parameter  :: dp = selected_real_kind(15)
    ! integer*8, parameter  :: nstar = 3, nel = 9
    integer*8, parameter  :: nstar = 1, nel = 1
    real(dp)            :: obs(nel)
    real(dp)            :: err(nel)
    real(dp)            :: abu(nstar, nel)
    real(dp)            :: c(nstar)
    real(dp)            :: f
    real(dp)            :: f1(nstar), f2(nstar)

    real(dp)            :: fup, fdo, f1up(nstar), f1do(nstar), f2up(nstar), f2do(nstar)
    real(dp)            :: cup(nstar), cdo(nstar)
    real(dp)            :: h
    integer*8             :: i

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

subroutine chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf)
    implicit none
    integer*8,        parameter           :: dp = selected_real_kind(15)

    real(dp),       parameter           :: ln10 = log(10.0d0)

    integer*8,        intent(in)          :: nstar, nel
    real(dp),       intent(in)          :: c(nstar)
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nstar, nel)
    integer*8,        intent(in)          :: icdf

    real(dp)        :: summed(nel)
    real(dp)        :: diff(nel)
    real(dp)        :: diffe(nel)
    real(dp)        :: errsum(nel)

    real(dp),       external            :: logcdf, logcdfp
    real(dp)                            :: fd1, fd2

    integer*8         :: i, j

    !f2py real(dp),       intent(out)       :: f
    !f2py real(dp),       intent(out)       :: f1(nstar), f2(nstar)
    real(dp),       intent(out)       :: f
    real(dp),       intent(out)       :: f1(nstar), f2(nstar)

    summed(:) = 0
    do i = 1, nstar
        do j = 1, nel
            summed(j) = summed(j) + c(i) * abu(i,j)
        enddo
    enddo

    diff(:) = obs(:) - log10(summed(:))
    diffe(:) = diff(:) / abs(err(:))
    errsum(:) = err(:) * summed(:)

    f = 0
    do i = 1, nel
        if (err(i) > 0) then
            f = f + diffe(i)**2
        else
            if (icdf == 1) then
                f = f - 2.d0*logcdf(diffe(i))
            else
                if (diff(i) < 0) then
                    f = f + diffe(i)**2
                endif
            endif
        endif
    enddo

    ! !Enable this before using the NR solver
    ! f1(:) = 0
    ! f2(:) = 0
    ! do i = 1, nstar
    !     do j = 1, nel
    !         if (err(j) > 0) then
    !             f1(i) = f1(i) + (diffe(j) * abu(i, j))/errsum(j)
    !             f2(i) = f2(i) + (abu(i, j) / errsum(j))**2 * (1/ln10 + diff(j))
    !         else
    !             if (icdf == 1) then
    !                 fd1 = logcdfp(diffe(j))
    !                 fd2 = -fd1 * (fd1 + diffe(j))
    !                 f1(i) = f1(i) - abu(i, j) * fd1/abs(errsum(j))
    !                 f2(i) = f2(i) - (1 / abs(err(j))) * (abu(i, j) / summed(j))**2 * (fd2/(ln10 * abs(err(j))) + fd1)
    !             else
    !                 if (diff(i) < 0) then
    !                     f1(i) = f1(i) + (diffe(j) * abu(i, j))/errsum(j)
    !                     f2(i) = f2(i) + (abu(i, j) / errsum(j))**2 * (1/ln10 + diff(j))
    !                 endif
    !             endif
    !         endif
    !     enddo
    ! enddo
    !
    ! f1(:) = f1(:) * (-2 / ln10)
    ! f2(:) = f2(:) * (2 / ln10)

end subroutine chisq

subroutine calfun(n, x, f)
    implicit none
    integer*8,      parameter           :: dp = selected_real_kind(15)
    integer*8,      parameter           :: nstarmax = 10, nelmax = 83

    integer*8,      intent(in)          :: n
    real(dp),       intent(in)          :: x(n)
    real(dp),       intent(out)         :: f

    real(dp)                            :: tanhx(n)

    integer*8                           :: nel_
    real(dp)                            :: obs_(nelmax)
    real(dp)                            :: err_(nelmax)
    real(dp)                            :: abu_(nstarmax, nelmax)
    integer*8                           :: icdf_

    real(dp)                            :: f1(n)
    real(dp)                            :: f2(n)

    integer*8                           :: i

    common /params/ nel_, abu_, obs_, err_, icdf_

    !Transform from -inf/inf to 0/1
    !x is used for the solver, tanhx is physically meaninful
    tanhx = 0.5d0 * (1.0d0 + tanh(x))

    !Calculate the chisq
    call chisq(f, f1, f2, tanhx, obs_, err_, abu_(1:n,1:nel_), n, nel_, icdf_)

    !Build a wall at zero
    do i = 1, n
        if (abs(x(i)) > 12.d0) then
            f = f* exp((abs(x(i)) - 12.d0)**2)
        endif
    enddo

end subroutine calfun

subroutine fitness(f, c, obs, err, abu, nstar, nel, nsol, ls, icdf)
    implicit none
    integer*8,        parameter         :: dp = selected_real_kind(15)

    integer*8,      intent(in)          :: nstar, nel, nsol
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nsol, nstar, nel)
    logical,        intent(in)          :: ls
    integer*8,      intent(in)          :: icdf

    real(dp)                            :: f1(nstar), f2(nstar), scale, summed(nel)

    integer*8                           :: i, j, k

    !f2py real(dp),       intent(out)       :: f(nsol)
    real(dp),       intent(out)       :: f(nsol)
    !f2py real(dp),       intent(in,out)          :: c(nsol, nstar)
    real(dp),       intent(inout)       :: c(nsol, nstar)

    !If localsearch is enabled, modify the offsets first
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
                call psolve(c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
            endif
            call chisq(f(i), f1, f2, c(i,:), obs, err, abu(i,:,:), nstar, nel, icdf)
        end do
    else
        do i = 1, nsol
            scale = 1.d0
            summed(:) = 0
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

subroutine newton(c, obs, err, abu, nstar, nel, icdf)
    implicit none
    integer*8,      parameter           :: dp = selected_real_kind(15)
    integer*8,      parameter           :: trials = 20
    integer*8,      intent(in)          :: nstar, nel

    real(dp),       intent(inout)       :: c(nstar)
    real(dp)                            :: cold(nstar)
    real(dp)                            :: starts(trials, nstar)
    real(dp)                            :: recordc(trials, nstar)
    real(dp)                            :: recordf(trials)
    real(dp)                            :: a(nstar)
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nstar, nel)
    integer*8,      intent(in)          :: icdf

    real(dp)                            :: f, fold
    real(dp)                            :: delta(nstar), deltaold(nstar)
    real(dp)                            :: f1(nstar), f2(nstar)
    integer*8                           :: i, j, jstop
    real(dp)                            :: maxstep, steplimit, deltalen
    real(dp)                            :: randoms(nstar)


    !Initial N-R values
    delta(:) = 0
    maxstep = 1.0d0
    fold = 1.d6
    f = 1.d3

    starts(1,:) = c(:)
    do j = 1, trials
        !Find random starting points
        if (j > 1) then
            call random_number(randoms)
            starts(j,:) = 10.0d0**(-5.0d0*randoms + 1.0d0)
        endif

        c(:) = starts(j,:)
        cold(:) = c(:)
        a(:) = log(c(:))
        do i = 1, 100
            fold = f
            !Evaluate function and derivatives
            call chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf)
            !If any NaNs are present, redo step with half the step size
            if (any(isnan(delta))) then
                delta(:) = 0.50d0*deltaold(:)
                deltaold(:) = delta(:)
                c(:) = cold(:)
                a(:) = log(c(:))
            !Calculate the step
            else
                deltaold(:) = delta(:)
                delta(:) = (f1(:) / (f1(:) + f2(:)*c(:)))
            endif
            !Calculate the length of the step
            deltalen = sqrt(dot_product(delta,delta))
            !Calculate how much the step needs to be scaled back
            steplimit = min(0.50d0, maxstep/deltalen)
            !Apply step
            ! print*,i, '|', f, '|', a, '|', c, '|', f1, '|', f2, '|', delta, '|', deltalen, maxstep, steplimit
            a(:) = a(:) - steplimit*delta(:)
            c(:) = exp(a(:))
            !Convergence detection
            if ((fold - f) < 1.d-6) then
                exit
            end if
        end do
        call chisq(f, f1, f2, c, obs, err, abu, nstar, nel, icdf)
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
    use Powell_Optimize
    implicit none

    integer*8,      parameter           :: dp = selected_real_kind(15)

    integer*8,      parameter           :: nstarmax = 10, nelmax = 83

    integer*8,      intent(in)          :: nstar, nel
    real(dp),       intent(inout)       :: c(nstar)
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nstar, nel)
    integer*8,      intent(in)          :: icdf

    integer*8                           :: nel_
    real(dp)                            :: abu_(nstarmax, nelmax)
    real(dp)                            :: obs_(nelmax)
    real(dp)                            :: err_(nelmax)
    integer*8                           :: icdf_

    real(dp)                            :: x(nstar)

    real(dp)                            :: rhobeg
    real(dp)                            :: rhoend
    integer*8                           :: calls
    integer*8                           :: iprint

    common /params/ nel_, abu_, obs_, err_, icdf_

    !Write to common block for calfun to get
    nel_ = nel
    abu_ = abu
    obs_ = obs
    err_ = err
    icdf_ = icdf

    !Options for the uobyqa solver
    rhobeg = 1.d0
    rhoend = 1.d-5
    iprint = 0
    calls = 50 * nstar

    !Convert offsets to solver space
    x = atanh(c * 2.d0 - 1.d0)

    !Call solver
    call uobyqa(nstar, x, rhobeg, rhoend, iprint, calls)
    !Convert solver space to offsets
    c = 0.5d0 * (1.d0 + tanh(x))

end subroutine psolve

subroutine singlesolve(c, obs, err, abu, nel, icdf)
    implicit none
    integer*8,      parameter           :: dp = selected_real_kind(15)

    integer*8,      intent(in)          :: nel
    real(dp),       intent(inout)       :: c
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nel)
    integer*8,      intent(in)          :: icdf

    real(dp)                            :: x, f, f1, f2, delta

    integer*8                           :: i

    x = log10(c)
    do i = 1, 20
        call prime(x, f1, f2, obs, err, abu, nel, icdf)
        delta = f1 / f2
        x = x - delta
        if (abs(delta) < 1.d-6) then
            c = 10**x
            exit
        endif
    enddo
end subroutine singlesolve

subroutine analyticsolve(c, obs, err, abu, nel)
    implicit none
    integer*8,      parameter           :: dp = selected_real_kind(15)

    integer*8,      intent(in)          :: nel
    real(dp),       intent(out)         :: c
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: err(nel)
    real(dp),       intent(in)          :: abu(nel)

    real(dp)                            :: x, newdiff
    real(dp)                            :: ei2(nel), diff(nel), logabu(nel)

    logical                             :: ignore(nel)
    logical                             :: change

    integer*8                           :: i

    logabu(:) = log10(abu(:))

    !First, ignore upper limits
    do i = 1, nel
        ei2(i) = 1 / err(i)**2
        diff(i) = (logabu(i) - obs(i)) * ei2(i)
    enddo
    x = -sum(diff) / sum(ei2)

    !Check position
    change = .false.
    ignore(:) = .false.
    do i = 1, nel
        if (err(i) < 0) then
            newdiff = logabu(i) + x - obs(i)
            if (newdiff < 0) then
                ignore(i) = .true.
                change = .true.
            endif
        endif
    enddo

    !Make adjustments until okay
    do while (change)
        !Recalculate with ignored elements
        do i = 1, nel
            if (ignore(i)) then
                ei2(i) = 0
                diff(i) = 0
            else
                ei2(i) = 1 / err(i)**2
                diff(i) = (logabu(i) - obs(i)) * ei2(i)
            endif
        enddo
        x = -sum(diff) / sum(ei2)

        !Check position
        change = .false.
        do i = 1, nel
            if ((err(i) < 0) .and. .not. ignore(i)) then
                newdiff = logabu(i) + x - obs(i)
                if (newdiff < 0) then
                    ignore(i) = .true.
                    change = .true.
                endif
            endif
        enddo
    enddo

    c = 10**x

end subroutine analyticsolve

subroutine prime(x, f1, f2, obs, error, abu, nel, icdf)
    implicit none
    !Returns the derivatives with respect to log(s)
    !For single star fits only!
    integer*8,      parameter           :: dp = selected_real_kind(15)

    real(dp)                            :: logcdfp

    integer*8,      intent(in)          :: nel
    real(dp),       intent(in)          :: abu(nel)
    real(dp),       intent(in)          :: obs(nel)
    real(dp),       intent(in)          :: error(nel)
    real(dp),       intent(in)          :: x
    integer*8,      intent(in)          :: icdf

    real(dp),       intent(out)         :: f1, f2

    real(dp)                            :: ei(nel), ei2(nel)
    real(dp)                            :: diff, de, df, e2

    integer*8                           :: i

    f1 = 0.d0
    f2 = 0.d0

    if (icdf == 0) then
        ei2 = 1.d0 / error**2
        do i = 1, nel
            diff = log10(abu(i)) + x - obs(i)
            !Upper limits
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
            !Upper limits if error < 0 (NOTE: changes sign of diff)
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
