module tests

  use fitting, only: &
       fitness, fitness_m

  use typedef, only: &
       real64, int64

  implicit none

contains

  subroutine test1(ls, icdf, idet, icov, nstar, flags)

    implicit none

    integer(kind=int64), intent(in) :: &
         icdf, ls, idet, icov, nstar, flags

    integer(kind=int64), parameter :: &
         nel = 8, &
         ncov = 2, &
         nsol = 1
    real(kind=real64), dimension(nel) :: &
         obs, err, det
    real(kind=real64), dimension(nel, ncov) :: &
         cov
    real(kind=real64), dimension(nsol, nstar, nel) :: &
         abu
    real(kind=real64), dimension(nsol, nstar) :: &
         c
    real(kind=real64), dimension(nsol) :: &
         f
    integer(kind=int64) :: &
         i, j

    c(:,:) = 0.1d0
    cov(:,:) = 0.d0
    do i=1, nel
       obs(i) = -2.d0 - (dble(i) - 3.5d0)**2 * 0.1d0
       err(i) = + dble(i) * 0.1d0 - 0.25d0
       det(i) = - 100.d0
       do j=1,nstar
          if (j == 1) then
             abu(1,j,i) = 10.d0 ** (-2.d0 + dble(i) * 0.1d0)
          else if (j == 2) then
             abu(1,j,i) = 10.d0 ** (-1.d0 - dble(i) * 0.1d0)
          else if (j == 3) then
             abu(1,j,i) = 10.d0 ** (-1.d0 - (dble(i) - 0.5d0*dble(nel))**2 * 0.1d0)
          else
             abu(1,j,i) = 10.d0 ** (-2.d0 + dble(i) * 0.1d0) * (1.1d0 + sin(dble((j-1)*i)))
          endif
       enddo
    enddo

    ! do j=1, nstar
    !    abu(1,j,j+3) = 1.d-80
    ! end do

    if (icov > 0) then
       do i = nel-4, nel-2
          do j = 1, ncov
             cov(i,j) = err(i) * dble(i * j) / dble(ncov * nel) * dble(1 - 2 * mod(i + j, 2))
          enddo
          !print*,'[err]', i, err(i)
          err(i) = sqrt(err(i)**2 - sum(cov(i, :)**2))
          !print*,'[err]', i, err(i)
       enddo

       !do i = 1, nel
       !   print*,i, err(i), cov(i,:)
       !enddo

    endif

    if (idet > 0) then
       det(nel) = obs(nel) - 1.d0
       det(nel-1) = obs(nel-1) + 1.d0
       det(nel-2) = obs(nel-2) - 1.d0
       det(nel-3) = obs(nel-3) - 2.d0
    endif

    print*
    print*,'[test1] ls = ', ls, 'icdf =', icdf, 'idet=', idet, 'icov=', icov, 'flags=', flags, 'nstar=', nstar
    call fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf, flags)
    print*, 'f=', f, ', c=', c

  end subroutine test1


  subroutine suite1()

    implicit none

    integer(kind=int64), parameter :: &
         nicdf = 2, &
         nls = 3, &
         ncov = 2, &
         ndet = 2, &
         nstar = 3, &
         nflags = 4
    integer(kind=int64), dimension(nls), parameter :: &
         xls = [-1, 0, 1]
    integer(kind=int64), dimension(nicdf), parameter :: &
         xicdf = [0, 1]
    integer(kind=int64), dimension(ncov), parameter :: &
         xcov = [0, 1]
    integer(kind=int64), dimension(ndet), parameter :: &
         xdet = [0, 1]
    integer(kind=int64), dimension(nstar), parameter :: &
         xstar = [1,2,3]
    integer(kind=int64), dimension(nflags), parameter :: &
         xflags = [0,1,2,3]

    integer(kind=int64) :: &
         iicdf, ils, icov, idet, istar, iflags

    do istar = 1, nstar
       do iflags=1, nflags
          do ils = 1, nls
             do iicdf = 1, nicdf
                do icov=1, ncov
                   do idet=1, ndet
                      call test1(xls(ils), xicdf(iicdf), xdet(idet), xcov(icov), xstar(istar), xflags(iflags))
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine suite1


  subroutine suite0()

    implicit none

    integer(kind=int64) :: &
         iicdf, ils, icov, idet, istar, iflags

    do istar = 3,3
       do iflags=1,1
          do ils = 1, 2
             do iicdf = 1, 1
                do icov = 1, 1
                   do idet=1, 1
                      call test1(ils, iicdf, idet, icov, istar, iflags)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine suite0




end module tests
