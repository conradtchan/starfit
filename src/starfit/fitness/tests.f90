module tests

  use fitting, only: &
       fitness

  use type_def, only: &
       real64, int64

  implicit none

contains

  subroutine test1(ls, icdf, idet, icov)

    implicit none

    integer(kind=int64), intent(in) :: &
         icdf, ls, idet, icov

    integer(kind=int64), parameter :: &
         nel = 8, &
         ncov = 2, &
         nsol = 1, &
         nstar = 1
    real(kind=real64), dimension(nel) :: &
         obs, err, det
    real(kind=real64), dimension(nel, ncov) :: &
         cov
    real(kind=real64), dimension(nsol, nstar, nel) :: &
         abu
    real(kind=real64), dimension(nstar) :: &
         c
    real(kind=real64), dimension(nsol) :: &
         f
    integer(kind=int64) :: &
         i, j

    c(:) = [0.1d0]
    cov(:,:) = 0.d0
    do i=1, nel
       obs(i) = -2.d0 - dble(i) * 0.1d0
       err(i) = + dble(i) * 0.1d0 - 0.25d0
       det(i) = - 100.d0
       abu(1,1,i) = 10.d0 ** (-2.d0 + dble(i) * 0.1d0)
    enddo

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
       det(nel-2) = obs(nel-2) - 1.d0
       det(nel-1) = obs(nel-1) + 1.d0
    endif

    print*
    print*,'[test1] ls = ', ls, 'icdf =', icdf, 'idet=', idet, 'icov=', icov
    call fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf)
    print*, 'f=', f, ', c=', c

  end subroutine test1


  subroutine suite1()

    implicit none

    integer(kind=int64), parameter :: &
         nicdf = 2, &
         nls = 3, &
         ncov = 2, &
         ndet = 2
    integer(kind=int64), dimension(nls), parameter :: &
         xls = [-1, 0, 1]
    integer(kind=int64), dimension(nicdf), parameter :: &
         xicdf = [0, 1]
    integer(kind=int64), dimension(ncov), parameter :: &
         xcov = [0, 1]
    integer(kind=int64), dimension(ndet), parameter :: &
         xdet = [0, 1]

    integer(kind=int64) :: &
         iicdf, ils, icov, idet

    do ils = 1, nls
       do iicdf = 1, nicdf
          do icov=1, ncov
             do idet=1, ndet
                call test1(xls(ils), xicdf(iicdf), xdet(idet), xcov(icov))
             end do
          end do
       end do
    end do

    ! do ils = 3, 3
    ! !do ils = 2, nls
    !    ! do iicdf = nicdf, nicdf
    !    do iicdf = 1, 1
    !       ! do icov=ncov, ncov
    !       do icov=1, 2
    !          ! do idet=ndet, ndet
    !          do idet=2, 2
    !             call test1(xls(ils), xicdf(iicdf), xdet(idet), xcov(icov))
    !          end do
    !       end do
    !    end do
    ! end do

  end subroutine suite1

end module tests
