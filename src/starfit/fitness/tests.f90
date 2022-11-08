module tests

  use fitting, only: &
       fitness

  use type_def, only: &
       real64, int64

  implicit none

contains

  subroutine test1(ls, icdf, idet, icov, nstar)

    implicit none

    integer(kind=int64), intent(in) :: &
         icdf, ls, idet, icov, nstar

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
       det(nel-2) = obs(nel-2) - 1.d0
       det(nel-1) = obs(nel-1) + 1.d0
    endif

    print*
    print*,'[test1] ls = ', ls, 'icdf =', icdf, 'idet=', idet, 'icov=', icov, 'nstar=', nstar
    call fitness(f, c, obs, err, det, cov, abu, nel, ncov, nstar, nsol, ls, icdf)
    print*, 'f=', f, ', c=', c

  end subroutine test1


  subroutine suite1()

    implicit none

    integer(kind=int64), parameter :: &
         nicdf = 2, &
         nls = 3, &
         ncov = 2, &
         ndet = 2, &
         nstar = 3
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

    integer(kind=int64) :: &
         iicdf, ils, icov, idet, istar

    do istar = 1, nstar
       do ils = 1, nls
          do iicdf = 1, nicdf
             do icov=1, ncov
                do idet=1, ndet
                   call test1(xls(ils), xicdf(iicdf), xdet(idet), xcov(icov), xstar(istar))
                end do
             end do
          end do
       end do
    enddo

  end subroutine suite1


  subroutine suite0()

    implicit none

    integer(kind=int64) :: &
         iicdf, ils, icov, idet, istar

    do istar = 3,3
       do ils = 1, 2
          do iicdf = 1, 1
             do icov = 1, 1
                do idet=1, 1
                   call test1(ils, iicdf, idet, icov, istar)
                end do
             end do
          end do
       end do
    enddo

  end subroutine suite0

end module tests
! f=   14.604249359153073      , c=   7.0420336228949054E-011   2.2653238429537481E-005   3.2096320643436704E-002
! f=   14.604249309232987      , c=   6.7158066014241258E-017   2.2653544306922120E-005   3.2096320469292339E-002
! f=   14.604133232385269      , c=   3.4405657758547667E-008   3.8946170688527796E-007   3.2101178367944717E-002

! f=   12.165745463860274      c=   4.6247117246878133E-011   7.2620212621110625E-008  0.11024695364678094
! f=   12.165735642608530      , c=   1.1541906478024271E-009   2.6135607039924539E-009  0.11025437963729502
! f=   12.165735784207287      , c=   1.3863142962285896E-009   2.3905531076758381E-009  0.11025437975892589
! f=   12.165734427445921      , c=   3.7751357595539048E-011   3.7751357595539048E-011  0.11025438112159469

! f=   3.1795739026386856      c=   1.1584081693882808E-008   1.3881175542351798E-008  0.11062983728596676
! f=   3.1795679850030267      , c=   3.7751357595539048E-011   3.7751357595539048E-011  0.11062966038664263
! f=   3.1795679850030218      , c=   3.7751357595539048E-011   3.7751357595539048E-011  0.11062966038664268
