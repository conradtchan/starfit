module mleqs

  use type_def, only: &
       int64, real64

  contains

  function leqs(a0, b0, n) result(x)

    implicit none

    save

    integer(kind=int64), intent(in) :: &
         n
    real(kind=real64), dimension(n, n), intent(in) :: &
         a0
    real(kind=real64), dimension(n), intent(in) :: &
         b0

    real(kind=real64), dimension(n) :: &
         x

    real(kind=real64), dimension(n, n) :: &
         a

    ! this subroutine solves the  linear system a*dy=b,

    integer(kind=int64) :: &
         k,i,j,n1
    real(kind=real64) :: &
         ajji,r

    a(:,:) = a0(:,:)
    x(:) = b0(:)

    n1=n-1

    ! no conditioning
    ! reduce matrix to upper triangular form

    do j=1,n-1
       ajji = 1.D0 / a(j,j)
       do i=j+1,n
          r = -a(i,j) * ajji
          do k=j+1,n
             a(i,k) = a(i,k) + r * a(j,k)
          enddo
          x(i) = x(i) + r * x(j)
       enddo
    enddo

    ! use back substitution to find the vector solution

    x(n) = x(n) / a(n,n)

    do i=n-1, 1, -1
       r = 0.d0
       do j=i+1,n
          r = r + a(i,j) * x(j)
       enddo
       ! r = dot_product(a(i,i+1:n), x(i+1:n))
       x(i) = (x(i) - r) / a(i,i)
    enddo

  end function leqs

end module mleqs


! program test

!   use typedef, only: &
!        int32, int64, real64

!   use mleqs, only: &
!        leqs

!   implicit none

!   integer(kind=int32), parameter :: &
!        n = 3

!   real(kind=real64), dimension(n,n) :: &
!        a = reshape([1.d0,.1d0,.2d0,0.3d0,2.d0,0.4d0,0.5d0,0.6d0,3.d0], [3,3])
!   real(kind=real64), dimension(n) :: &
!        b = [1., 3., 7.], &
!        x

!   print*, 'a', a(1,:)
!   print*, '-', a(2,:)
!   print*, '-', a(3,:)
!   print*, 'b', b

!   print*,''
!   print*,'----------------'
!   print*,''

!   x = leqs(a, b, n)
!   print*, 'x', x

! end program test
