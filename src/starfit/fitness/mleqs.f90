module mleqs

  use type_def, only: &
       int64, real64

  contains

  function leqs(a_, b_, n) result(x)

    implicit none

    save

    integer(kind=int64), intent(in) :: &
         n
    real(kind=real64), dimension(n, n), intent(in) :: &
         a_
    real(kind=real64), dimension(n), intent(in) :: &
         b_

    real(kind=real64), dimension(n) :: &
         x

    real(kind=real64), dimension(n, n) :: &
         a

    ! this subroutine solves the  linear system a*dy=b,

    integer(kind=int64) :: &
         k,i,j,n1
    real(kind=real64) :: &
         ajji,r

    if (n == 1) then
       x(1) = b_(1) / a_(1,1)
       return
    endif

    a(:,:) = a_(:,:)
    x(:) = b_(:)

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


  function inverse(a_, n) result(x)

    implicit none

    save

    integer(kind=int64), intent(in) :: &
         n
    real(kind=real64), dimension(n, n), intent(in) :: &
         a_

    real(kind=real64), dimension(n, n) :: &
         x

    real(kind=real64), dimension(n, n) :: &
         a

    ! this subroutine solves the  linear system a*dy=b,

    integer(kind=int64) :: &
         k,i,j,l,n1
    real(kind=real64) :: &
         ajji,r

    if (n == 1) then
       x(1,1) = 1.d0 / a_(1,1)
       return
    endif

    a(:,:) = a_(:,:)
    x(:,:) = 0.d0
    do j = 1, n
       x(j,j) = 1.d0
    end do

    n1 = n-1

    ! no conditioning
    ! reduce matrix to upper triangular form

    do j=1,n-1
       ajji = 1.D0 / a(j,j)
       do i=j+1,n
          r = -a(i,j) * ajji
          do k=j+1,n
             a(i,k) = a(i,k) + r * a(j,k)
          enddo
          do k=1,n
             x(i,k) = x(i,k) + r * x(j,k)
          enddo
       end do
    end do

    ! use back substitution to find the vector solution

    do l=1, n
       x(n,l) = x(n,l) / a(n,n)
    end do

    do l=1, n
       do i=n-1, 1, -1
          r = 0.d0
          do j=i+1,n
             r = r + a(i,j) * x(j,l)
          end do
          ! r = dot_product(a(i,i+1:n), x(i+1:n))
          x(i,l) = (x(i,l) - r) / a(i,i)
       end do
    end do

  end function inverse

end module mleqs


! program test_leqs

!   use type_def, only: &
!        int64, real64

!   use mleqs, only: &
!        leqs, invert

!   implicit none

!   integer(kind=int64), parameter :: &
!        n = 3

!   real(kind=real64), dimension(n,n) :: &
!        a = reshape([1.d0,.1d0,.2d0,0.3d0,2.d0,0.4d0,0.5d0,0.6d0,3.d0], [3,3]), &
!        m, q
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

!   print*,''
!   print*,'----------------'
!   print*,''

!   m = invert(a, n)

!   print*, 'm', m(1,:)
!   print*, '-', m(2,:)
!   print*, '-', m(3,:)

!   x = matmul(m, b)

!   print*,''
!   print*,'----------------'
!   print*,''

!   x = leqs(a, b, n)
!   print*, 'x', x

!   q = matmul(m, a)

!   print*,''
!   print*,'----------------'
!   print*,''

!   print*, 'q', q(1,:)
!   print*, '-', q(2,:)
!   print*, '-', q(3,:)

! end program test_leqs
