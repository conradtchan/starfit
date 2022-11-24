module linalg

  use typedef, only: &
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

    ! this subroutine solves the linear system a*dy=b,

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

    ! this subroutine finds the inverse of matrix a

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


  function sqrtm(a, n) result(x1)

    use utils, only: &
         signan

    implicit none

    save

    integer(kind=int64), parameter :: &
         itmax = 100
    real(kind=real64), parameter :: &
         tol = 1.d-12, &
         tol2 = 1.d-8

    integer(kind=int64), intent(in) :: &
         n
    real(kind=real64), dimension(n, n), intent(in) :: &
         a

    real(kind=real64), dimension(n, n) :: &
         x1

    real(kind=real64), dimension(n, n) :: &
         x0

    ! this subroutine finds the square root of a matrix

    integer(kind=int64) :: &
         i
    real(kind=real64) :: &
         d

    if (n == 1) then
       x1(1,1) = sqrt(a(1,1))
       return
    endif

    x0(:,:) = 0.d0
    do i = 1, n
       x0(i,i) = 1.d0
    end do

    x0(:,:) = 0.5d0 * (x0(:,:) + a(:,:))

    do i=1, itmax
       x1(:,:) = 0.5d0 * (x0(:,:) + matmul(a(:,:),inverse(x0(:,:), n)))
       d = sqrt(sum((x0(:,:) - x1(:,:))**2) / sum(x0(:,:)**2))
       if (d < tol) goto 1000
       x0(:,:) = x1(:,:)
    end do

    d = sqrt(sum((matmul(x1, x1) - a(:,:))**2) / sum(a(:,:)**2))
    if (d < tol2) goto 1000

    ! print*, '[sqrtm]', d
    ! error stop '[sqrtm] No convergence.'
    x1(:,:) = signan()
    return

1000 continue

  end function sqrtm

end module linalg
