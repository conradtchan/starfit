program test_leqs

  use typedef, only: &
       int64, real64

  use linalg, only: &
       leqs, inverse, sqrtm

  implicit none

  integer(kind=int64), parameter :: &
       n = 3

  real(kind=real64), dimension(n,n) :: &
       a = reshape([1.d0,.1d0,.2d0,0.3d0,2.d0,0.4d0,0.5d0,0.6d0,3.d0], [3,3]), &
       m, q
  real(kind=real64), dimension(n) :: &
       b = [1., 3., 7.], &
       x

  print*, 'a', a(1,:)
  print*, '-', a(2,:)
  print*, '-', a(3,:)
  print*, 'b', b

  print*,''
  print*,'----------------'
  print*,''

  x = leqs(a, b, n)
  print*, 'x', x

  print*,''
  print*,'----------------'
  print*,''

  m = inverse(a, n)

  print*, 'm', m(1,:)
  print*, '-', m(2,:)
  print*, '-', m(3,:)

  x = matmul(m, b)

  print*,''
  print*,'----------------'
  print*,''

  x = leqs(a, b, n)
  print*, 'x', x

  q = matmul(m, a)

  print*,''
  print*,'----------------'
  print*,''

  print*, 'q', q(1,:)
  print*, '-', q(2,:)
  print*, '-', q(3,:)

  print*,''
  print*,'----------------'
  print*,''

  q = sqrtm(a, n)

  print*, 'q', q(1,:)
  print*, '-', q(2,:)
  print*, '-', q(3,:)

  print*,''
  print*,'----------------'
  print*,''

  q = matmul(q, q)

  print*, 'q', q(1,:)
  print*, '-', q(2,:)
  print*, '-', q(3,:)

end program test_leqs
