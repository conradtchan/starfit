module star_data

  use type_def, only: &
       real64, int64

  implicit none

  save

  integer(kind=int64) :: &
       nel, ncov
  real(kind=real64), dimension(:), allocatable :: &
       obs, err
  real(kind=real64), dimension(:, :), allocatable :: &
       cov
  integer(kind=int64) :: &
       icdf

  logical, dimension(:), allocatable :: &
       upper, covar, uncor
  integer(kind=int64) :: &
       nupper, ncovar, nuncor
  integer(kind=int64), dimension(:), allocatable :: &
       iupper, icovar, iuncor

  real(kind=real64), dimension(:, :), allocatable :: &
       m
  real(kind=real64), dimension(:), allocatable :: &
       zp, z1, z
  real(kind=real64) :: &
       mp, zz, mm


contains

  subroutine set_star_data(obs_, err_, cov_, nel_, ncov_, icdf_)

    implicit none

    integer(kind=int64), intent(in) :: &
         nel_, ncov_, icdf_
    real(kind=real64), dimension(:), intent(in) :: &
         obs_, err_
    real(kind=real64), dimension(:,:), intent(in) :: &
         cov_

    if (allocated(obs)) then
       deallocate(obs, err, cov)
    endif

    if (.not. (size(obs_, 1) == nel_)) then
       print*, '[set_star_data] obs', size(obs_, 1),  nel_
       error stop '[set_star_data] obs dimension mismatch with nel'
    endif

    if (.not. (size(err_, 1) == nel_)) then
       print*, '[set_star_data] err', size(obs_, 1),  nel_
       error stop '[set_star_data] err dimension mismatch with nel'
    endif

    if (.not. (size(cov_, 1) == nel_)) then
       print*, '[set_star_data] cov', size(cov_, 1),  nel_
       error stop '[set_star_data] cov dimension 1 mismatch with nel'
    endif

    if (.not. (size(cov_, 2) == ncov_)) then
       print*, '[set_star_data] cov', size(cov_, 2),  ncov_
       error stop '[set_star_data] cov dimension 2 mismatch with ncov'
    endif

    nel = nel_
    ncov = ncov_
    icdf = icdf_

    ! allocate implicitly

    obs = obs_
    err = err_
    cov = cov_

    call init_domains()
    call init_covariance_matrix()

  end subroutine set_star_data

  subroutine init_domains

    implicit none

    integer(kind=int64), dimension(:), allocatable :: &
         ii

    integer(kind=int64) :: &
         i

    if (allocated(upper)) then
       deallocate(upper, covar, uncor)
       deallocate(iupper, icovar, iuncor)
    endif

    allocate(upper(nel), covar(nel), uncor(nel))

    ! find masks for correlated, uncorreclated, and limit errors

    upper(:) = err < 0.d0
    covar(:) = any(cov /= 0, 2)
    uncor(:) = .not.(upper.or.covar)

    nupper = count(upper)
    ncovar = count(covar)
    nuncor = count(uncor)

    allocate(ii(nel))
    do i=1, nel
       ii(i) = i
    enddo
    iupper = pack(ii, upper)
    icovar = pack(ii, covar)
    iuncor = pack(ii, uncor)

  end subroutine init_domains

  subroutine init_covariance_matrix

    use mleqs, only: &
         leqs

    implicit none

    integer(kind=int64) :: &
         i, j, k, ie, je

    if (allocated(m)) then
       deallocate(m)
    endif

    allocate(m(ncovar, ncovar))

    do i=1, ncovar
       ie = icovar(i)
       m(i,i) = err(ie)**2
       do j=1,i
          if (j < i) then
             m(i,j) = 0.0d0
          endif
          je = icovar(j)
          do k = 1, ncov
             m(i,j) = m(i,j) + cov(ie,k) * cov(je,k)
          enddo
          if (j < i) then
             m(j,i) = m(i,j)
          endif
       enddo
    enddo

  end subroutine init_covariance_matrix

  subroutine init_covaricance_const

    use mleqs, only: &
         leqs

    implicit none

    if (allocated(zp)) then
       deallocate(zp, z1)
    endif

    allocate(zp(ncovar), z1(ncovar))

    z1(:) = 1.d0
    zp(:) = leqs(m, z1, ncovar)
    mp = sum(zp(:) * z1(:))

  end subroutine init_covaricance_const

  subroutine init_abu_covariance(abu)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         abu

    real(kind=real64), dimension(:), allocatable :: &
         v
    if (size(abu, 1) /= nel) then
       error stop '[init_abu_covariance] abu dimension mismatch'
    endif

    if (allocated(z)) then
       deallocate(z)
    endif

    allocate(z(ncovar), v(ncovar))

    v(:) = obs(icovar) - abu(icovar)

    z(:) = leqs(m, v, ncovar)
    mm = sum(v(:) * z(:))
    zz = sum(z(:))

    deallocate(v)

  end subroutine init_abu_covariance

  function abu_covariance(abu) result(xcov)

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         abu

    real(kind=real64) :: &
         xcov

    real(kind=real64), dimension(:), allocatable :: &
         diff

    if (size(abu, 1) /= nel) then
       error stop '[abu_covariance] abu dimension mismatch'
    endif

    diff = obs(icovar) - abu(icovar)
    xcov = diff_covariance(diff)

    deallocate(diff)

  end function abu_covariance

  function diff_covariance(diff) result(xcov)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64) :: &
         xcov

    real(kind=real64), dimension(:), allocatable :: &
         part1, part2

    if (size(diff, 1) /= nel) then
       error stop '[diff_covariance] diff dimension mismatch'
    endif

    part1 = diff(icovar)
    part2 = leqs(m, part1, ncovar)
    xcov = sum(part1(:) * part2(:))

    deallocate(part1, part2)

  end function diff_covariance

end module star_data
