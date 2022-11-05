module star_data

  use type_def, only: &
       real64, int64

  implicit none

  save

  real(kind=real64), parameter :: &
       det_lim = -80.d0

  integer(kind=int64) :: &
       nel, ncov
  real(kind=real64), dimension(:), allocatable :: &
       obs, err, det
  real(kind=real64), dimension(:, :), allocatable :: &
       cov
  integer(kind=int64) :: &
       icdf

  logical, dimension(:), allocatable :: &
       upper, covar, uncor, measu, detec, nocov, erinv
  integer(kind=int64) :: &
       nupper, ncovar, nuncor, nmeasu, ndetec, nnocov, nerinv
  integer(kind=int64), dimension(:), allocatable :: &
       iupper, icovar, iuncor, imeasu, idetec, inocov, ierinv

  real(kind=real64), dimension(:, :), allocatable :: &
       mm
  real(kind=real64), dimension(:), allocatable :: &
       zvp, zv1, zv, &
       erri, erri2
  real(kind=real64) :: &
       mp, z, m

contains

  subroutine set_star_data(obs_, err_, det_, cov_, nel_, ncov_, icdf_)

    implicit none

    integer(kind=int64), intent(in) :: &
         nel_, ncov_, icdf_
    real(kind=real64), dimension(:), intent(in) :: &
         obs_, err_, det_
    real(kind=real64), dimension(:,:), intent(in) :: &
         cov_

    if (allocated(obs)) then
       deallocate(obs, err, det, cov)
    endif

    if (.not. (size(obs_, 1) == nel_)) then
       print*, '[set_star_data] obs', size(obs_, 1),  nel_
       error stop '[set_star_data] obs dimension mismatch with nel'
    endif

    if (.not. (size(err_, 1) == nel_)) then
       print*, '[set_star_data] err', size(err_, 1),  nel_
       error stop '[set_star_data] err dimension mismatch with nel'
    endif

    if (.not. (size(det_, 1) == nel_)) then
       print*, '[set_star_data] det', size(det_, 1),  nel_
       error stop '[set_star_data] det dimension mismatch with nel'
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
    det = det_
    cov = cov_

    call init_domains()
    call init_erri()
    call init_covariance_matrix()

  end subroutine set_star_data

  subroutine init_domains

    implicit none

    integer(kind=int64), dimension(:), allocatable :: &
         ii

    integer(kind=int64) :: &
         i

    if (allocated(upper)) then
       deallocate(upper, covar, uncor, measu, detec, nocov, erinv)
       deallocate(iupper, icovar, iuncor, imeasu, idetec,inocov, ierinv)
    endif

    ! find masks for correlated, uncorreclated, and limit errors

    upper = err < 0.d0
    covar = any(cov /= 0, 2)
    uncor = .not.(upper.or.covar)
    measu = .not.upper
    detec = measu(:) .and. (det(:) > det_lim)
    nocov = .not.covar
    erinv = upper.or.detec.or.uncor

    nupper = count(upper)
    ncovar = count(covar)
    nuncor = count(uncor)
    nmeasu = nel - nupper
    ndetec = count(detec)
    nnocov = nel - ncov
    nerinv = count(erinv)

    allocate(ii(nel))
    do i=1, nel
       ii(i) = i
    enddo
    iupper = pack(ii, upper)
    icovar = pack(ii, covar)
    iuncor = pack(ii, uncor)
    imeasu = pack(ii, measu)
    idetec = pack(ii, detec)
    inocov = pack(ii, nocov)
    ierinv = pack(ii, erinv)

  end subroutine init_domains


  subroutine init_covariance_matrix

    use mleqs, only: &
         leqs

    implicit none

    integer(kind=int64) :: &
         i, j, k, ie, je

    if (allocated(mm)) then
       deallocate(mm)
    endif

    allocate(mm(ncovar, ncovar))

    do i=1, ncovar
       ie = icovar(i)
       mm(i,i) = err(ie)**2
       do j=1,i
          if (j < i) then
             mm(i,j) = 0.0d0
          endif
          je = icovar(j)
          do k = 1, ncov
             mm(i,j) = mm(i,j) + cov(ie,k) * cov(je,k)
          enddo
          if (j < i) then
             mm(j,i) = mm(i,j)
          endif
       enddo
    enddo

  end subroutine init_covariance_matrix


  subroutine init_covaricance_const

    use mleqs, only: &
         leqs

    implicit none

    if (allocated(zvp)) then
       deallocate(zvp, zv1)
    endif

    allocate(zvp(ncovar), zv1(ncovar))

    zv1(:) = 1.d0
    zvp(:) = leqs(mm, zv1, ncovar)
    mp = sum(zvp(:) * zv1(:))

  end subroutine init_covaricance_const

  subroutine init_nocov_erri2

    implicit none

    if (allocated(erri2)) then
       deallocate(erri2)
    endif

    allocate(erri2(nel))

    erri2(inocov) = 1.d0 / err(inocov)**2

  end subroutine init_nocov_erri2


  subroutine init_erri

    implicit none

    if (allocated(erri)) then
       deallocate(erri)
    endif

    allocate(erri(nel))

    erri(ierinv) = 1.d0 / err(ierinv)

  end subroutine init_erri


  subroutine init_abu_covariance(abu)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         abu

    real(kind=real64), dimension(:), allocatable :: &
         vv
    if (size(abu, 1) /= nel) then
       error stop '[init_abu_covariance] abu dimension mismatch'
    endif

    if (allocated(zv)) then
       deallocate(zv)
    endif

    allocate(zv(ncovar), vv(ncovar))

    vv(:) = abu(icovar) - obs(icovar)

    zv(:) = leqs(mm, vv, ncovar)
    mm = sum(vv(:) * zv(:))
    z = sum(zv(:))

    deallocate(vv)

  end subroutine init_abu_covariance


  subroutine init_abu_z(abu)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         abu

    if (size(abu, 1) /= nel) then
       error stop '[init_abu_covariance] abu dimension mismatch'
    endif

    z = sum(zvp(:) * abu(icovar))

  end subroutine init_abu_z


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
    part2 = leqs(mm, part1, ncovar)
    xcov = sum(part1(:) * part2(:))

    deallocate(part1, part2)

  end function diff_covariance


  function diff_z(diff) result(xz)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64) :: &
         xz

    if (size(diff, 1) /= nel) then
       error stop '[diff_covariance] diff dimension mismatch'
    endif

    xz = sum(zvp(:) * diff(icovar))

  end function diff_z

end module star_data
