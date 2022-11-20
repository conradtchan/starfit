module star_data

  use typedef, only: &
       real64, int64

  implicit none

  save

  ! use of inverse may be less accurate but a lot faster if comparing
  ! to many models.

  logical, parameter :: &
       use_inverse = .true., &
       warn_subthreshold_detection = .false., &
       stop_subthreshold_detection = .false., &
       warn_zero_error = .true., &
       stop_zero_error = .true.

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
       upper, covar, uncor, measu, detec, nocov, erinv, ernoi
  integer(kind=int64) :: &
       nupper, ncovar, nuncor, nmeasu, ndetec, nnocov, nerinv, nernoi
  integer(kind=int64), dimension(:), allocatable :: &
       iupper, icovar, iuncor, imeasu, idetec, inocov, ierinv, iernoi

  real(kind=real64), dimension(:, :), allocatable :: &
       mm, mm1
  real(kind=real64), dimension(:), allocatable :: &
       zvp, zv, &
       ert, eri, ei2
  real(kind=real64) :: &
       mp

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

    call init_check_errors()
    call init_domains()
    call init_ert()
    call init_eri()
    call init_covariance_matrix()
    call init_inverse()
    call init_check_thresholds()

  end subroutine set_star_data


  subroutine init_check_thresholds()

    implicit none

    integer(kind=int64) :: &
         i, i1

    if (any(det(idetec) > obs(idetec))) then
       if (warn_subthreshold_detection) then
          do i1 = 1, ndetec
             i = idetec(i1)
             if (det(i) > obs(i)) then
                print*,'[set_star_data] WARNING i=',i,'det=',det(i), 'obs=', obs(i)
             endif
          end do
       endif

       if (stop_subthreshold_detection) then
          error stop '[set_star_data] observation below detection limit'
       endif
    end if

  end subroutine init_check_thresholds


  subroutine init_check_errors()

    use utils, only: &
         signan

    implicit none

    integer(kind=int64) :: &
         i

    ! check errors

    if (any(err == 0.d0)) then
       if (warn_zero_error) then
          do i = 1, nel
             if (err(i) == 0.d0) then
                print*,'[init_check_error] WARNING i=',i,'err=',err(i), 'obs=', obs(i)
                err(i) = signan()
             endif
          end do
       endif
       if (stop_zero_error) then
          error stop '[init_check_error] uncorrelated errors must not be zero'
       endif
    endif

  end subroutine init_check_errors


  subroutine init_domains

    implicit none

    integer(kind=int64), dimension(:), allocatable :: &
         ii
    integer(kind=int64) :: &
         i

    if (allocated(upper)) then
       deallocate(upper, covar, uncor, measu, detec, nocov, erinv, ernoi)
       deallocate(iupper, icovar, iuncor, imeasu, idetec,inocov, ierinv, iernoi)
    endif

    ! find masks for correlated, uncorreclated, and limit errors

    upper = err < 0.d0
    measu = .not.upper
    covar = any(cov /= 0.d0, 2).and.measu
    uncor = .not.(upper.or.covar)
    ! detec = measu(:).and.(det(:) > det_lim)
    ! temporary fix
    detec = measu(:).and.(det(:) > det_lim).and.uncor
    nocov = .not.covar
    erinv = upper.or.detec.or.uncor
    ernoi = .not.erinv

    nupper = count(upper)
    ncovar = count(covar)
    nuncor = count(uncor)
    nmeasu = nel - nupper
    ndetec = count(detec)
    nnocov = nel - ncovar
    nerinv = count(erinv)
    nernoi = nel - nerinv

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
    iernoi = pack(ii, ernoi)

  end subroutine init_domains


  subroutine init_covariance_matrix

    use mleqs, only: &
         leqs

    implicit none

    integer(kind=int64) :: &
         i, j, k, ie, je

    if (ncovar == 0) then
       return
    endif

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


  subroutine init_inverse

    use mleqs, only: &
         inverse

    implicit none

    if (.not.use_inverse) then
       return
    endif

    if (ncovar == 0) then
       return
    endif

    if (.not.allocated(mm)) then
       error stop '[init_inverse] matrix mm not present'
       return
    endif

    if (allocated(mm1)) then
       deallocate(mm1)
    endif

    mm1 = inverse(mm, ncovar)

  end subroutine init_inverse


  subroutine init_covaricance_const

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), allocatable :: &
         zv1

    if (ncovar == 0) then
       mp = 0.d0
       return
    endif

    if (allocated(zvp)) then
       deallocate(zvp)
    endif

    allocate(zvp(ncovar))
    if (use_inverse) then
       zvp(:) = sum(mm1(:,:), 2)
    else
       allocate(zv1(ncovar))
       zv1(:) = 1.d0
       zvp(:) = leqs(mm, zv1, ncovar)
    endif
    mp = sum(zvp(:))

  end subroutine init_covaricance_const

  subroutine init_ei2

    use utils, only: &
         signan

    implicit none

    if (allocated(ei2)) then
       deallocate(ei2)
    endif

    allocate(ei2(nel))

    ei2(iernoi) = signan()
    ei2(ierinv) = 1.d0 / ert(ierinv)**2

  end subroutine init_ei2


  subroutine init_eri

    use utils, only: &
         signan

    implicit none

    if (allocated(eri)) then
       deallocate(eri)
    endif

    allocate(eri(nel))

    eri(iernoi) = signan()
    eri(ierinv) = 1.d0 / ert(ierinv)

  end subroutine init_eri


  subroutine init_ert

    use utils, only: &
         signan

    implicit none

    if (allocated(ert)) then
       deallocate(ert)
    endif

    allocate(ert(nel))

    if (nupper > 0) then
       ert(iupper) = err(iupper)
    endif
    if (ncovar > 0) then
       ert(icovar) = sqrt(sum(cov(icovar,:)**2, 2) + err(icovar)**2)
    endif
    if (nuncor > 0) then
       ert(iuncor) = err(iuncor)
    endif

  end subroutine init_ert


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

    if (ncovar == 0) then
       xcov = 0.d0
       return
    endif

    allocate(diff(ncovar))

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

    real(kind=real64), dimension(:), allocatable  :: &
         part1, part2

    if (size(diff, 1) /= nel) then
       print*, '[diff_covariance] size(diff, 1) = ', size(diff, 1), &
            'expected nel = ', nel
       error stop '[diff_covariance] diff dimension mismatch'
    endif

    if (ncovar == 0) then
       xcov = 0.d0
       return
    endif

    allocate(part1(ncovar), part2(ncovar))
    part1 = diff(icovar)
    if (use_inverse) then
       part2 = matmul(mm1, part1)
    else
       part2 = leqs(mm, part1, ncovar)
    endif
    xcov = sum(part1(:) * part2(:))

    deallocate(part1, part2)

  end function diff_covariance


  function diff_covariance_m(diff) result(xcov)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64), dimension(ncovar, ncovar) :: &
         xcov

    real(kind=real64), dimension(:), allocatable  :: &
         part1, part2

    integer(kind=int64) :: &
         i

    if (size(diff, 1) /= nel) then
       print*, '[diff_covariance] size(diff, 1) = ', size(diff, 1), &
            'expected nel = ', nel
       error stop '[diff_covariance] diff dimension mismatch'
    endif

    if (ncovar == 0) then
       return
    endif

    allocate(part1(ncovar))
    part1 = diff(icovar)
    if (use_inverse) then
       do i=1, ncovar
          xcov(i,:) = mm1(i,:) * part1(i) * part1(:)
       enddo
    else
       allocate(part2(ncovar))
       do i=1, ncovar
          part2(:) = 0.d0
          part2(i) = part1(i)
          part2(:) = leqs(mm, part2, ncovar)
          xcov(i,:) = part2(:) * part1(:)
       enddo
       deallocate(part2)
    endif

    deallocate(part1)

  end function diff_covariance_m


  function diff_z(diff) result(xz)

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64) :: &
         xz

    if (ncovar == 0) then
       xz = 0.d0
       return
    endif

    if (size(diff, 1) /= nel) then
       print*, '[diff_z] size(diff, 1) = ', size(diff, 1), &
            'expected nel = ', nel
       error stop '[diff_z] diff dimension mismatch'
    endif

    xz = sum(zvp(:) * diff(icovar))

  end function diff_z


  function diff_zv(diff) result(xzv)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64), dimension(ncovar) :: &
         xzv

    real(kind=real64), dimension(ncovar)  :: &
         part

    if (ncovar == 0) then
       return
    endif

    if (size(diff, 1) /= nel) then
       print*, '[diff_zv] size(diff, 1) = ', size(diff, 1), &
            'expected nel = ', nel
       error stop '[diff_zv] diff dimension mismatch'
    endif

    part(:) = diff(icovar)
    if (use_inverse) then
       xzv(:) = matmul(mm1, part)
    else
       xzv(:) = leqs(mm, part, ncovar)
    endif

  end function diff_zv


  function reduced_diff_zv(diff) result(xzv)

    use mleqs, only: &
         leqs

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         diff

    real(kind=real64), dimension(ncovar) :: &
         xzv

    if (ncovar == 0) then
       return
    endif

    if (size(diff, 1) /= ncovar) then
       print*, '[reduced_diff_zv] size(diff, 1) = ', size(diff, 1), &
            'expected ncovar = ', ncovar
       error stop '[reduced_diff_zv] diff dimension mismatch'
    endif

    if (use_inverse) then
       xzv(:) = matmul(mm1, diff)
    else
       xzv(:) = leqs(mm, diff, ncovar)
    endif

  end function reduced_diff_zv

end module star_data
