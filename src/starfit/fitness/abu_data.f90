module abu_data

  use typedef, only: &
       real64, int64

  implicit none

  save

  public

  integer(kind=int64) :: &
       nel, nstar
  real(kind=real64), dimension(:, :), allocatable :: &
       abu

  real(kind=real64), dimension(:), allocatable :: &
       abumix, logabumix, logabu

  real(kind=real64), parameter, private :: &
       ln10 = log(10.0d0), &
       ln10i = 1.d0 / ln10

contains

  subroutine set_abu_data(abu_, nel_, nstar_)

    implicit none

    integer(kind=int64), intent(in) :: &
         nel_, nstar_
    real(kind=real64), dimension(..), intent(in), target :: &
         abu_

    if (allocated(abu)) then
       deallocate(abu)
    endif

    if (rank(abu_) == 2) then
       if (.not. (size(abu_, 1) == nstar_)) then
          print*, '[set_abu_data] abu', size(abu_, 1),  nstar_
          error stop '[set_abu_data] abu dimension 1 mismatch with nstar'
       endif
       if (.not. (size(abu_, 2) == nel_  )) then
          print*, '[set_abu_data] abu', size(abu_, 2),  nel_
          error stop '[set_abu_data] abu dimension 2 mismatch with nel'
       endif

    else if (rank(abu_) == 1) then
       if (.not. (1 == nstar_)) then
          print*, '[set_abu_data] nstar', nstar_
          error stop '[set_abu_data] nstar needs to be 1 for abu rank 1'
       endif
       if (.not. (size(abu_, 1) == nel_  )) then
          print*, '[set_abu_data] abu', size(abu_, 2),  nel_
          error stop '[set_abu_data] abu dimension 2 mismatch with nel'
       endif
    else
       print*, '[set_abu_data] rank(abu_)', rank(abu_)
       error stop '[set_abu_data] rank(abu_) needs to be 1 or 2.'
    endif

    nel = nel_
    nstar = nstar_

    select rank(abu_)
       rank(1)
          allocate(abu(1,nel))
          abu(1,:) = abu_(:)
       rank(2)
          allocate(abu(nstar,nel))
          abu(:,:) = abu_(:,:)
       rank default
          error stop '[set_abu_data] bad rank (you should never see this message)'
    end select

  end subroutine set_abu_data


  subroutine init_logabu()

    implicit none

    if (nstar /= 1) then
       print*, '[init_logabu] nstar', nstar
       error stop '[init_abu] nstar should be 1'
    endif

    if (allocated(logabu)) then
       deallocate(logabu)
    endif

    ! allocate implicitly

    logabu = log(abu(1,:)) * ln10i

  end subroutine init_logabu


  subroutine init_abu_mix(offsets)

    implicit none

    real(kind=real64), dimension(:), intent(in) :: &
         offsets

    integer(kind=int64) :: &
         i

    if (allocated(abumix)) then
       deallocate(abumix, logabumix)
    endif

    allocate(abumix(nel))

    do i = 1, nel
       abumix(i) = sum(offsets(:) * abu(:,i))
    enddo

    ! allocate implicitly

    logabumix = log(abumix) * ln10i

  end subroutine init_abu_mix

end module abu_data
