module abu_data

  use type_def, only: &
       real64, int64

  implicit none

  save

  integer(kind=int64) :: &
       nel, nstar
  real(kind=real64), dimension(:, :), allocatable :: &
       abu

contains

  subroutine set_abu_data(abu_, nel_, nstar_)

    implicit none

    integer(kind=int64), intent(in) :: &
         nel_, nstar_
    real(kind=real64), dimension(:,:), intent(in) :: &
         abu_

    if (allocated(abu)) then
       deallocate(abu)
    endif

    if (.not. (size(abu_, 1) == nstar_)) then
       print*, '[set_abu_data] abu', size(abu_, 1),  nstar_
       error stop '[set_abu_data] abu dimeension 1 mismatch with nstar'
    endif
    if (.not. (size(abu_, 2) == nel_  )) then
       print*, '[set_abu_data] abu', size(abu_, 2),  nel_
       error stop '[set_abu_data] abu dimeension 2 mismatch with nel'
    endif

    nel = nel_
    nstar = nstar_

    ! allocate implicitly

    abu = abu_

  end subroutine set_abu_data

end module abu_data
