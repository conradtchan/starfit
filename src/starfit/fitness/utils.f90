module utils

  use typedef, only: &
       real64

  implicit none

contains

  function signan() result(nan)

    use, intrinsic :: &
         ieee_arithmetic, only: &
         ieee_value, &
         ieee_signaling_nan

    implicit none

    real(kind=real64) :: &
         nan

    nan = ieee_value(0.d0,ieee_signaling_nan)

  end function signan

end module utils
