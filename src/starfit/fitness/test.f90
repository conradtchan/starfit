program test

  use type_def, only: &
       real64, int64

  use tests, only: &
       suite1, suite0

  implicit none

  call suite1()


end program test
