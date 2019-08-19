program main
  use iso_fortran_env, only: r8 => real64
  use m_function
  use m_solve
  use m_timer
  implicit none
  type(timer) :: clock
  real(kind=r8) :: T, entropy, deltaT
  integer :: i

  call clock%start()
  do i = 1, 50
    T = real(i, kind=r8)
    call sovle_entropy(T, entropy, deltaT)
    write(*,"(F7.2,'   ',F8.4,'   ',F12.10)") T, entropy, deltaT
  end do
  call clock%stop()

end program