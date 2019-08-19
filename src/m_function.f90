module m_function
  use iso_fortran_env, only: r8 => real64
  use m_global, only: pi
  implicit none
  private
  public :: f32, f12

  interface
    function inner(x)
      import :: r8
      real(kind=r8) :: inner
      real(kind=r8), intent(in) :: x
    end function
  end interface

  contains
  subroutine Simpson(func, up, down, num, int)
    use omp_lib ! OpenMP lib
    !---------------------------------------------------------------------------
    procedure(inner) :: func
    real(kind=r8), intent(in) :: up, down
    integer, intent(in) :: num   
    real(kind=r8), intent(out) :: int
    real(kind=r8) :: step, tmp
    integer :: i, tid
    !---------------------------------------------------------------------------
    ! number of threads, should manually set it to 1 when compiling codes in serial
    integer, parameter :: num_threads = 4
    
    step = 0.5_r8 * (up - down) / num
    tmp = down
    int = 0.0_r8
    
    !$OMP PARALLEL private(tmp,i,tid) reduction(+:int)
    tid = omp_get_thread_num() + 1
    do i = tid, num-1, num_threads
      tmp = step * (2*i - 1)
      int = int + 4.0_r8 * func(tmp) + 2.0_r8 * func(tmp+step)
    end do
    !$OMP END PARALLEL

    int = step * (func(down) + func(up) + 4.0_r8 * func(tmp + step) + int) / 3.0_r8

  end subroutine

  real(kind=r8) function f32(a)
    real(kind=r8), intent(in) :: a
    real(kind=r8), parameter :: N = 4._r8 / (3._r8 * sqrt(pi))
    real(kind=r8), parameter :: down = 0._r8
    real(kind=r8), parameter :: up = 0.999999_r8
    integer, parameter :: num_step = 500000
    real(kind=r8) :: int

    call Simpson(f, up, down, num_step, int)
    f32 = N * int

    contains
    pure real(kind=r8) function f(x)
      real(kind=r8), intent(in) :: x
      real(kind=r8) :: tmp
      tmp = x / (1._r8 - x)
      f = tmp * sqrt(x) / (1._r8 - x)**2.5 / (1 + exp(tmp) / a)
    end function
  end function

  real(kind=r8) function f12(a)
    real(kind=r8), intent(in) :: a
    real(kind=r8), parameter :: N = 2._r8 / sqrt(pi)
    real(kind=r8), parameter :: down = 0._r8
    real(kind=r8), parameter :: up = 0.999999_r8
    integer, parameter :: num_step = 1000000
    real(kind=r8) :: int

    call Simpson(f, up, down, num_step, int)
    f12 = N * int

    contains
    pure real(kind=r8) function f(x)
      real(kind=r8), intent(in) :: x
      real(kind=r8) :: tmp
      tmp = x / (1._r8 - x)
      f = sqrt(x) / (1._r8 - x)**2.5 / (1 + exp(tmp) / a)
    end function
  end function

end module