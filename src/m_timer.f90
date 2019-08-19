module m_timer
  use iso_fortran_env, only: r8 => real64
  use omp_lib, only: OMP_GET_WTIME
  implicit none
  private
  public :: timer

  type :: timer
    private
    real(kind=r8) :: t0
    logical :: reset = .false. ! set .true. when start
    contains
      procedure, public :: start
      procedure, public :: stop
  end type

  contains

    subroutine start(this)
      class(timer) :: this
  
      if (this%reset) stop "timer%start: start a started timer"
      this%t0 = OMP_GET_WTIME()
      this%reset = .true.
  
    end subroutine
  
    subroutine stop(this)
      class(timer) :: this
      real(kind=r8) :: t, time

      if (.not. this%reset) stop "timer%stop: calling timer%start is required before used"
      t = OMP_GET_WTIME()
      time = t - this%t0

      write(*,"(A,F10.3,A)") "Total Run Time: ", time, " s"
      this%reset = .false.

    end subroutine

end module