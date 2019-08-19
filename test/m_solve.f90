module m_solve
  use iso_fortran_env, only: r8 => real64
  use m_global
  use m_function
  implicit none
  private
  public :: sovle_entropy, solve_a

  interface
    function inner(x)
      import :: r8
      real(kind=r8) :: inner
      real(kind=r8), intent(in) :: x
    end function
  end interface

  contains
  subroutine set_bisection_bound(T, lbound, ubound)
  !---------------------------------------------------------
  ! set lower bound and up bound for bisection method.
  ! values are determined with reference to
  ! J. Phys. Chem. 1994, 98, 6420−6424.
  !   T/K        a
  !   1.0   0.74552E+207
  !   2.0   0.28457E+104
  !   5.0   0.25952E+42
  !  10.0   0.52295E+21
  !  20.0   0.21157E+11
  !  50.0   10854.024
  ! 100.0   70.502286
  !    :         :
  !
  ! the relationship between T and a can be discribed by
  ! log(a) = A * T^B
  ! with A = 206.88202919, B = -0.99995677 when T < 100 K.
  ! 
  ! So when T < 100 K,
  ! log(lbound) = A * T^B - M
  ! log(ubound) = A * T^B + M
  !
  ! when T > 100 K
  ! lbound = 0, ubound = 100
  !---------------------------------------------------------
    real(kind=r8), intent(in) :: T
    real(kind=r8), intent(out) :: lbound, ubound
    real(kind=r8), parameter :: A = 206.88202919_r8, B = -0.99995677_r8
    real(kind=r8), parameter :: M = 1.5_r8

   if (T < 100._r8) then
      lbound = 10**(A * T**B - M)
      ubound = 10**(A * T**B + M)
    else
      lbound = 0._r8
      ubound = 100._r8
    end if

  end subroutine

  subroutine bisection(func, result, T)
    procedure(inner) :: func
    real(kind=r8), intent(out) :: result
    real(kind=r8), intent(in) :: T
    real(kind=r8) :: f_up, f_low, f_mid
    real(kind=r8) :: up, low, mid
    real(kind=r8) :: h
    real(kind=r8), parameter :: error = 1.E-9_r8
    integer :: counter

    call set_bisection_bound(T, low, up)

    counter = 0
    f_up = func(up)
    f_low = func(low)
    
    if ((f_up * f_low) > 0) stop "Bisection: Unreasonable guess!"

    do while (.true.)
      mid = 0.5_r8 * (up + low)
      f_mid = func(mid)

      if (f_mid == 0) then
        result = mid
        exit
      end if

      if ((f_mid * f_low) < 0) then
        up = mid
      else
        low = mid
      end if
      
      h = 0.5_r8 * (up - low)
      if (h/mid < error) then
        result = mid
        exit
      end if

      counter = counter + 1
      if (counter == 129) stop "Bisection: step exceeds"
      
    end do
    
  end subroutine

  subroutine solve_a(T, a)
    real(kind=r8), intent(in) :: T
    real(kind=r8), intent(out) :: a
    real(kind=r8) :: tmp

    ! J. Chem. Theory Comput. 2013, 9, 3165−3169
    associate (p => atm, J => degenerate, &
               h => planck_const, m => mass_e, &
               k => Boltzmann_const)

    tmp = p * h**3 / (J * (2*pi*m)**1.5 * k**2.5)

    end associate

    call bisection(f, a, T=T)

    contains
    real(kind=r8) function f(x)
      real(kind=r8), intent(in) :: x
      f = f32(x) - tmp / T**2.5
    end function
  end subroutine

  subroutine sovle_entropy(T, entropy, deltaT)
    real(kind=r8), intent(in) :: T
    real(kind=r8), intent(out) :: entropy, deltaT
    real(kind=r8) :: a
    real(kind=r8) :: calcT

    call solve_a(T, a)
    !a = 2.5456062752067030E207_r8
    !write(*,*) f32(a), f12(a)
    associate (p => atm, J => degenerate, &
               h => planck_const, m => mass_e, &
               k => Boltzmann_const)

    calcT = (p/(J*f32(a)))**0.4 * (h**2/(2*pi*m))**0.6 / k

    end associate
    deltaT = calcT - T

    entropy = gas_const * ( (5*f32(a))/(2*f12(a)) - log(a) )

  end subroutine


end module

program main
  use iso_fortran_env, only: r8 => real64
  use m_solve
  implicit none
  real(kind=r8) :: T, S, dT, a

  T = 1000._r8
  call sovle_entropy(T, S, dT)
  write(*,*) T, S, dT
  !T = 100._r8
  !call solve_a(T, a)
  !write(*,*) a
!
  !T = 50._r8
  !call solve_a(T, a)
  !write(*,*) a
!
  !T = 1._r8
  !call solve_a(T, a)
  !write(*,*) a

end program
