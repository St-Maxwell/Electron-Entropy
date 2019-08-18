module m_solve
  use iso_fortran_env, only: r8 => real64
  use m_function
  implicit none
  private

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
  ! with A = 211.47, B = -1.018 when T < 100 K. R^2 = 0.9997
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
    real(kind=r8), parameter :: A = 211.47_r8, B = -1.018_r8, M = 1.5_r8

    if (T < 100._r8) then
      lbound = 10**(A * T**B - M)
      ubound = 10**(A * T**B + M)
    else
      lbound = 0._r8
      ubound = 100._r8
    end if

  end subroutine

  subroutine bisection(T, KK, a)
    real(kind=r8), intent(in) :: T, KK
    real(kind=r8), intent(out) :: a
    real(kind=r8) :: temp
    real(kind=r8) :: f32_up, f32_low, f32_mid
    real(kind=r8) :: up, low, mid
    real(kind=r8) :: h
    real(kind=r8), parameter :: error = 1.E-5_r8
    integer :: counter

    temp = KK / T**2.5_r8
    
    call set_bisection_bound(T, low, up)

    counter = 0
    f32_up = f32(up) - temp
    f32_low = f32(low) - temp
    
    if ((f32_up * f32_low) > 0) stop "Bisection: Unreasonable guess!"

    do while (.true.)
      mid = 0.5_r8 * (up + mid)
      f32_mid = f32(mid) - temp
    
      if (f32_mid == 0) then
        a = mid
        exit
      end if

      if ((f32_mid * f32_low) < 0) then
        up = mid
      else
        low = mid
      end if
      
      h = 0.5_r8 * (up - low)
      if (h/mid < error) then
        a = mid
        exit
      end if

      counter = counter + 1
      if (counter == 129) stop "Bisection: step exceeds"
      
    end do
    
  end subroutine

end module