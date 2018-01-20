PROGRAM electron_entropy
    implicit none
    real ( kind = 8 ) :: T = 250.0D0
    real ( kind = 8 ) :: T0 = 0.0D0
    real ( kind = 8 ) :: step = 5.0D1
    integer :: i
    write(*,*) "T(K)     S(J/mol·K)  A"
    
    do i = 1, 100
        T = T0 + step * dble( i )
        call calcS ( T )
    end do
    
    stop
END PROGRAM

subroutine calcS ( T )
    implicit none
    real ( kind = 8 ), parameter :: h = 6.6260699D0 ! D-34
    real ( kind = 8 ), parameter :: kB = 1.3806485D0 ! D-23
    real ( kind = 8 ), parameter :: mass_e = 9.10938356D0 ! D-31
    real ( kind = 8 ), parameter :: pi = 3.14159265358979D0
    real ( kind = 8 ), parameter :: J = 2.0D0 ! degree of degeneration
    real ( kind = 8 ), parameter :: atm = 1.0D0 ! D5
    real ( kind = 8 ), parameter :: R = 8.314D0
    real ( kind = 8 ) :: fp32, fp12, KK, T, Sm
    real ( kind = 8 ) :: a 
    
    KK = ( atm * h**3.0D0 * 1.0D7 ) / ( J * ( 2.0D0 * pi * mass_e )**1.5D0 * kB**2.5D0 )
    call bisection ( T, KK, a )

    call calcfp( a, 32, fp32 )
    call calcfp( a, 12, fp12 )
    
    T = ( KK / fp32 )**0.4D0
    Sm = R * ( ( 5.0D0 * fp32 ) / ( 2.0D0 * fp12 ) - log( a ) )

    write(*,"(F7.2,'   ',F8.4,'   ',F10.6)") T, Sm, a

end subroutine

subroutine calcfp( a, p, res )
  implicit none
  external func
  integer , parameter :: n = 50000
  real ( kind = 8 ) , parameter :: low = 0.0D0
  real ( kind = 8 ) , parameter :: up = 6.0D1
  real ( kind = 8 ) :: a
  integer :: p
  real ( kind = 8 ) :: res

  call simpson ( func, res, low, up, n, a, p )

  select case ( p )
    case ( 12 )
      res = res * 1.128379167D0
    case ( 32 )
      res = res * 0.752252778D0
    case default
      stop
  end select

end subroutine

subroutine func ( f, x, a, p )
  implicit none
  real ( kind = 8 ) :: f, x, a
  integer :: p

  select case ( p )
    case ( 12 )
      f = x**0.5D0 / ( 1.0D0 + exp ( x ) / a )
    case ( 32 )
      f = x**1.5D0 / ( 1.0D0 + exp ( x ) / a )
    case default
      stop
  end select

end subroutine

subroutine simpson (func,s,low,up,n,a,p)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  复合Simpson公式计算数值积分
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.         func 为外部子程序
!       2.         
!       3.         a,b积分区间
!       4.         n 区间划分个数
!  Output parameters  :
!       1.         s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
  implicit none
  external func
  integer :: n, k, p
  real ( kind = 8 ) :: a, s
  real ( kind = 8 ) :: up, low, h, t1, t2
  real ( kind = 8 ) :: f1, f2, f3, f4

  s = 0D0
  h = ( up - low ) / n / 2d0
  
  call func ( f1, low, a, p )
  call func ( f2, up, a, p )
  
  s = f1 + f2
  
  !k=0 情况
  call func ( f1, low + h, a, p )
  s = s + 4D0 * f1
  
  do k = 1, n-1
  
    t1 = low + ( 2D0 * k + 1 ) * h
    t2 = low + 2D0 * k * h
    
    call func ( f3, t1, a, p )
    call func ( f4, t2, a, p )
    
    s = s + f3 * 4D0 + f4 * 2D0  
  
  end do
  
  s = s * h / 3D0

end subroutine

subroutine bisection( T, KK, a )
    implicit none
    real ( kind = 8 ), intent ( in ) :: T, KK
    real ( kind = 8 ) :: temp
    real ( kind = 8 ) :: a, a_up, a_low, a_mid, hda
    real ( kind = 8 ) :: f32_up, f32_low, f32_mid
    real ( kind = 8 ), parameter :: error = 1.0D-9

    temp = KK / T**2.5D0
    
    select case ( ceiling( T ) )
      case ( 100: )
        a_up = 1.0D2
        a_low = 1.0D-10
      case ( 50:99 )
        a_up = 2.0D5
        a_low = 5.0D1
      case default
        write(*,*) " Temperature is out of range!"
        write(*,*) " Please set temperature higher than 50 K."
        stop
    end select

    do while ( .true. )

        call calcfp ( a_up, 32, f32_up )
        call calcfp ( a_low, 32, f32_low )
        f32_up = f32_up - temp
        f32_low = f32_low - temp
        
        if ( (f32_up * f32_low) < 0.0D0 ) then
            a_mid = 0.5D0 * ( a_up + a_low )
        else
            write(*,*) " Initial guess is unreasonable!"
            stop
        end if
    
        call calcfp ( a_mid, 32, f32_mid )
        f32_mid = f32_mid - temp
    
        if ( (f32_mid * f32_low) < 0.0D0 ) then
            a_up = a_mid
        else
            a_low = a_mid
        end if
    
        hda = 0.5D0 * ( a_up - a_low )
        
        if ( hda < error ) exit
        
    end do
    
    a = 0.5D0 * ( a_up + a_low )

end subroutine
