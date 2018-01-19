PROGRAM electron_entropy
    implicit none
    real ( kind = 8 ) :: T
    real ( kind = 8 ) :: T0 = 250.0D0
    real ( kind = 8 ) :: step = 1.0D0
    integer :: i
    write(*,*) "T(K)     S(J/mol·K)  A"
    
    do i = 1,100
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

    write(*,"(F7.2,'   ',F7.4,'   ',F10.6)") T, Sm, a

end subroutine

subroutine calcfp( a, p, res )
  implicit none
  integer,parameter :: n = 40000
  real ( kind = 8 ) , external :: fp
  real ( kind = 8 ) , parameter :: low = 0.0D0
  real ( kind = 8 ) , parameter :: up = 1.2D2
  real ( kind = 8 ) :: a
  integer :: p
  real ( kind = 8 ) :: res

  call integral ( fp, low, up, res, n, a, p )

  select case ( p )
    case( 12 )
      res = res * 1.128379167D0
    case( 32 )
      res = res * 0.752252778D0
    case default
      stop
  end select

end subroutine

function fp ( x, a, p )
  implicit none
  real ( kind = 8 ) :: fp, x, a
  integer :: p
  real ( kind = 8 ), parameter :: pi = 3.14159265358979D0

  select case ( p )
    case( 12 )
      fp = x**0.5D0 / ( 1.0D0 + exp ( x ) / a )
    case( 32 )
      fp = x**1.5D0 / ( 1.0D0 + exp ( x ) / a )
    case default
      stop
  end select

  return
end function

subroutine integral ( fp, low, up, res, n, a, p )
  implicit none
  real ( kind = 8 ), external :: fp
  real ( kind = 8 ) :: low, up, res, a, l, h
  real ( kind = 8 ) :: temp
  integer :: i, n, p

  l = up - low
  h = l / dble( n )
    
  temp = 0.0D0
  do i=1,n-1
    temp = temp + 2.0D0 * fp( low + dble(i) * h, a, p )
  end do

  res = 0.5D0 * h * ( fp( low, a, p ) + fp( up, a, p ) + temp )

end subroutine

subroutine bisection( T, KK, a )
    implicit none
    real ( kind = 8 ), intent ( in ) :: T, KK
    real ( kind = 8 ) :: temp
    real ( kind = 8 ) :: a, a_up, a_low, a_mid, hda
    real ( kind = 8 ) :: f32_up, f32_low, f32_mid
    real ( kind = 8 ), parameter :: error = 1.0D-7

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
        
        if ( hda < error ) then
            a = 0.5D0 * ( a_up + a_low )
            exit
        else
            cycle
        end if
        
    end do

end subroutine
