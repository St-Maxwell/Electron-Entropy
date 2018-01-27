module global
    implicit none
    integer,parameter :: inputfile = 10
    integer,parameter :: outputfile = 11
    character ( len = 4 ) :: option
    real ( kind = 8 ) :: startT, step
    real ( kind = 8 ), allocatable :: groupT (:)
    integer :: n
end module

PROGRAM electron_entropy
    use global
    implicit none
    integer :: i

    call read()

    open ( outputfile, file = "result.txt" )
    
    write(outputfile,*) "Compute the entropy of free electron gas"
    write(outputfile,*) "Programmed by St Maxwell"
    write(outputfile,*) "Reference: J. Chem. Theory Comput. 2013, 9, 3165−3169"
    write(outputfile,*) 
    write(outputfile,*) "T(K)     S(J/mol-K)  dT(K)"

    do i = 1, n
        call calcS ( groupT ( i ) )
    end do
    
    stop
END PROGRAM
    
subroutine calcS ( T0 )
    use global
    implicit none
    real ( kind = 8 ), parameter :: h = 6.626070040D0 ! D-34
    real ( kind = 8 ), parameter :: kB = 1.38064852D0 ! D-23
    real ( kind = 8 ), parameter :: mass_e = 9.10938356D0 ! D-31
    real ( kind = 8 ), parameter :: pi = 3.14159265358979D0
    real ( kind = 8 ), parameter :: J = 2.0D0 ! degree of degeneration
    real ( kind = 8 ), parameter :: atm = 1.0D0 ! D5
    real ( kind = 8 ), parameter :: R = 8.314D0
    real ( kind = 8 ) :: fp32, fp12, KK, T0, T, Sm, dT
    real ( kind = 8 ) :: a 
    
    KK = ( atm * h**3.0D0 * 1.0D7 ) / ( J * ( 2.0D0 * pi * mass_e )**1.5D0 * kB**2.5D0 )
    call bisection ( T0, KK, a )

    call calcfp( a, 32, fp32 )
    call calcfp( a, 12, fp12 )
    
    T = ( KK / fp32 )**0.4D0
    Sm = R * ( ( 5.0D0 * fp32 ) / ( 2.0D0 * fp12 ) - log( a ) )
    dT = abs ( T - T0 )

    write(outputfile,"(F7.2,'   ',F8.4,'   ',F12.10)") T, Sm, dT

end subroutine

subroutine calcfp ( a, p, res )
    implicit none
    external fp
    real ( kind = 8 ) :: low
    real ( kind = 8 ) :: up
    real ( kind = 8 ) :: tol
    real ( kind = 8 ) :: a
    integer :: p
    real ( kind = 8 ) :: res
    real ( kind = 8 ), parameter :: pi = 3.14159265358979D0
    real ( kind = 8 ), parameter :: const12 = 2.0D0 / sqrt ( pi )
    real ( kind = 8 ), parameter :: const32 = 4.0D0 / ( 3.0D0 * sqrt ( pi ) )
    
    if ( a < 1.0D21 ) then
        low = 0.0D0
        up = 1.0D2
        tol = 1.0D-9
    else if ( a < 1.0D42 ) then
        low = 0.0D0
        up = 1.5D2
        tol = 1.0D-9
    else if ( a < 1.0D104 ) then
        low = 0.0D0
        up = 3.0D2
        tol = 1.0D-8
    else if ( a < 1.0D160 ) then
        low = 0.0D0
        up = 4.0D2
        tol = 1.0D-7
    else
        low = 0.0D0
        up = 5.0D2
        tol = 1.0D-7
    end if
  
    call auto_gauss ( fp, res, low, up, tol, a, p )
    
    select case ( p )
        case ( 12 )
            res = res * const12
        case ( 32 )
            res = res * const32
        case default
            stop
    end select
    
end subroutine

subroutine fp ( f, x, a, p )
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

subroutine auto_gauss ( func, s, low, up, tol, a, p )
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-31
!-----------------------------------------------------
!  Purpose   :  变步长高斯积分方法
!-----------------------------------------------------
!  Input  parameters  :
!       1.  func 外部函数（待计算的函数）
!       2.  s 积分结果
!       3.  a,b 积分区间
!       4.  tol  误差容限
!  Output parameters  :
!       1.  s  积分结果
!       2.  M  实际区间划分个数
!  Common parameters  :
!       1.    a     fp函数参数
!       2.    p     fp函数参数
!----------------------------------------------------
!  Post Script :
!       1.  需要调用复合高斯积分函数
!----------------------------------------------------
    implicit none
    real ( kind = 8 ) :: s, low, up, tol, s1, s2, del, a
    external :: func
    integer :: m, i, p
    
    m = 2

    !最大允许重新划分20次，
    do i = 1, 60
        
      call com_Gauss ( func, s1, low, up, m, a, p )
      
      !划分细度加倍
      m = m * 2
      call com_Gauss ( func, s2, low, up, m, a, p )  
      
      !前后两次积分值之差
      del = abs ( s2 - s1 )
      
      if ( del < tol ) exit
      
    end do
    
    s = s2
    
end subroutine

subroutine com_Gauss ( func, s, low, up, n, a, p )
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  复合5点高斯勒让德积分函数
!-----------------------------------------------------
!  Input  parameters  :
!       1. func 外部函数
!       2. a,b 积分区间
!       3. n 区间划分个数
!  Output parameters  :
!       1.  s 积分结果
!  Common parameters  :
!       1.    a     fp 函数参数
!       2.    p     fp 函数参数
!----------------------------------------------------
!  Post Script :
!       1.   需要调用单区间高斯勒让德公式函数
!----------------------------------------------------
    implicit none
    real ( kind = 8 ) :: s, s1, low, up, hstep, c, d, a
    external :: func
    integer :: n, i, p
    
    hstep = ( up - low ) / n
    s = 0.0D0
    
    do i = 1, n
        c = low + ( i - 1 ) * hstep
        d = low + i * hstep
        
        call GL ( func, s1, c, d, a, p )
        
        s = s + s1
    end do
    
end subroutine

subroutine GL ( func, s, low, up, a, p )
!---------------subroutine  comment-------------------
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :5点 Gauss-Legendre积分
!-----------------------------------------------------
!  Input  parameters  :
!       1.    func  外部函数
!       2.    a,b   积分区间
!  Output parameters  :
!       1.    s 积分结果
!  Common parameters  :
!       1.    a     fp函数参数
!       2.    p     fp函数参数
!----------------------------------------------------
!  Post Script :
!       1.  选用5节点的公式
!       2.  计算时先把一般区间变换到[-1,1]区间上利用公式
!----------------------------------------------------
    implicit none
    real ( kind = 8 ) :: s, fx, low, up, a
    external :: func
    integer :: i, p
    real ( kind = 8 ) :: node ( 5 ), w ( 5 ), t ( 5 )
    
    node = (/ -0.9061798459, -0.5384693101, 0.0000000000, &
               0.5384693101, 0.9061798459 /)
    w = (/ 0.2369268851, 0.4786286705, 0.5688888889, &
           0.4786286705, 0.2369268851 /)
    
    t = ( low + up ) / 2.0D0 + ( up - low ) * node / 2.0D0  
    
    s = 0.0D0
    
    do i = 1, 5
        call func ( fx, t(i), a, p )
        s = s + fx * w ( i )
    end do
    
    s = s * ( up - low ) / 2.0D0
    
end subroutine

subroutine bisection( T, KK, a )
    implicit none
    real ( kind = 8 ), intent ( in ) :: T, KK
    real ( kind = 8 ) :: temp
    real ( kind = 8 ) :: a, a_up, a_low, a_mid, hda
    real ( kind = 8 ) :: f32_up, f32_low, f32_mid
    real ( kind = 8 ) :: error = 1.0D-10
    integer :: counter, eff_num

    temp = KK / T**2.5D0
    
    if ( T < 5.0D1 ) then
        eff_num = ceiling ( 1.22164D0 + 193.43663D0 / ( T - 0.06002D0 ) ) - 10
        error = 10**dble(eff_num)
    end if
    
    select case ( ceiling( T ) )
        case ( 101: )
            a_up = 1.0D3
            a_low = 0.0D0
        case ( 51:100 )
            a_up = 10.0D0 ** ( -0.07D0 * T + 9.0D0 )
            a_low = 10.0D0 ** ( -0.07D0 * T + 7.4D0 )
        case ( 21:50 )
            a_up = 10.0D0 ** ( -0.15D0 * T + 14.0D0 )
            a_low = 10.0D0 ** ( -0.15D0 * T + 11.0D0 )  
        case ( 11:20 )
            a_up = 10.0D0 ** ( -1.0D0 * T + 31.0D0 )
            a_low = 10.0D0 ** ( -1.0D0 * T + 26.0D0 )
        case ( 6:10 )
            a_up = 10.0D0 ** ( -4.2D0 * T + 64.0D0 )
            a_low = 10.0D0 ** ( -4.2D0 * T + 53.0D0 )
        case ( 3:5 )
            a_up = 10.0D0 ** ( -21.0D0 * T + 149.0D0 )
            a_low = 10.0D0 ** ( -21.0D0 * T + 128.0D0 )
        case ( 1:2 )
            a_up = 10.0D0 ** ( -100.0D0 * T + 315.0D0 )
            a_low = 10.0D0 ** ( -100.0D0 * T + 275.0D0 )
        case default
            write(*,*) " Temperature is out of range!"
            read(*,*)
            stop
    end select
        
    counter = 0
    
    do while ( .true. )
        
        call calcfp ( a_up, 32, f32_up )
        call calcfp ( a_low, 32, f32_low )
        
        f32_up = f32_up - temp
        f32_low = f32_low - temp
        
        if ( (f32_up * f32_low) < 0.0D0 ) then
            a_mid = 0.5D0 * ( a_up + a_low )
        else if ( (f32_up * f32_low) == 0.0D0 ) then
            exit
        else
            write(*,*) " Unreasonable guess!"
            read(*,*)
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
        counter = counter + 1
        
        if ( hda < error .or. counter > 128 ) then
            if ( counter == 129 ) write(*,"(' Bisection step exceeds')") 
            a = 0.5D0 * ( a_up + a_low )
            exit
        end if
        
    end do
    
    if ( f32_up == 0.0D0 ) then
        a = a_up
    else if ( f32_up == 0.0D0 ) then
        a = a_low
    end if
    
end subroutine

subroutine read ()
    use global
    implicit none
    integer :: i
    logical :: alive
    integer :: error

    inquire ( file = "eS_input.txt", exist = alive )
    if ( .not. alive ) then
        write(*,*) "eS_input.txt dosen't exist."
        read(*,*)
        stop
    end if

    open ( inputfile, file = "eS_input.txt", status = "old", iostat = error )
    if ( error /= 0 ) then
        write(*,*) "Fail to open eS_input.txt."
        read(*,*)
        stop
    end if

    read(inputfile,*) option

    if ( option == "scan" ) then
        read(inputfile,*) startT, step, n
        allocate(groupT(n))
        do i = 1, n
            groupT ( i ) = startT + step * dble( i - 1 )
        end do
    else if ( option == "1by1" ) then
        read(inputfile,*) n
        read(inputfile,*) 
        allocate(groupT(n))
        do i = 1, n
            read(inputfile,*) groupT ( i )
        end do
    else
        write(*,*) "Invalid option!"
        read(*,*)
    end if

end subroutine
