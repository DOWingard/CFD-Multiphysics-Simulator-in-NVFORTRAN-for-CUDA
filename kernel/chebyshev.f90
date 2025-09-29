

module chebyshev

    use numtype
    implicit none
    integer, parameter :: maxch = 50
    real(dp), dimension(0:maxch) :: cheb, chder, chder2
    real(dp), dimension(maxch) :: z0
    integer :: iz0
	
    contains
        
        subroutine chebyex(func,n,a,ya,yb) 
        !   func([ya,yb]) = sum_{i=0}^n  a_i T_i
    
            real(dp), external :: func
            integer :: n
            real(dp), dimension(0:maxch) :: f, a
            real(dp) :: ya, yb, aa, bb, x, ss
            integer :: i, j
    
            if ( n > maxch ) stop '  n > maxch '
            aa = (yb-ya)/2; bb = (yb+ya)/2
            do i = 0, n
                x = cos(pi/(n+1)*(i+0.5_dp))
                f(i) = func(aa*x+bb)
            end do
            do j = 0, n
                ss = 0._dp
                do i = 0, n
                    ss = ss + &
                        f(i)*cos((pi/(n+1))*j*(i+0.5_dp))
                end do
                a(j) = 2._dp*ss/(n+1)
            end do
            a(0) = 0.5_dp*a(0)
            
        end subroutine chebyex
      
        subroutine chebyderiv(n,a,der,ya,yb) ! 

            integer :: n
            real(dp) :: ya, yb, a(0:maxch), der(0:maxch)
            integer :: j
    
            der(n) = 0._dp; der(n-1) = 2*n*a(n)
            do j = n-1, 1, -1
                der(j-1) = der(j+1)+2*j*a(j)
            end do
            der(0) = der(0)/2
            der(0:n-1) = der(0:n-1)*2/(yb-ya)

        end subroutine chebyderiv
      
        function cheby(y,n,a,ya,yb) result(t) 
        ! func(y) =  sum_{i=0}^n  a_i T_i (x)

            implicit none
            integer :: n
            real(dp) :: y, ya, yb
            real(dp) :: a(0:maxch)
            real(dp) :: aa, bb, x, t, y0, y1
            integer :: k
    
            aa = (yb-ya)/2; bb = (yb+ya)/2
            x = (y-bb)/aa
            y1 = 0._dp; y0 = a(n)
            do k = n-1, 0, -1
                t = y1; y1 = y0
                y0 = a(k)+2*x*y1-t
            end do
            t = y0-x*y1

        end function cheby
                
        subroutine chebyzero(n,a,ya,yb,z0,iz0) 
        !   find zero by using Boyd's method
        
            integer :: n, iz0
            real(dp), dimension(0:maxch) :: a
            integer :: j
            real(dp), dimension(maxch) :: wr0, wi0, z0, wwr0
            real(dp) :: ya, yb
            
            call boyd(n,a,wr0,wi0)
            wwr0(1:n) = wr0(1:n)*(yb-ya)/2+(yb+ya)/2
            
            iz0 = 0
            do j = 1, n
                if( wi0(j) == 0._dp .and. &
                  -1 <= wr0(j) .and. wr0(j) <= 1 ) then 
                        iz0 = iz0+1;  z0(iz0) = wwr0(j)
                end if
            end do
                        
            contains
                
                subroutine boyd(n,a,wr,wi)
                
                    integer :: n, j, ie
                    real(dp) :: a(0:maxch)
                    real(dp) :: wr(maxch), wi(maxch)
                    integer, parameter :: lwork=4*maxch
                    real(dp) :: aamat(maxch,maxch),  &
                        work(lwork), rwork(lwork), &
                        vl(1), vr(1)
                    
                    if (abs(a(n)) == 0._dp) stop 'a(n)=0'
                    aamat(1:n,1:n) = 0._dp
                    aamat(1,2) = 1._dp
                    do j = 2, n-1
                        aamat(j,j-1) = 0.5_dp
                        aamat(j,j+1) = 0.5_dp
                    end do
                    aamat(n,1:n) = -a(0:n-1)/(2*a(n))
                    aamat(n,n-1) = aamat(n,n-1) + 0.5_dp
                    
                    ie = 0
                    call dgeev('n','n',n,aamat,maxch,wr,&
                        wi,vl,1,vr,1,work,lwork,rwork,ie)                      
                    if( ie /= 0 ) stop ' boyd: ie /= 0 ' 
                            
                end subroutine boyd
                
        end subroutine chebyzero
        
        subroutine root_polish(func,zz,dz,eps,maxf)
        
            complex(dp), external :: func
            real(dp) :: eps
            complex(dp) :: zz, dz, z1, z2, z3, &
                f1, f2, f3, a12, a23, a31
            integer :: i, maxf
            
            z1 = zz+dz;   f1 = func(z1)
            z2 = zz-dz;   f2 = func(z2)
            z3 = zz;      f3 = func(z3)

            do i = 1,maxf
                a23 = (z2-z3)*f2*f3
                a31 = (z3-z1)*f1*f3
                a12 = (z1-z2)*f1*f2
                zz = (z1*a23+z2*a31+z3*a12)/(a23+a31+a12)
                if ( abs(zz-z3) < eps ) exit
                z1 = z2;  f1 = f2
                z2 = z3;  f2 = f3
                z3 = zz;  f3 = func(z3)
            end do
        
        end subroutine root_polish
      
end module chebyshev


