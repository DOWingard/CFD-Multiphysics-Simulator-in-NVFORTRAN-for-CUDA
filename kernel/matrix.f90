
module matrix

    use NumType, only : dp, nb
    implicit none
    
    contains
        
        subroutine dsydet(n,matr,lda,det)
        
            integer :: n, lda
            real(dp) :: matr(lda,lda), det(2)
            real(dp), dimension(nb*lda) :: work
            integer :: info, ipvt(lda), k
            real(dp) :: t, d
        
            info = 0
            call dsytrf('u',n,matr,lda,ipvt,work,nb*lda,info)
        
            det(1) = 1.0_dp
            det(2) = 0.0_dp
            t = 0.0_dp
            do k = 1, n
                d = matr(k,k)
                if ( ipvt(k) <= 0 ) then
                    if ( t == 0.0_dp ) then
                        t = matr(k,k+1)
                        d = (d/t)*matr(k+1,k+1) - t
                    else
                        d = t
                        t = 0.0_dp
                    end if
                end if
                det(1) = d*det(1)
                if ( det(1) /= 0.0_dp ) then
                    do while ( abs(det(1)) <= 1 )
                        det(1) = 10*det(1)
                        det(2) = det(2) - 1
                    end do
                    do while ( abs(det(1)) >= 10 )
                        det(1) = det(1)/10
                        det(2) = det(2) + 1
                    end do
                end if
            end do
            
        end subroutine dsydet
        
        subroutine zsydet(n,matr,lda,det)
        
            integer :: n, lda
            complex(dp) :: matr(lda,lda), det(2)
            complex(dp), dimension(nb*lda) :: work
            integer :: info, ipvt(lda), k
            complex(dp) :: t, d
        
            info = 0
            call zsytrf('u',n,matr,lda,ipvt,work,nb*lda,info)
        
            det(1) = 1.0_dp
            det(2) = 0.0_dp
            t = 0.0_dp
            do k = 1, n
                d = matr(k,k)
                if ( ipvt(k) <= 0 ) then
                    if ( t == 0.0_dp ) then
                        t = matr(k,k+1)
                        d = (d/t)*matr(k+1,k+1) - t
                    else
                        d = t
                        t = 0.0_dp
                    end if
                end if
                det(1) = d*det(1)
                if ( det(1) /= 0.0_dp ) then
                    do while ( abs(det(1)) <= 1 )
                        det(1) = 10*det(1)
                        det(2) = det(2) - 1
                    end do
                    do while ( abs(det(1)) >= 10 )
                        det(1) = det(1)/10
                        det(2) = det(2) + 1
                    end do
                end if
            end do
            
        end subroutine zsydet
        
        subroutine dgei(n,matr,lda)
        
            integer, intent(in) :: n, lda
            real(dp) :: matr(lda,lda)
            real(dp), dimension(nb*lda) :: work
            integer :: info, ipvt(lda)
            
            info=0
            call dgetrf(n,n,matr,lda,ipvt,info)
                if(info /= 0) stop ' dgetrf '
            call dgetri(n,matr,lda,ipvt,work,nb*lda,info)
                if(info /= 0) stop ' dgetri '
        
        end subroutine dgei
        
        subroutine zgei(n,matr,lda)
        
            integer, intent(in) :: n, lda
            complex(dp) :: matr(lda,lda)
            complex(dp), dimension(nb*lda) :: work
            integer :: info, ipvt(lda)
            
            info=0
            call zgetrf(n,n,matr,lda,ipvt,info)
                if(info /= 0) stop ' dgetrf '
            call zgetri(n,matr,lda,ipvt,work,nb*lda,info)
                if(info /= 0) stop ' dgetri '
        
        end subroutine zgei
        
        subroutine dsyi(n,matr,lda)
        
            integer, intent(in) :: n, lda
            real(dp) :: matr(lda,lda)
            real(dp), dimension(nb*lda) :: work
            integer :: i, j, info, ipvt(lda)
            
            info=0
            call dsytrf('u',n,matr,lda,ipvt,work,nb*lda,info)
                if(info /= 0) stop ' dgetrf '
            call dsytri('u',n,matr,lda,ipvt,work,info)
                if(info /= 0) stop ' dgetri '
            forall ( i=1:n, j=1:n, i < j ) matr(j,i) = matr(i,j)
        
        end subroutine dsyi
        
        
        subroutine zsyi(n,matr,lda)
        
            integer, intent(in) :: n, lda
            complex(dp) :: matr(lda,lda)
            complex(dp), dimension(nb*lda) :: work
            integer :: i, j, info, ipvt(lda)
            
            info=0
            call zsytrf('u',n,matr,lda,ipvt,work,nb*lda,info)
                if(info /= 0) stop ' dgetrf '
            call zsytri('u',n,matr,lda,ipvt,work,info)
                if(info /= 0) stop ' dgetri '
            forall ( i=1:n, j=1:n, i < j ) matr(j,i) = matr(i,j)
        
        end subroutine zsyi
        
end module matrix
