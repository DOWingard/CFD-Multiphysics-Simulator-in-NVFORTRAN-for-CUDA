program cfd
!
!   test the cfd system
!
    use fluid
    implicit none

    ! instance the cfd manifold
    type(manifold) :: water

    ! parameters 
    integer  :: N, i
    real(dp) :: gamma, volume, dt, ds

    N      = 50
    dt     = 1._dp
    ds     = 1._dp
    gamma  = 7._dp
    volume = 1._dp


    ! set up the manifold
    call water%fillManifoldCPU(N,volume)
    call water%allocBundleCPU(N)

    do i = 1,3
        print*, water%fluid(i+3*i,3,4)%U(1)
    end do


    call water%initDensityGaussianCPU(N) ! adjust density to give pressure gradient

    do i = 1,3
        print*, water%fluid(i+3*i,3,4)%U(1)
    end do

    call water%loadBundleCPU(dt, ds, gamma, N)
    call water%updateCPU(dt, ds)

    do i = 1,3
        print*, water%fluid(i+3*i,3,4)%U(1)
    end do





end program cfd