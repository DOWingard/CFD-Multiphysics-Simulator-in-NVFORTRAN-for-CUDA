program cfdtest
!
!   test the cfd system
!
    use fluid
    implicit none

    ! instance the cfd manifold
    type(manifold) :: water

    ! parameters 
    integer  :: N, i
    real(dp) :: gamma, volume, dt

    N      = 50
    dt     = 0.01_dp
    gamma  = 7._dp  ! like in water config
    volume = 1._dp


    ! set up the manifold
    call water%fillManifoldCPU(N,volume)
    call water%allocBundleCPU(N)

    if (isDebug .eq. 1) then
        print*, "\nInit Values"
        do i = 1,3
            write(*,'(A,3I4,A,F12.6)') 'site: ', water%fluid(i+N/2,N/2,N/2)%location, '.  density: ', water%fluid(i+N/2,N/2,N/2)%U(1)
        end do
    end if

    call water%initDensityGaussianCPU(N) ! adjust density to give pressure gradient

    if (isDebug .eq. 1) then
        print*, "\nPost Density Init Values"
        do i = 1,3
            write(*,'(A,3I4,A,F12.6)') 'site: ', water%fluid(i+N/2,N/2,N/2)%location, '.  density: ', water%fluid(i+N/2,N/2,N/2)%U(1)
        end do
    end if

    call water%loadBundleCPU(gamma, N)
    call water%updateCPU(dt)

    if (isDebug .eq. 1) then
        print*, "\nPost time step incriment"
        do i = 1,3
            write(*,'(A,3I4,A,F12.6)') 'site: ', water%fluid(i+N/2,N/2,N/2)%location, '.  density: ', water%fluid(i+N/2,N/2,N/2)%U(1)
        end do
    end if





end program cfdtest