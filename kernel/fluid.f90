module fluid
!
!   module to house type that contains cells 
!   and related logic
!
    use cell
    use tangent
    implicit none

    type :: manifold

        class(cell),   allocatable :: fluid(:,:,:)
        class(bundle), allocatable :: tangentBundle
        integer                    :: coords(3), N


    contains

        procedure :: initCellWater
        procedure :: allocFluid
        procedure :: assignLocation
        ! CPU
        procedure :: fillManifoldCPU
        procedure :: allocBundleCPU
        procedure :: loadBundleCPU
        procedure :: initDensityGaussianCPU
        procedure :: updateCPU

        ! CPU TODO
        !procedure :: computeAllBasis  ! instructs cells to compute basis

        ! GPU

        !procedure :: updateManifold

        final     :: destructorManifold__

    end type  



contains

    subroutine initCellWater(this, location, volume)
    !
    !   initialize the cell
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: location(3)
        real(dp), intent(in)           :: volume
        ! private
        real(dp)                       :: config(6)


        config = (/ &
            1.0e5_dp,  &! (p      : pressure)
            4.186_dp,  &! (\gamma : specific heat ratio)  
            1000.0_dp, &! (\rho   : density)
            0.0_dp,    &! (       : x-velocity)
            0.0_dp,    &! (       : y-velocity)
            0.0_dp     &! (       : z-velocity)
         /)

        call this%fluid(location(1),location(2),location(3))%init(config, volume)

    end subroutine initcellwater


    subroutine allocFluid(this, N)
    !
    !   allocate the Fluid    
    !
        class(manifold), intent(inout) :: this 
        integer, intent(in)            :: N



        if (allocated(this%fluid)) deallocate(this%fluid)
        allocate(this%fluid(N,N,N))

    end subroutine

    subroutine assignLocation(this, location)
    !
    !   assign location to the cell in the manifold at this location
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: location(3)   

        ! call setLocation on the cell at (i,j,k)
        call this%fluid(location(1), location(2), location(3))%setLocation(location)

    end subroutine assignLocation


    subroutine fillManifoldCPU(this, N, volume)
    !
    !   populates manifold with cells
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: N
        real(dp), intent(in)           :: volume
        ! private
        integer                        :: i, j, k, location(3)
        
        if (allocated(this%fluid)) deallocate(this%fluid)
        allocate(this%fluid(N,N,N))

        this%N = N

        ! allocate separately solved memory bug
        do i = 1,N
        do j = 1,N
        do k = 1,N

            location =  (/ i,j,k /)
            allocate(this%fluid(i,j,k))

        end do
        end do
        end do


        !=========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,location)
        !=========================================================

        do i = 1,N
        do j = 1,N
        do k = 1,N

            location =  (/ i,j,k /)

            call assignLocation(this, location)
            call initCellWater(this, location,  volume)
 

        end do
        end do
        end do

        !=========================================================
        !$OMP END PARALLEL DO
        !=========================================================

        
    end subroutine fillManifoldCPU



    subroutine allocBundleCPU(this, N)
    !
    !   allocate entire tangent bundle
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: n


        if (allocated(this%tangentBundle)) then 
            deallocate(this%tangentBundle)
        end if

        allocate(this%tangentBundle)
        call this%tangentBundle%alloc(N, this%fluid(1,1,1)%size)

    end subroutine


    subroutine loadBundleCPU(this, gamma, N)
    !
    !   loads new values into the tangent bundle
    !   and smooths them (roughly) to maintain 
    !   numerical stability
    !
        class(manifold), intent(inout) :: this 
        real(dp), intent(in)           :: gamma
        integer, intent(in)            :: N
        ! private
        integer :: i, j, k, f
        integer :: ni, nj, nk
        integer :: nhat(3)
        real(dp) :: flux(5), oldFlux(5)
        real(dp) :: threshold
        class(cell), pointer :: cellL, cellR
        integer :: component

        threshold = 10.0_dp   ! spike threshold for smoothing

        if (.not. allocated(this%tangentBundle)) then
            print*, "tangent Bundle not allocated"
            return
        end if

        !=====================================================================================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,f,ni,nj,nk,nhat,flux,oldFlux,cellL,cellR,component)
        !=====================================================================================================================
        do i = 1, N
            do j = 1, N
                do k = 1, N

                    cellR => this%fluid(i,j,k)

                    do f = 1,6
                        ! set neighbors
                        select case(f)
                        case(1)
                            ni = max(i-1,1); nj = j; nk = k; nhat = (/-1,0,0/)
                        case(2)
                            ni = min(i+1,N); nj = j; nk = k; nhat = (/1,0,0/)
                        case(3)
                            ni = i; nj = max(j-1,1); nk = k; nhat = (/0,-1,0/)
                        case(4)
                            ni = i; nj = min(j+1,N); nk = k; nhat = (/0,1,0/)
                        case(5)
                            ni = i; nj = j; nk = max(k-1,1); nhat = (/0,0,-1/)
                        case(6)
                            ni = i; nj = j; nk = min(k+1,N); nhat = (/0,0,1/)
                        end select

                        ! free-flow boundary
                        if (ni < 1 .or. ni > N .or. nj < 1 .or. nj > N .or. nk < 1 .or. nk > N) then
                            cellL => cellR
                        else
                            cellL => this%fluid(ni,nj,nk)
                        end if

                        flux = rusFlux(cellL, cellR, gamma, nhat)

                        oldFlux = this%tangentBundle%mesh(i,j,k)%fluxes(f,:)

                        ! smooth unrealistic jumps and handle NaN/Inf
                        do component = 1,5

                            if (ieee_is_nan(flux(component))) then

                                flux(component) = 1.5*oldFlux(component)

                            else if (.not. ieee_is_finite(flux(component))) then

                                flux(component) = sign(1._dp, flux(component)) * max(abs(oldFlux(component))*10._dp, 1._dp)

                            else if (abs(flux(component) - oldFlux(component)) > threshold*max(abs(oldFlux(component)),1._dp)) then

                                flux(component) = oldFlux(component) * 1.1_dp * sign(flux(component) - oldFlux(component), 1._dp)

                            end if
                        end do

                        this%tangentBundle%mesh(i,j,k)%fluxes(f,:) = flux

                    end do

                end do
            end do
        end do
        !=====================================================================================================================
        !$OMP END PARALLEL DO
        !=====================================================================================================================
    end subroutine loadBundleCPU



    subroutine initDensityGaussianCPU(this, N)
    !
    !   components the densities by scaling mesh
    !
        class(manifold), intent(inout)    :: this
        integer, intent(in)               :: N
        ! private
        real(dp), allocatable             :: scalingMesh(:,:,:)
        integer                           :: i, j, k


        scalingMesh = pressureAtCenter(N)


        !=========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k)
        !=========================================================

        do i = 1,N
        do j = 1,N
        do k = 1,N

            this%fluid(i,j,k)%U = this%fluid(i,j,k)%U(:) * scalingMesh(i,j,k)

        end do
        end do
        end do

        !=========================================================
        !$OMP END PARALLEL DO
        !=========================================================

    end subroutine initDensityGaussianCPU



    subroutine updateCPU(this, dt)
    !
    !   Updates all cells in the manifold using
    !   the 6 stored fluxes in the tangent bundle.
    !
        class(manifold), intent(inout) :: this
        real(dp), intent(in)           :: dt
        integer                        :: i,j,k,f
        real(dp)                       :: dU(5)
        integer                        :: N

        N = this%N

        !=========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,f,dU)
        !=========================================================
        do i = 1,N
        do j = 1,N
        do k = 1,N

            dU = 0.0_dp

            do f = 1,6
                dU = dU + this%tangentBundle%mesh(i,j,k)%fluxes(f,:)  
            end do

            dU = - ( dt / (this%fluid(i,j,k)%volume )) * dU
            !print*, dU

            ! Update cell's conserved variables
            this%fluid(i,j,k)%U = this%fluid(i,j,k)%U + dU   
 
        end do
        end do
        end do
        !=========================================================
        !$OMP END PARALLEL DO
        !=========================================================

    end subroutine updatecpu

    !subroutine computeAllBasis(this)
    !
    !   computes the 3d basis functions
    !   of each cell
    !
        !type(manifold), intent(inout) :: this 
        
        ! for all i,j,k call cell%computeBasis

    !end subroutine computeAllBasis


    subroutine destructorManifold__(this)
        !
        !   Finalizer for manifold: safely deallocate all internal memory
        !
        type(manifold), intent(inout) :: this
        integer :: i,j,k


        if (allocated(this%fluid)) then
            deallocate(this%fluid)
        end if

        if (allocated(this%tangentBundle)) then
            if (allocated(this%tangentBundle%mesh)) then
                do k = 1, size(this%tangentBundle%mesh,3)
                do j = 1, size(this%tangentBundle%mesh,2)
                do i = 1, size(this%tangentBundle%mesh,1)
                    if (allocated(this%tangentBundle%mesh(i,j,k)%fluxes)) then
                        deallocate(this%tangentBundle%mesh(i,j,k)%fluxes)
                    end if
                end do
                end do
                end do
                deallocate(this%tangentBundle%mesh)
            end if
            deallocate(this%tangentBundle)
        end if


        this%N      = -1
        this%coords =  0

    end subroutine destructorManifold__




! =================== Helpers =========================

    function rusFlux(cellL, cellR, gamma, nhat) result(flux)
    !
    !   very loose stability management, could use a lot more thought
    !   with physicial justifications
    !
        class(cell), intent(in) :: cellL, cellR
        real(dp), intent(in)    :: gamma
        integer, intent(in)     :: nhat(3)
        real(dp)                :: flux(5)
        
        real(dp) :: F_L(5), F_R(5)
        real(dp) :: F_xL(5), F_yL(5), F_zL(5)
        real(dp) :: F_xR(5), F_yR(5), F_zR(5)
        real(dp) :: alpha
        real(dp) :: rhoL, uL, vL, wL, pL, EL
        real(dp) :: rhoR, uR, vR, wR, pR, ER
        integer :: i

        rhoL = max(cellL%U(1), 1e-12_dp)
        uL   = cellL%U(2)/rhoL
        vL   = cellL%U(3)/rhoL
        wL   = cellL%U(4)/rhoL
        EL   = cellL%U(5)/rhoL
        pL   = max((gamma-1.0_dp)*(rhoL*EL - 0.5_dp*rhoL*(uL**2+vL**2+wL**2)), 1e-12_dp)

        F_xL = [rhoL*uL, rhoL*uL**2 + pL, rhoL*uL*vL, rhoL*uL*wL, uL*(rhoL*EL + pL)]
        F_yL = [rhoL*vL, rhoL*uL*vL, rhoL*vL**2 + pL, rhoL*vL*wL, vL*(rhoL*EL + pL)]
        F_zL = [rhoL*wL, rhoL*uL*wL, rhoL*vL*wL, rhoL*wL**2 + pL, wL*(rhoL*EL + pL)]

        rhoR = max(cellR%U(1), 1e-12_dp)
        uR   = cellR%U(2)/rhoR
        vR   = cellR%U(3)/rhoR
        wR   = cellR%U(4)/rhoR
        ER   = cellR%U(5)/rhoR
        pR   = max((gamma-1.0_dp)*(rhoR*ER - 0.5_dp*rhoR*(uR**2+vR**2+wR**2)), 1e-12_dp)

        F_xR = [rhoR*uR, rhoR*uR**2 + pR, rhoR*uR*vR, rhoR*uR*wR, uR*(rhoR*ER + pR)]
        F_yR = [rhoR*vR, rhoR*uR*vR, rhoR*vR**2 + pR, rhoR*vR*wR, vR*(rhoR*ER + pR)]
        F_zR = [rhoR*wR, rhoR*uR*wR, rhoR*vR*wR, rhoR*wR**2 + pR, wR*(rhoR*ER + pR)]

        F_L = nhat(1)*F_xL + nhat(2)*F_yL + nhat(3)*F_zL
        F_R = nhat(1)*F_xR + nhat(2)*F_yR + nhat(3)*F_zR

        alpha = max(abs(uL*nhat(1)+vL*nhat(2)+wL*nhat(3)) + sqrt(gamma*pL/rhoL), &
                    abs(uR*nhat(1)+vR*nhat(2)+wR*nhat(3)) + sqrt(gamma*pR/rhoR))

        do i = 1,5
            flux(i) = 0.5_dp*(F_L(i) + F_R(i)) - 0.5_dp*alpha*(cellR%U(i) - cellL%U(i))
        end do

        flux(1) = max(flux(1), 1e-12_dp)  
        flux(5) = max(flux(5), 1e-12_dp)  

    end function rusFlux




    ! create demo 3d Point Spread Function 
    ! type pressure scaling mesh

    function pressureAtCenter(N) result(p)
    !
    !   gaussian pressure bulge at center
    !
        implicit none
        integer, intent(in) :: N
        real(dp), allocatable :: p(:,:,:)
        real(dp) :: component, sigma, r2
        integer :: i, j, k, ic, jc, kc

        ic = N/2
        jc = N/2
        kc = N/2

        component = 0.1_dp       ! small bump 
        sigma = real(N, dp)/6

        allocate(p(N,N,N))

        !==========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,r2) SHARED(p)
        !==========================================================
        do i = 1, N
            do j = 1, N
                do k = 1, N
                    r2 = real((i-ic)**2 + (j-jc)**2 + (k-kc)**2, dp)
                    p(i,j,k) = max(1.0_dp + component*exp(-r2/(2.0_dp*sigma**2)), 1e-12_dp)
                end do
            end do
        end do
        !==========================================================
        !$OMP END PARALLEL DO
        !==========================================================
    end function pressureAtCenter






end module fluid