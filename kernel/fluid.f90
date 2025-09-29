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
        procedure :: assignLocation
        ! CPU
        procedure :: fillCPU
        procedure :: allocBundleCPU
        ! GPU

        !procedure :: updateManifold

        final     :: destructor

    end type  



contains

    subroutine initCellWater(this, cell, volume)
    !
    !   initialize the cell
    !
        class(manifold), intent(inout) :: this
        class(cell), intent(inout)     :: cell
        real(dp), intent(in)           :: volume
        ! private
        real(dp)                       :: config(6)


        config = (/ &
            1.0e5_dp,  &! (p      : pressure)
            7.00_dp,   &! (\gamma : specific heat ratio)  
            1000.0_dp, &! (\rho   : density)
            0.0_dp,    &! (       : x-velocity)
            0.0_dp,    &! (       : y-velocity)
            0.0_dp     &! (       : z-velocity)
         /)

        call cell%init(config,volume)

    end subroutine initcellwater

    subroutine assignLocation(this, location)
    !
    !   assign location to the cell in the manifold at this location
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: location(3)   

        ! call setLocation on the cell at (i,j,k)
        call this%fluid(location(1), location(2), location(3))%setLocation(location)

    end subroutine assignLocation


    subroutine fillCPU(this, N, volume)
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

        !parallel computing in fortran is truly the best

        !=========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,location)
        !=========================================================

        do i = 1,N
        do j = 1,N
        do k = 1,N

            location =  (/ i,j,k /)

            allocate(this%fluid(i,j,k))

            call assignLocation(this, location)
            call initCellWater(this, this%fluid(i,j,k),  volume)
 

        end do
        end do
        end do

        !=========================================================
        !$OMP END PARALLEL DO
        !=========================================================

        
    end subroutine fillCPU



    subroutine allocBundleCPU(this, N)
    !
    !   allocate entire tangent bundle
    !
        class(manifold), intent(inout) :: this
        integer, intent(in)            :: n

        allocate(this%tangentBundle)
        call this%tangentBundle%alloc(N)

    end subroutine



    subroutine destructor(this)
    !
    !   possible cleanup stuff
    !
        type(manifold), intent(inout) :: this

        ! IDK YET

    end subroutine destructor




! =================== Helpers =========================

function rusFlux(cellL, cellR, gamma, nhat) result(flux)
!
!   computes flux across faces using
!   Rusanov 
!
    class(cell), intent(in) :: cellL, cellR
    real(dp), intent(in)    :: gamma         ! specific heat ratio
    integer, intent(in)     :: nhat(3)       ! (+/- i^,j^,k^)
    real(dp)                :: flux(5)
    ! private            
    real(dp) :: F_L(5), F_R(5)
    real(dp) :: F_xL(5), F_yL(5), F_zL(5)
    real(dp) :: F_xR(5), F_yR(5), F_zR(5)
    real(dp) :: alpha
    real(dp) :: rhoL, uL, vL, wL, pL, EL
    real(dp) :: rhoR, uR, vR, wR, pR, ER
    integer :: i

    !  cellL 
    rhoL = cellL%U(1)
    uL   = cellL%U(2)/rhoL
    vL   = cellL%U(3)/rhoL
    wL   = cellL%U(4)/rhoL
    EL   = cellL%U(5)/rhoL
    pL   = (gamma-1.0_dp)*(rhoL*EL - 0.5_dp*rhoL*(uL**2+vL**2+wL**2))
    

    ! cellR 
    rhoR = cellR%U(1)
    uR   = cellR%U(2)/rhoR
    vR   = cellR%U(3)/rhoR
    wR   = cellR%U(4)/rhoR
    ER   = cellR%U(5)/rhoR
    pR   = (gamma-1.0_dp)*(rhoR*ER - 0.5_dp*rhoR*(uR**2+vR**2+wR**2))


    ! flux vector components
    F_xL = [rhoL*uL, rhoL*uL**2 + pL, rhoL*uL*vL, rhoL*uL*wL, uL*(rhoL*EL + pL)]
    F_yL = [rhoL*vL, rhoL*uL*vL, rhoL*vL**2 + pL, rhoL*vL*wL, vL*(rhoL*EL + pL)]
    F_zL = [rhoL*wL, rhoL*uL*wL, rhoL*vL*wL, rhoL*wL**2 + pL, wL*(rhoL*EL + pL)]

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

end function rusFlux

function totalChange(location, dt, ds, gamma, N) result(dU)
!
!   computes the total change to a cell
!   using sum of Rusanov FLuxes
!
!   assumed dx,dy,dz = ds 
    real(dp), intent(in)   :: dt, ds, gamma
    integer, intent(in)    :: N, location(3)
    real(dp)               :: dU
    ! private
    real(dp)               :: fluxes(6,5) ! 6 size(5) flux vectors
    integer                :: i, j, k



    !=========================================================
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,location)
    !=========================================================

    !do i = location(1), location(1) + 1
    !do j = location(2), location(2) + 1
    !do k = location(3), location(3) + 1

    !    location =  (/ i,j,k /)

    !    allocate(this%fluid(i,j,k))

    !    call assignLocation(this, this%fluid(i,j,k), location)
    !    call initCellWater(this, this%fluid(i,j,k),  volume)


    !end do
    !end do
    !end do

    !=========================================================
    !$OMP END PARALLEL DO
    !=========================================================



end function totalchange


! =====================================================
end module fluid