module element 
!
!   creates cell type to store
!   conserved values and 
!   provide update functionality
!
    use numtype
    implicit none


    type :: cell
    !
    !   finite element 
    !   
        real(dp)              :: U(5)        ! conserved quantities
        integer               :: location(3) ! (i,j,k)
        real(dp)              :: volume
        integer               :: size

    contains

        procedure :: init
        procedure :: update
        procedure :: setLocation

    end type cell


contains 

    subroutine init(this, config, vol)
    !
    !   build initialization for 
    !   compressible fluid kernel
    !
        class(cell), intent(inout) :: this
        real(dp), intent(in)       :: config(6) ! [ p, \gamma , \rho , v_x, v_y, v_z ]
        real(dp), intent(in)       :: vol
        integer                    :: size
        real(dp)                   :: e_, Energy

        ! compressible fluid kernel

        e_     = config(1) / (config(3) * (config(2) - 1)) 
        Energy = e_ + 0.5_dp * sum(config(4:6)**2)

        this%U(1) = config(3)             ! rho
        this%U(2) = config(3) * config(4) ! rho * u
        this%U(3) = config(3) * config(5) ! rho * v
        this%U(4) = config(3) * config(6) ! rho * w
        this%U(5) = config(3) * Energy    ! rho * E = rho*e

        this%size   = 5
        this%volume = vol

    end subroutine init

    subroutine update(this, dU)
    !
    !   update U values
    !
        class(cell), intent(inout) :: this
        real(dp), intent(in)       :: dU(:)

        if (size(dU) /= this%size) then
            print *, "Update vs kernal size mismatch"
            stop
        end if


        this%U = this%U + dU

    end subroutine update

    subroutine setLocation(this, coords)
    !
    !   set location of cell (i,j,k)
    !
        class(cell), intent(inout) :: this
        integer, intent(in)        :: coords(3)

        this%location = coords

    end subroutine setLocation



end module element