module cell 
!
!   creates cell type to store
!   conserved values and 
!   provide update functionality
!
    use numtype
    use, intrinsic :: ieee_arithmetic
    implicit none


    type :: cell
    !
    !   finite element 
    !   
        real(dp)              :: U(5)        ! conserved quantities
        integer               :: location(3) ! (i,j,k)
        real(dp)              :: volume, spHeat
        integer               :: size

    contains

        procedure :: init
        procedure :: update
        procedure :: setLocation
        ! functions
        procedure :: cellPressure   ! harcoded as 20-100 C value for water
                                    ! TODO: update it dynamically WRT temperature

        final     :: destructorCell__  
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
        this%spHeat = config(2)

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


    function cellPressure(this) result(pressure)
    !
    !   Returns the pressure of the cell safely
    !
        class(cell), intent(in) :: this
        real(dp)                :: pressure
        ! private
        real(dp) :: kinetic, rho, e_internal

        ! Default to NaN if things go really bad
        pressure = ieee_value(1.0_dp, ieee_quiet_nan)

        rho = this%U(1)

        if (this%volume <= 0 .or. this%size <= 0 .or. any(this%location <= 0)) then
            print *, "Misconfigured cell"
            error stop 1
        end if

        if (rho <= tiny(1.0_dp)) then
            if (isDebug .eq. 1) print *, "Invalid density in cell (rho ~ 0)"
            return
        end if

        ! Kinetic energy per volume
        kinetic = 0.5_dp * sum(this%U(2:4)**2) / rho

        ! Internal energy per volume
        e_internal = this%U(5) - kinetic

        if (e_internal <= 0.0_dp) then
            if (isDebug .eq. 1) print *, "Non-physical state: negative internal energy"
            return
        end if

        ! Pressure = (gamma - 1) * e_internal
        pressure = (this%spHeat - 1.0_dp) * e_internal

    end function cellPressure



    subroutine destructorCell__(this)
    !
    !   final for cell
    !
        type(cell), intent(inout) :: this

        this%volume = -1.0_dp
        this%size   = -1
        this%location = 0

    end subroutine destructorCell__



end module cell