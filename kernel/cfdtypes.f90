module cfdtypes
!
!   provides possible types for CFD sims
!   scalable for Plasma, etc.
!
    use numtype
    implicit none

    type :: cfdtype
        character(len=20), allocatable :: fluids(:)  ! dynamic array of fluid names
    contains
        procedure :: addFluid
        procedure :: init

        final     :: destructorTypes__
    end type cfdtype

contains

    subroutine init(this)
        class(cfdtype), intent(inout) :: this
        ! initialize with one fluid type: "water"
        allocate(this%fluids(1))
        this%fluids(1) = "water"
    end subroutine init

    subroutine addFluid(this, name)
        class(cfdtype), intent(inout) :: this
        character(len=*), intent(in)  :: name
        integer :: n
        character(len=20), allocatable :: tmp(:)

        if (.not. allocated(this%fluids)) then
            allocate(this%fluids(1))
            this%fluids(1) = name
            return
        end if

        n = size(this%fluids)

        allocate(tmp(n+1))
        tmp(1:n) = this%fluids
        tmp(n+1) = name

        deallocate(this%fluids)
        call move_alloc(tmp, this%fluids)

    end subroutine addFluid



    subroutine destructorTypes__(this)
    !
    !   Finalizer for cfdtype: deallocate dynamic arrays
    !
        type(cfdtype), intent(inout) :: this

        if (allocated(this%fluids)) then
            deallocate(this%fluids)
        end if

    end subroutine destructorTypes__

end module cfdtypes