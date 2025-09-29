module noether
    use numtype
    implicit none

    public :: noether

    !-------------------------------
    ! Abstract kernel type
    !-------------------------------
    type, abstract :: setKernel
        integer :: size
        real(dp), allocatable :: U(:)
    contains
        procedure(noether_func), deferred :: init
    end type setKernel

    !-------------------------------
    ! Concrete kernel
    !-------------------------------
    type :: noether
        real(dp), allocatable :: U(:)
        procedure(noether_func), pointer :: init_ => null()
    contains
        procedure :: constructor
        procedure :: init
    end type noether

    !-------------------------------
    ! Abstract interface for init
    !-------------------------------
    abstract interface
        subroutine noether_func(self, config)
            class(setKernel), intent(inout) :: self
            real(dp), intent(in)            :: config(:)
        end subroutine noether_func
    end interface

contains

    !=======================
    subroutine constructor(this, kernel)
        class(noether), intent(inout) :: this
        class(setKernel), intent(in)  :: kernel

        if (.not. allocated(this%U)) allocate(this%U(kernel%size))
        this%U = 0._dp

        this%init_ => kernel%init
    end subroutine constructor

    !=======================
    subroutine init(this, config)
        class(noether), intent(inout) :: this
        real(dp), intent(in)          :: config(:)

        if (associated(this%init_)) then
            call this%init_(this, config)
        else
            print *, "Error: init procedure not assigned."
        end if
    end subroutine init

end module noether
