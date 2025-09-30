module tangent
!
!   simple module to contain
!   the type for tangent bundle
!
    use numtype
    implicit none 

    type flux 
        real(dp), allocatable :: fluxes(:,:)
    contains
        procedure :: setSize
        
        final     :: destructorTangent__
    end type flux


    type bundle 
        class(flux), allocatable :: mesh(:,:,:)
    contains
        procedure :: alloc 
    end type bundle

contains 

    subroutine alloc(this, size, nvars)
    !
    !   allocate whole bundle
    !
        class(bundle), intent(inout) :: this
        integer, intent(in) :: size
        integer, intent(in) :: nvars   ! number of conserved variables
        integer :: i,j,k

        allocate(this%mesh(size,size,size))

        !=========================================================
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k)
        !=========================================================
        do k = 1, size
        do j = 1, size
        do i = 1, size
            call this%mesh(i,j,k)%setSize(nvars)
        end do
        end do
        end do
        !=========================================================
        !$OMP END PARALLEL DO
        !=========================================================

    end subroutine alloc

    subroutine setSize(this, nvars)
    !
    !   allocate size
    !
        class(flux), intent(inout) :: this
        integer, intent(in)        :: nvars

        if (allocated(this%fluxes)) deallocate(this%fluxes)
        allocate(this%fluxes(6,nvars))
        this%fluxes = 0.0_dp      ! initialize elements at each point

    end subroutine setsize

    subroutine destructorTangent__(this)
    !
    !   deallocate the cells
    !
        type(flux), intent(inout) :: this

        if (allocated(this%fluxes)) then
            deallocate(this%fluxes)
        end if
    end subroutine destructortangent__

end module tangent
