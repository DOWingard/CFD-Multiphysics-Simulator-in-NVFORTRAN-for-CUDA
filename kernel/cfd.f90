module cfd
!
!   core module to access the fluid
!   module and calculate the 
!   CPU vs GPU methods
!
    use fluid
    use cfdtypes
    implicit none


    type :: CFD

        class(manifold), allocatable   :: fluid
        class(cfdtype), allocatable    :: cfdtypes 
        character(len=20)              :: type

        ! sim parameters
        integer   :: meshSize
        real(dp)  :: timestep, volume, specfHeat
        real(dp)  :: clock = 0
        logical   :: useGPU = .false.
    contains

        !procedure :: addCFDtype
        procedure :: setup                  ! Usage: (str fluidtype, int N, double V, double dt, double gamma, bool isCUDA)
        procedure :: init                   ! Usage: ()
        procedure :: setState               ! Usage: ()
        procedure :: update                 ! Usage: ()
        procedure :: writePressureState     ! Usage: ()
        ! TO DO
        !procedure :: setFEbasis             ! Usage: ()


        ! private methods
        procedure :: writePressureStateCPU 
        procedure :: loadTangentBundleCPU
        procedure :: updateFromTangentBundleCPU

        final     :: destructorCFD__

    end type cfd


contains

    subroutine setup(this, fluidtype, N, V, dt, gamma, isCUDA)
        !
        !   set up the simulation for 
        !   CFD type
        !
        class(cfd), intent(inout)    :: this
        real(dp), intent(in)         :: gamma, V, dt
        integer, intent(in)          :: N
        character(len=*), intent(in) :: fluidtype
        logical, intent(in)          :: isCUDA
        ! private
        logical :: found
        integer :: i

        ! Ensure CFD type object is allocated
        if (.not. allocated(this%cfdtypes)) then
            allocate(this%cfdtypes)
            call this%cfdtypes%init()
        end if

        found = .false.
        do i = 1, size(this%cfdtypes%fluids)
            if (trim(this%cfdtypes%fluids(i)) == trim(fluidtype)) then
                found = .true.
                exit
            end if
        end do


        if (.not. found) then
            write(*,*) "Error: CFD type '", trim(fluidtype), &
                       "' is not recognized."
            stop 1
        end if

        if (mod(N,2) /= 0) then
            write(*,*) "Error: meshSize N must be even"
            stop 1
        end if

        if (dt <= 0) then
            write(*,*) "Error: timestep must be >0"
            stop 1
        end if

        this%type      = fluidtype
        this%meshSize  = N
        this%volume    = V
        this%specfHeat = gamma
        this%useGPU    = isCUDA
        this%timestep  = dt

    end subroutine setup



    subroutine init(this)
        !
        !   initialize the CFD system
        !
        class(cfd), intent(inout) :: this
        ! private 
        logical  :: valid
        integer  :: i

        if (len_trim(this%type) == 0) then
            write(*,*) "Error: CFD type has not been set. Call setup() first."
            stop 1
        end if

        
        valid = .false.
        do i = 1, size(this%cfdtypes%fluids)
            if (trim(this%cfdtypes%fluids(i)) == trim(this%type)) then
                valid = .true.
                exit
            end if
        end do


        if (.not. valid) then
            write(*,*) "Error: CFD type '", trim(this%type), &
                       "' is not recognized."
            stop 1
        end if

        if (this%meshSize <= 0) then
            write(*,*) "Error: mesh size N has not been set. Call setup() first."
            stop 1
        end if

        if (this%specfHeat <= 0._dp) then
            write(*,*) "Error: specific heat gamma has not been set. Call setup() first."
            stop 1
        end if

        ! Ensure fluid object is allocated
        if (.not. allocated(this%fluid)) then
            allocate(this%fluid)
        end if

        if (.not. this%useGPU) then 
            if (isDebug .eq. 1) print*, "\nUSING CPU for CFD%init"
            call this%fluid%fillManifoldCPU(this%meshSize, this%volume)
            call this%fluid%allocBundleCPU(this%meshSize)
        else
            if (isDebug .eq. 1) print*, "\nCUDA ACTIVE: USING GPU for CFD%init (not yet implemented)"
        end if

        write(*,*) "\nCFD system initialized for type:", trim(this%type), "\n"


    end subroutine init


    subroutine setState(this)
        !
        !   sets state of the system
        !   >CURRENTLY: used for initializing pressure gradient
        !
        class(cfd), intent(inout) :: this

        if (this%meshSize <= 0) then
            write(*,*) "Error: mesh size N has not been set. Call setup() first."
            stop 1
        end if

        call this%fluid%initDensityGaussianCPU(this%meshSize)

    end subroutine setState

    subroutine update(this)
    !
    !   complete update, use private_ methods for 
    !   building/testing/debugging
    !
        class(cfd), intent(inout) :: this

        if (.not. this%useGPU) then
            if (isDebug .eq. 1) print*, "\nBefore update: ", this%fluid%fluid(this%meshSize/2, this%meshSize/2, this%meshSize/2)%U(1)
            if (isDebug .eq. 1) print*, "Using CPU for CFD%update"

            ! Build and apply update
            call this%loadTangentBundleCPU()
            call this%updateFromTangentBundleCPU()

            if (isDebug .eq. 1) print*, "After update: ", this%fluid%fluid(this%meshSize/2, this%meshSize/2, this%meshSize/2)%U(1)
        else
            if (isDebug .eq. 1) print*, "\nCUDA ACTIVE: USING GPU for CFD%update (not yet implemented)"
        end if

    end subroutine update



    subroutine writePressureState(this)
    !
    !   write pressure state out for 
    !   plotting
    !   
        class(cfd), intent(in)         :: this

        if (.not. this%useGPU) then
            if (isDebug .eq. 1) print*, "CPU for cell pressures"
            call writePressureStateCPU(this)
        end if 

        if ( this%useGPU ) then
            if (isDebug .eq. 1) print*, "CUDA + GPU for cell pressures" 
        ! TODO :   this%writePressureStateGPU_()
        end if


    end subroutine writepressurestate

    subroutine writePressureStateCPU(this)
    !
    !   write pressure state of system to data/
    !
        class(cfd), intent(in) :: this
        ! private
        integer :: timestep_idx
        integer :: i,j,k
        integer :: unit
        integer :: ios
        character(len=64) :: filename
        real(dp), allocatable :: pressures(:,:,:)
        
        timestep_idx = nint(this%clock / this%timestep)
        write(filename,'(A,I0)') 'data/fort.', timestep_idx

        allocate(pressures(this%meshSize,this%meshSize,this%meshSize))

        
        !=============================================================
        !$OMP PARALLEL DO collapse(3) private(i,j,k) shared(pressures)
        !=============================================================
        do k = 1, this%meshSize
            do j = 1, this%meshSize
                do i = 1, this%meshSize
                    pressures(i,j,k) = this%fluid%fluid(i,j,k)%cellPressure()
                end do
            end do
        end do
        !=============================================================
        !$OMP END PARALLEL DO
        !=============================================================



        open(newunit=unit, file=filename, status='replace', action='write', form='formatted', iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error opening file:", filename
            return
        end if

        ! thread safe to do this in series :(
        do k = 1, this%meshSize
            do j = 1, this%meshSize
                do i = 1, this%meshSize
                    write(unit,*) i, j, k, pressures(i,j,k)
                end do
            end do
        end do

        close(unit)
        deallocate(pressures)

    end subroutine writePressureStateCPU



    subroutine loadTangentBundleCPU(this)
        !
        !   loads tangent bundle with 
        !   flux values
        !
        class(cfd), intent(inout) :: this

        if (this%meshSize <= 0) then
            write(*,*) "Error: mesh size N has not been set. Call setup() first."
            stop 1
        end if

        if (this%specfHeat <= 0._dp) then
            write(*,*) "Error: specific heat gamma has not been set. Call setup() first."
            stop 1
        end if

        call this%fluid%loadBundleCPU(this%specfHeat, this%meshSize)

    end subroutine loadTangentBundleCPU


    subroutine updateFromTangentBundleCPU(this)
        !
        !   updates the manifold from the tangent Bundle
        !
        class(cfd), intent(inout) :: this

        call this%fluid%updateCPU(this%timestep)
        this%clock = this%clock + this%timestep

    end subroutine updateFromTangentBundleCPU


    !subroutine setFEbasis(this) 
    !
    !   calls manifold%computeBasis 
    !
    !    class(cfd), intent(inout) :: this

    !    call this%computeBasis()

    !end subroutine setfebasis






    subroutine destructorCFD__(this)
        !
        !   IDK if necessary but best practice
        !   
        type(cfd), intent(inout) :: this

        if (allocated(this%fluid)) deallocate(this%fluid)
        if (allocated(this%cfdtypes)) then
            if (allocated(this%cfdtypes%fluids)) deallocate(this%cfdtypes%fluids)
            deallocate(this%cfdtypes)
        end if

    end subroutine destructorcfd__ 


end module cfd
