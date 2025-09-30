program simTest
!
!   main program to run compressible fluid FIE simulation
!
    use cfd
    implicit none

    ! water system test parameters
    class(CFD), allocatable   :: model
    character(5)              :: type = "water"
    integer                   :: modelSize = 50!= N <=> manifold(N,N,N)
    real(dp)                  :: volume = 1._dp, timestep = 0.01_dp, spHeat = 4.186_dp
    logical                   :: useGPU = .false.


    !type :: CFD
    !    class(manifold), allocatable   :: fluid
    !    class(cfdtype), allocatable    :: cfdtypes 
    !    character(len=20)              :: type
    !    ! sim parameters
    !    integer   :: meshSize
    !    real(dp)  :: timestep, volume, specfHeat
    !    real(dp)  :: clock = 0
    !    logical   :: useGPU = .false.
    !contains
    !    procedure :: setup         ! Usage: (str fluidtype, int N, double V, double dt, double gamma, bool isCUDA)
    !    procedure :: init          ! Usage: ()
    !    procedure :: setState      ! Usage: ()
    !    procedure :: update        ! Usage: ()
    !end type cfd


    allocate(model)

    call model%setup(type, modelSize, volume, timestep, spHeat, useGPU)
    call model%init()
    call model%setState()
    call model%update()

    if (allocated(model)) deallocate(model)



end program simtest