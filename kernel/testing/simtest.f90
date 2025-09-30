program simTest
!
!   main program to run compressible fluid FIE simulation
!
    use cfd
    use logger
    implicit none

    ! water system test parameters
    class(CFD), allocatable   :: model
    character(5)              :: type = "water"
    integer                   :: modelSize = 40!= N <=> manifold(N,N,N)
    real(dp)                  :: volume = 1._dp, timestep = 0.002_dp, spHeat = 4.186_dp
    logical                   :: useGPU = .false.
    integer                   :: i, tmax
  

    !
    !   timestep here is below CFL condition estimation of ~0.003-0.004
    !

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

    ! test 3 separate time
    
    ! 1: t = 1  (initial state in fortran)
    call model%writePressureStateCPU()
    
    tmax = 100


    call logger_start()
    do i = 1,tmax
        call model%update()
        if (mod(i,2) .eq. 0) call model%writePressureStateCPU()
    end do
    call logger_stop()

    print*, "\n"
    write(*,'(A25,F12.6,A)') "Total Computation Time: ", get_total_time(), " seconds"
    write(*,'(A25,F12.6,A)') "Simulation Timesteps: ", model%clock / model%timestep, " steps"
    write(*,'(A25,F12.6,A)') "Simulation Duration: ", model%clock, " seconds"
    print*, "\n"

    call logger_reset()

    if (allocated(model)) deallocate(model)



end program simtest