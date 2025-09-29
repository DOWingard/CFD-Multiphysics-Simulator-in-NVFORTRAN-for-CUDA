program CFluid
!
!   main program to run compressible fluid FIE simulation
!
use cell
implicit none





! construct the cell for the simulation
type(cell) :: FEcell
real(dp)   :: volume
real(dp)   :: dU(5)
real(dp)   :: config(6)
integer    :: location(3)
! specify the fluid initial state
config = (/ &
            1.0e5_dp,  &! (p      : pressure)
            7.00_dp,   &! (\gamma : specific heat ratio)  1.4 for air
            1000.0_dp, &! (\rho   : density)
            0.0_dp,    &! (x-velocity)
            0.0_dp,    &! (y-velocity)
            0.0_dp     &! (z-velocity)
         /)

volume   = 1._dp
dU       = 1._dp
location = (/1,1,1/)

call fecell%setLocation(location)
call fecell%init(config,volume)
!call fecell%update(dU)



print*, fecell%location, fecell%U





end program cfluid