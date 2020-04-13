module parameters
implicit none


!Define float precision
!integer, parameter :: dp = selected_real_kind(33,4931)
integer, parameter :: dp = selected_real_kind(15, 307)
!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!Orbital parameters

real(kind=dp), parameter :: semi_major = 20.0_dp 
real(kind=dp), parameter :: eccentricity = 0.80_dp !Orbital eccentricity
real(kind=dp), parameter :: iota = 35.0_dp !Inclination w.r.t equatorial plane in degrees
real(kind=dp), parameter :: N_orbit = 30.50_dp !Number of orbits to integrate


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 4.31d6!BH mass in solar masses
real(kind=dp), parameter :: a= +0.9980_dp !BH spin parameter Now set later

!PSR intrinsic parameters
real(kind=dp), parameter :: MPSR = 1.40_dp !pulsar mass in solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !pulsar radius in km
real(kind=dp), parameter :: stheta = PI/4.0_dp, sphi = PI/4.0_dp


real(kind=dp), parameter :: lambda = 1.0_dp !Turn on/off spin-curvature coupling (1 = on)
real(kind=dp), parameter :: p0 = 1.0e-3 !spin period
real(kind=dp), parameter :: eta = 3.0_dp*PI/12.0_dp !Oreintation of initial momentum

!IO location
character(len=200) :: IO_path = '/Users/tomkimpson/Data/ThesisData/MPD/'


!Debugging
integer(kind=dp), parameter :: print_status = 1 !Turns on/off 1/0 print commands 



end module parameters
