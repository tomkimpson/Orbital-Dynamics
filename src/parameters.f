module parameters
implicit none


!Define float precision
!integer, parameter :: dp = selected_real_kind(33,4931) !quad precision
integer, parameter :: dp = selected_real_kind(15, 307) !double precision
!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!Orbital parameters

real(kind=dp), parameter :: semi_major = 200.0_dp 
real(kind=dp), parameter :: eccentricity = 0.10_dp !Orbital eccentricity
real(kind=dp), parameter :: iota = 0.0_dp !Inclination w.r.t equatorial plane in degrees
real(kind=dp), parameter :: N_orbit = 1.0_dp !Number of orbits to integrate. Used if duration = 'long' (below)
real(kind=dp), parameter :: N_spins = 2.0_dp !Number of spin periods to integrate. Used if duration = 'short' (below)


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 4.3100d6!BH mass in solar masses
real(kind=dp), parameter :: a= 0.60_dp !BH spin parameter Now set later

!PSR intrinsic parameters
real(kind=dp), parameter :: MPSR = 1.40_dp !pulsar mass in solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !pulsar radius in km
real(kind=dp), parameter :: stheta = PI/6.0_dp, sphi = 0.0_dp !Initial oreintatin of spin axis
real(kind=dp), parameter :: psi = PI/6.0_dp !Polar angle of radiation beam w.r.t spin axis
real(kind=dp), parameter :: p0 = 1.0e-3 !spin period

!Obaserver Location
real(kind=dp), parameter :: ThetaObs = PI/4.0_dp , PhiObs = 0.0_dp

!Some additional settings
real(kind=dp), parameter :: lambda = 1.0_dp !Turn on/off spin-curvature coupling (1 = on)
real(kind=dp), parameter :: eta = 3.0_dp*PI/12.0_dp !Oreintation of initial momentum


!Integration settings
integer(kind=dp), parameter :: adaptive = 1 !turn on/off adaptive stepsize
character(len=20), parameter :: duration = 'long' !long, short. Long integrates for Norbits * period. !Short integrates for N_spins * p0


!IO settings
character(len=200) :: IO_path = '/Users/tomkimpson/Data/ThesisData/MPD/'
integer(kind=dp) :: N_targets = 20 !Number of target points to extract to use with Ray Tracing

!Debugging
integer(kind=dp), parameter :: print_status = 1 !Turns on/off 1/0 print commands 



end module parameters
