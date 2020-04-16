module IO


use parameters
use constants


implicit none


private

public IO_array


contains


subroutine IO_array(array,rows,cols)
!Arguments
integer(kind=dp), intent(in) :: rows,cols
real(kind=dp), intent(in) :: array(rows,cols)

!Other

!Create row vectors to load the data
real(kind=dp), dimension(rows) :: r, theta, phi,s1,s2,s3,tau

!Cartesian components 
real(kind=dp), dimension(rows) :: sx, sy,sz

!Nutation and precession angles
real(kind=dp), dimension(rows) :: thetaSL, phiSL

!Create an integer for iterating
integer(kind=dp) :: j

!Critical phase angle
real(kind=dp),dimension(rows) :: top, bot,critical_phase 



!Load the data - vectorised
r = array(:,2) ; theta = array(:,3); phi = array(:,4)
s1 = array(:,10) ; s2 = array(:,11); s3=array(:,12)
tau = array(:,13) / convert_s


!Some calculations with the spin components
sx = s1*sin(theta)*cos(phi) + s2*r*cos(theta)*cos(phi) - s3*r*sin(theta)*sin(phi)
sy = s1*sin(theta)*sin(phi) + s2*r*cos(theta)*sin(phi) + s3*r*sin(theta)*cos(phi)
sz = s1*cos(theta) - s2*r*sin(theta)

thetaSL = atan2(sqrt(Sx**2 + Sy**2),Sz)
phiSL = atan2(Sy,Sx)


!Critical phase angle for specific case
top = cos(phiSL)*cos(thetaSL) - sin(thetaSL)
bot = Sqrt(cos(phiSL)**2 * cos(thetaSL)**2 + sin(phiSL)**2 + sin(thetaSL)**2 - cos(phiSL)*sin(2.0_dp*thetaSL))

critical_phase=acos(top/bot)

!Timing delay due to shift of centre of pulse profile
!Only holds for psi = PI/4, 



!Write to file
open(unit=10,file=savefile2,status='replace',form='formatted')

do j = 1,rows


write(10,*) tau(j),thetaSL(j),phiSL(j), sx(j), sy(j), critical_phase(j), tau(j)/(PeriodEst / convert_s)
enddo


close(10)



end subroutine IO_array



end module IO
