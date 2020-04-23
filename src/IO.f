module IO


use parameters
use constants
use tensors

implicit none


private

public IO_array


contains


subroutine IO_array(array,array_derivatives,rows,cols)
!Arguments
integer(kind=dp), intent(in) :: rows,cols
real(kind=dp), intent(in) :: array(rows,cols)
real(kind=dp), intent(in) :: array_derivatives(rows,4)

!Other

!Create row vectors to load the data
real(kind=dp), dimension(rows) :: r, theta, phi,s1,s2,s3,tau,u0,u1,u2,u3

!Cartesian components 
real(kind=dp), dimension(rows) :: sx, sy,sz,x,y,z

!Nutation and precession angles
real(kind=dp), dimension(rows) :: thetaSL, phiSL

!Create an integer for iterating
integer(kind=dp) :: j

!Critical phase angle
real(kind=dp),dimension(rows) :: top, bot,critical_phase 

!Create an integer for extracting the right number of target points
integer(kind=dp) :: div


!Create a column vector to calculate 4-velocities
real(kind=dp), dimension(12) :: yvector,dyvector


!Create an array to hold the metric
real(kind=dp), dimension(4,4) :: metric



!Load the data - vectorised
r = array(:,2) ; theta = array(:,3); phi = array(:,4)
s1 = array(:,10) ; s2 = array(:,11); s3=array(:,12)
tau = array(:,13) / convert_s

!And the 4-velocity
u0 = array_derivatives(:,1) ; u1 = array_derivatives(:,2) ; u2 = array_derivatives(:,3) ; u3 = array_derivatives(:,4)


!Calculate the cartesian cpts
x = sqrt(r**2 + a**2)*sin(theta)*cos(phi)
y = sqrt(r**2 + a**2)*sin(theta)*sin(phi)
z = r * cos (theta)




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



!Write to file - spin stuff
open(unit=10,file=savefile2,status='replace',form='formatted')
do j = 1,rows
write(10,*) tau(j),thetaSL(j),phiSL(j), sx(j), sy(j), critical_phase(j), tau(j)/(PeriodEst / convert_s)
enddo
close(10)



!Write target points to file
div = rows/N_targets
print *, 'div = ',div
open(unit=20,file=savefile_targets,status='replace',form='formatted')
do j = 1,rows


if (mod(j , div) .eq. 0.0_dp)  then
write(20,*) array(j,1), x(j), y(j),z(j),&
            u0(j), u1(j), u2(j),u3(j), &
            r(j), theta(j), phi(j)


call calculate_covariant_metric(r(j),theta(j),metric)



print *, '----------------'
print *, 'Check points:'
print *, r(j), theta(j)
print *, metric(1,1), metric(2,2),metric(3,3),metric(4,4),metric(4,1)
print *, u0(j), u1(j), u2(j),u3(j)


print *, 'check:', metric(1,1) * u0(j)**2 &
                 + metric(2,2) * u1(j)**2 &
                 + metric(3,3) * u2(j)**2 &
                 + metric(4,4) * u3(j)**2 &
                 +2.0_dp*metric(4,1)*u0(j)*u3(j)


endif

enddo
close(20)



end subroutine IO_array



end module IO
