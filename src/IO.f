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

integer(kind=dp) :: hack_switch, start_index, periapsis_index, dI, p1,p2



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


!First, do a run through to get the 'start point'
do j=1,rows

    if (r(j) .GT. semi_major) then

    start_index = j
    exit

    endif


enddo

!periapsis_index = minloc(r(start_index:))




print *, 'minloc', minloc(r(start_index:))

periapsis_index = 14820
!This is an articifial estimate of the index at periapsis
!Could read this each time,but we want the times to be the same between runs
!To be refined later

dI =  rows - periapsis_index

p1 = periapsis_index - dI/5
p2 = periapsis_index + dI/5




open(unit=20,file=savefile_targets,status='replace',form='formatted')





div = (p1-start_index + rows - p2)/(N_targets/4.0_dp)!division if using even sampling

do j = start_index,p1

    if (mod(j,div) .eq. 0.0_dp) then

    write(20,*) array(j,1), x(j), y(j),z(j),&
                u0(j), u1(j), u2(j),u3(j), &
                r(j), theta(j), phi(j),tau(j)
    endif

enddo



div = (p2-p1)/(3.0_dp*N_targets/4.0_dp)
do j = p1,p2

    if (mod(j,div) .eq. 0.0_dp) then

    write(20,*) array(j,1), x(j), y(j),z(j),&
                u0(j), u1(j), u2(j),u3(j), &
                r(j), theta(j), phi(j),tau(j)
    endif

enddo



div = (p1-start_index + rows - p2)/(N_targets/4.0_dp)!division if using even sampling
do j = p2,rows

    if (mod(j,div) .eq. 0.0_dp) then
    write(20,*) array(j,1), x(j), y(j),z(j),&
                u0(j), u1(j), u2(j),u3(j), &
                r(j), theta(j), phi(j),tau(j)
    endif

enddo


!do j = 1,rows



!if (r(j) .GT. semi_major) then
!hack_switch = 1
!endif

! hack to ignore the start of the orbit
!if (mod(j , div) .eq. 0.0_dp .and. hack_switch .EQ. 1)  then





!write(20,*) array(j,1), x(j), y(j),z(j),&
 !           u0(j), u1(j), u2(j),u3(j), &
  !          r(j), theta(j), phi(j),tau(j)


!call calculate_covariant_metric(r(j),theta(j),metric)



!print *, '----------------'
!print *, 'Check points:'
!print *, r(j), theta(j)
!print *, metric(1,1), metric(2,2),metric(3,3),metric(4,4),metric(4,1)
!print *, u0(j), u1(j), u2(j),u3(j)


!print *, 'check:', metric(1,1) * u0(j)**2 &
 !                + metric(2,2) * u1(j)**2 &
  !               + metric(3,3) * u2(j)**2 &
   !              + metric(4,4) * u3(j)**2 &
    !             +2.0_dp*metric(4,1)*u0(j)*u3(j)


!endif

!enddo




close(20)



end subroutine IO_array



end module IO
