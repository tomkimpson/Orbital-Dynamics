module rungekutta

use parameters
use constants
use derivatives
use analysis
use IO

implicit none

private


public rk

contains


subroutine rk(y0)

 
!Arguments
real(kind=dp),intent(IN),dimension(entries) :: y0 !initial conditions

!Other
real(kind=dp), dimension(size(y0)) :: y, y1,dy !Some useful vectors used during integration
real(kind=dp), dimension(nrows,ncols) :: AllData !Big array to save all data. 12 coordinates + tau +
real(kind=dp), dimension(nrows,4) :: AllDerivs !Big array to save all data. 12 coordinates + tau +
real(kind=dp), dimension(:,:),allocatable :: output, outputDerivs !smaller array which will be outout
integer(kind=dp) :: i,j,NSteps !,nsteps !index for saving to array


real(kind=dp) :: tau, ur
real(kind=dp) :: mm, xC, yC, zC !Cartesian components
real(kind=dp),dimension(3) :: xi, xi_B !Cartesian cpts



!Create some unit tangent vector for testing in place of the ray tracing solution
!And its comoving version
real(kind=dp), dimension(4) :: k, k_tetrad

!And the cartesian 3 vector
real(kind=dp), dimension(3) :: k_cartesian

!Create some vectors for reading in position, velocity etc
real(kind=dp), dimension(4) :: vector_position, vector_velocity

!The normal vector
real(kind=dp),dimension(3) :: vector_normal
real(kind=dp),dimension(4) :: n,n_BL, n_tetrad


!different pitch angles
real(kind=dp) :: pitch_coordinate, pitch_tetrad, pitch_coordinate2,mag

!Logicals
logical :: condition



!Get the cartesian vector in the global frame
k_cartesian(1) = sin(ThetaObs)*cos(PhiObs) 
k_cartesian(2) = sin(ThetaObs)*sin(PhiObs) 
k_cartesian(3) = cos(ThetaObs) 

!Set the integration tolerance
if (dp .EQ. 8) then
escal = 1.0d15
else if (dp .EQ. 16) then
escal = 1.0d19
endif


!Assign y0 to vecor y
y  = y0
tau = 0.0_dp

!Save the first row to array
i = 1
AllData(i,1:12) = y
AllData(i,13) = tau !tau

!Get some derivative info
call derivs(y,dy)
AllDerivs(i,1:4) = dy(1:4)

!Integrate

if (duration .EQ. 'long') then
    condition = y(1) .LT. time_cutoff
elseif (duration .EQ. 'short') then
    condition = tau .LT. N_spins
else
    print *, 'Integration condition not specified'
    stop
endif



!do while condition

do while (condition)

    !Update
    call RKF(y,y1)
    y = y1
 
    !Save the output
    i = i + 1
 

    if (i .GT. nrows) then
    print *, 'i GT nrows'
    exit
    endif


    AllData(i,1:12) = y
    tau = tau + h
    AllData(i,13) = tau


    call derivs(y,dy)
    AllDerivs(i,1:4) = dy(1:4)



!Update condition
if (duration .EQ. 'long') then
    condition = y(1) .LT. time_cutoff
elseif (duration .EQ. 'short') then
    condition = tau .LT. N_spins
else
    print *, 'Integration condition not specified'
    stop
endif




enddo
NSteps = i




print *, 'Total number of steps is = ', i
print *, 'Runge Kutta completed. Start data I/O'
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!



!Binary format. See discussion at https://stackoverflow.com/questions/24395686/best-way-to-write-a-large-array-to-file-in-fortran-text-vs-other



!First reallocate to create a smaller array
allocate(output(NSteps,ncols))
output = AllData(1:i, :)

allocate(outputDerivs(NSteps,4))
outputDerivs = AllDerivs(1:i, :)


!Save the spatial trajectory
open(unit=30,file=savefile1,status='replace',form='formatted')
do j=1,NSteps

    !Read in the 4-poisiton and 4-velocity
    vector_position(1:4) = output(j,1:4)
    vector_velocity(1:4) = outputDerivs(j,1:4)


    !Modify to zamo velocity if you want, just for testing
    !call zamo(vector_position,vector_velocity)

    !Calculate the beam vector and normalise it to be a unit vector
    call transform_to_radiation_point(output(j,1:4) , output(j,9:12), output(j,13)/convert_s, xi, xi_B)
    vector_normal = xi_B - xi

    mag = sqrt(vector_normal(1)**2 + vector_normal(2)**2 + vector_normal(3)**2)
    vector_normal = vector_normal / mag !normalise

    !Pitch angle in the coordinate frame
    call pitch_angle(vector_normal,k_cartesian,pitch_coordinate)


    !Transform the vector normal to BL coordinates
    n(1) =1.0_dp
    n(2:4) = vector_normal
    call change_coordinates(n,n_BL,vector_position,'contravariant', 'BL')

    !Change it back to check transform is ok
    !call change_coordinates(n_BL,n,vector_position,'contravariant', 'cartesian')


    !Can also change it to a covariant vector
    !call lower_index(n_BL,n_BL,vector_position)


 
    !And transform it to the comoving frame
   call transform_to_comoving_frame(n_BL,vector_velocity,vector_position,n_tetrad,'contravariant')


   !Create a ray at the position of the radiation point
    call create_ray(xi,k)

   !Calculate the comoving transform of k
   call transform_to_comoving_frame(k,vector_velocity,vector_position,k_tetrad,'covariant')

    
    !Make the spatial parts of both the comoving normal and the comoving ray unit vectors
    mag = sqrt(n_tetrad(2)**2 + n_tetrad(3)**2 + n_tetrad(4)**2)
    n_tetrad = n_tetrad / mag


    mag = sqrt(k_tetrad(2)**2 + k_tetrad(3)**2 + k_tetrad(4)**2)
    k_tetrad = k_tetrad / mag

    !Calculate the comoving pitch angle
    call pitch_angle(n_tetrad(2:4),k_tetrad(2:4),pitch_tetrad)


write(30,*) output(j,1)/convert_s,xi(1),xi(2),xi(3), &
            xi(1) / convert_m, xi(2)/convert_m, xi(3)/convert_m, &
            xi_B(1), xi_B(2), xi_B(3), &
            xi_B(1)/convert_m, xi_B(2)/convert_m, xi_B(3)/convert_m, &
            pitch_coordinate, pitch_tetrad, &
            vector_position
enddo
close(30)


print *, output(Nsteps,:)

call IO_array(output,Nsteps,ncols)

!open(unit=40,file=savefile2,status='replace',form='formatted')

!do j=1,NSteps

!enddo

!close(40)





!Save the time delays


end subroutine rk


end module rungekutta
