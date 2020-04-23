module analysis

use parameters
use constants
use tensors

implicit none

private 


public transform_to_comoving_frame,&
       PostKeplerianDelays,&
       transform_BL_to_XY,&
       transform_to_radiation_point,&
       covariant_transform,&
       pitch_angle,&
       create_ray,&
       zamo,&
       change_coordinates, &
       lower_index

contains


!transform_BL_to_XY. Takes a 3 vector in BL coordinates and transforms to cartesian coordinates

!create_ray. Creates an artificial ray using the initial conditions mathematics from ray tracing. Uses observerlocation from parameters.f


!zamo. creates a zero angular momentum contravariant 4 velocity


!transform_to_comoving frame. Takes a 4-vector in the global frame, transforms it to the comovign frame at position x, 4 vecloty v


!pitch_angle. Takes two spatial 3-vectors to compute the angle between them





subroutine zamo(x,u)
!Create zamo at 4-position x
!Outputs contravariant 4-velocity u
real(kind=dp), intent(in),dimension(4) :: x
real(kind=dp), intent(out),dimension(4) :: u

!Other
real(kind=dp) ::  r,theta,phi
real(kind=dp),dimension(4,4) :: metric
real(kind=dp) :: omega, norm,delta


!Load coordinates
r = x(2) ; theta = x(3) ; phi = x(4)

!Calculate the covariant metric
call calculate_covariant_metric(r,theta,metric)


!Construct the zamo
omega = -metric(1,4) / metric(4,4)
delta = r**2 -2.0_dp*r + a**2
norm = sqrt(metric(4,4)/ (delta*sin(theta)**2))

u(1) = 1.0_dp ; u(2) = 0.0_dp ; u(3) = 0.0_dp ; u(4) = omega

u = u*norm


end subroutine





subroutine lower_index(v,w,x)
!Arguments
real(kind=dp), dimension(4), intent(in) :: v !input contravariant vector
real(kind=dp), dimension(4), intent(out) :: w !output covariant vector
real(kind=dp), dimension(4), intent(in) :: x !input position

!other
real(kind=dp), dimension(4,4) :: metric
real(kind=dp) :: r,theta,phi


r = x(2) ; theta = x(3) ; phi = x(4)


call calculate_covariant_metric(r,theta,metric)


w(1) = metric(1,1) * w(1) + metric(1,2) * w(2) + metric(1,3) * w(3) + metric(1,4) * w(4)
w(2) = metric(2,1) * w(1) + metric(2,2) * w(2) + metric(2,3) * w(3) + metric(2,4) * w(4)
w(3) = metric(3,1) * w(1) + metric(3,2) * w(2) + metric(3,3) * w(3) + metric(3,4) * w(4)
w(4) = metric(4,1) * w(1) + metric(4,2) * w(2) + metric(4,3) * w(3) + metric(4,4) * w(4)






end subroutine lower_index


subroutine change_coordinates(v,w,x,vector_type,goal)

!Change the coordinates of a vector v, of type (covar/contra) vector_type.
!Goal specifies direction e.g. goal = 'BL' ---> transform from cartesian to BL 

!Arguments
real(kind=dp), dimension(4), intent(in) :: v
real(kind=dp), dimension(4), intent(in) :: x !position in BL
real(kind=dp), dimension(4), intent(out) :: w
character(len=*),intent(in) :: vector_type, goal

!Some useful stuff
real(kind=dp) :: r,theta,phi,sigma, mm

w(1) = v(1) !Time components does not change

!Read in position
r = x(2); theta = x(3); phi = x(4)

!Define some useful quantities
sigma = r**2 + (a*cos(theta))**2
mm = sqrt(r**2 + a**2)


if (vector_type .EQ. 'contravariant') then

    if (goal .EQ. 'BL') then
    !Transform a contravariant vector from cartesian to BL coordinates
    w(2) = (mm*r/sigma) *sin(theta)*cos(phi) * v(2) + (mm*r/sigma) * sin(theta)*sin(phi) * v(3) + mm**2/sigma *cos(theta)*v(4)

    w(3) = mm/sigma * cos(theta)*cos(phi) * v(2) + mm/sigma * cos(theta)*sin(phi) * v(3) -r/sigma * sin(theta) * v(4)   
    
    w(4) = -sin(phi)/(mm*sin(theta)) * v(2) + cos(phi)/(mm*sin(theta)) * v(3)


    elseif (goal .EQ. 'cartesian') then

    
    w(2) = r/mm * sin(theta)*cos(phi) * v(2) + mm*cos(theta)*cos(phi)*v(3) - mm*sin(theta)*sin(phi)*v(4)
    w(3) = r/mm * sin(theta)*sin(phi) * v(2) + mm*cos(theta)*sin(phi)*v(3) + mm*sin(theta)*cos(phi)*v(4)
    w(4) = cos(theta)*v(2) -r*sin(theta)*v(3)


    endif

endif





end subroutine change_coordinates





subroutine transform_BL_to_XY(vBL, vXY,r,theta,phi)
!Arguments
real(kind=dp), dimension(3),intent(in) :: vBL
real(kind=dp), dimension(3),intent(inout) :: vXY
real(kind=dp), intent(in) :: r,theta,phi

!Other
real(kind=dp) :: m

m = sqrt(r**2 + a**2)



end subroutine transform_BL_to_XY


subroutine create_ray(xV,p)
!Input = 3 position vector
!Output = tangent ray 4 vector, p

!Arguments
real(kind=dp),intent(in), dimension(3) :: xV
real(kind=dp),intent(out), dimension(4) :: p

!Other
real(kind=dp) :: xdot, ydot, zdot
real(kind=dp) :: rdot0, thetadot0, phidot0
real(kind=dp) :: mm, sigma
real(kind=dp) :: sig, del
real(kind=dp) :: s1, Eobs,Eprime,E2,En,pr,ptheta,Lz
real(kind=dp) :: r,theta,phi
real(kind=dp) :: x,y,z,w    



!Load initial position from parameters file
x= xV(1); y = xV(2) ; z = xV(3)

!Cartesian direction of ray from angles in parameters.f
xdot = sin(ThetaObs)*cos(PhiObs)
ydot = sin(ThetaObs)*sin(PhiObs)
zdot = cos(ThetaObs)


!Convert both to BL
w = x**2+y**2+z**2 -a**2



r = sqrt((w+sqrt(w**2 + 4.0_dp*a**2*z**2))/(2.0_dp))
theta = acos(z/r)
phi = atan2(y,x)


sigma = r**2.0_dp +(a*cos(theta))**2.0_dp
mm = sqrt(r**2 + a**2)



rdot0 = mm*r*sin(theta)*cos(phi)*xdot/sigma &
       +mm*r*sin(theta)*sin(phi)*ydot/sigma &
       +mm**2*cos(theta)*zdot/sigma



thetadot0 = (mm*cos(theta)*cos(phi) * xdot &
           +mm*cos(theta)*sin(phi) * ydot &
           -r*sin(theta)* zdot&
           )/sigma


phidot0 = (-sin(phi)*xdot + cos(phi)*ydot)/(mm*sin(theta))


!Define some useful quantitites
sig = r**2.0_dp +(a*cos(theta))**2.0_dp
del = r**2.0_dp - 2.0_dp*r +a**2



!Calculate ray velocity magnitude
s1 = sig-2.0*r
Eobs = 1.0_dp !2.0_dp*PI*nu_obs*1e9
!Eprime = sqrt( (Eobs**2 ) / (s1*(rdot0**2 /del + thetadot0**2) + del*sin(theta)**2*phidot0**2) )





!... and correct velocity components
!rdot0 = rdot0 * Eprime ; thetadot0 = thetadot0 * Eprime ; phidot0 = phidot0 * Eprime

!...and double check the derived energy matches the one you want
E2 = s1*(rdot0**2.0/del +thetadot0**2.0) + del*sin(theta)**2.0*phidot0**2.0
En  = sqrt(E2)


!Momenta
pr = rdot0 *sig/del
ptheta = sig*thetadot0


!Angular momentum
Lz = ((sig*del*phidot0 - 2.0*a*r*En) * sin(theta)**2.0/s1)

!Normalize to E=1
pr = pr/Eobs
ptheta = ptheta/Eobs
Lz = Lz/Eobs




p(1) = 1.0_dp
p(2) = pr
p(3) = ptheta
p(4) = Lz



end subroutine create_ray







subroutine covariant_transform(v,vprime,x,coords)
!Caclulates the coordinate transform between covariatn vectors
!Takes input vector v, output vprime
! 3-position x and input coordinate, coords

!Arguments
real(kind=dp), dimension(4),intent(in) :: v
real(kind=dp), dimension(4),intent(out) :: vprime
real(kind=dp), dimension(3),intent(in) :: x
character(len=2) :: coords !Cartesian CT or Boyer Lindquist BL

!Other 
real(kind=dp) :: m,r,theta,phi

r = x(1) ; theta = x(2) ; phi = x(3)
m = sqrt(r**2 + a**2)




vprime(1) = v(1) !Time component

if (coords .EQ. 'CT') then

vprime(2) = r/m *sin(theta)*cos(phi) * v(2) + r/m * sin(theta)*sin(phi) * v(3) + cos(theta) * v(4)

vprime(3) = m*cos(theta)*cos(phi) * v(2) + m*cos(theta)*sin(phi) * v(3) - r*cos(theta) * v(4)

vprime(4) = -m*sin(theta)*sin(phi) * v(2) + m*sin(theta)*cos(phi) * v(3)


endif




end subroutine covariant_transform








subroutine transform_to_radiation_point(x,s,tau,xi,xi_B)

!Returns the centre of maxx xi, and the radiation point xi_B in carteisan coordinates. The normal is then the differnece


!arguments
real(kind=dp), intent(in), dimension(4) :: x,s !BL coordinates going in
real(kind=dp), intent(in) :: tau
real(kind=dp), intent(inout), dimension(3) :: xi, xi_B

!other
!Data to load
real(kind=dp) :: r,theta,phi,s1,s2,s3
!Cartesian components of spin
real(kind=dp) :: sx,sy,sz
!Spin angles
real(kind=dp) :: thetaSPIN, phiSPIN
!Chi - beam phase
real(kind=dp) :: chi
!Rotation matrices
real(kind=dp), dimension(3,3) :: Rz, Ry
!Radiation beam vector
real(kind=dp), dimension(3) :: Bvector
!Cartesian components of COM and B Vector
real(kind=dp) :: m


!Load the data from the vectors
r = x(2) ; theta = x(3) ; phi = x(4)
s1 = s(2) ; s2 = s(3) ; s3 = s(4)


!Transform to cartesian
Sx = S1*sin(theta)*cos(phi) + S2*r*cos(theta)*cos(phi) - S3*r*sin(theta)*sin(phi)
Sy = S1*sin(theta)*sin(phi) + S2*r*cos(theta)*sin(phi) + S3*r*sin(theta)*cos(phi)
Sz = S1*cos(theta) - S2*r*sin(theta)



!Get the angles
thetaSPIN = atan2(sqrt(Sx**2 + Sy**2),Sz)
phiSPIN = atan2(Sy,Sx)

!Get the beam phase
chi = 2.0_dp*PI*tau/p0

!Construct the rotation matrices

Rz(1,1) = cos(phiSPIN) ; Rz(1,2) = -sin(phiSPIN) ; Rz(1,3) = 0.0_dp
Rz(2,1) = sin(phiSPIN) ; Rz(2,2) = cos(phiSPIN) ;  Rz(2,3) = 0.0_dp
Rz(3,1) = 0.0_dp ;       Rz(3,2) =0.0_dp ;         Rz(3,3) = 1.0_dp



Ry(1,1) = cos(thetaSPIN) ; Ry(1,2) = 0.0_dp ; Ry(1,3) = sin(thetaSPIN)
Ry(2,1) = 0.0_dp ;         Ry(2,2) = 1.0_dp ; Ry(2,3) = 0.0_dp
Ry(3,1) = -sin(thetaSPIN); Ry(3,2) =0.0_dp  ; Ry(3,3) = cos(thetaSPIN)


!Define the vectors in cartesian cpts
Bvector(1) = sin(psi)*cos(chi) ; Bvector(2) =sin(psi)*sin(chi) ; Bvector(3) = cos(psi)
m = sqrt(r**2 +a**2)
xi(1) = m*sin(theta)*cos(phi) ; xi(2) = m*sin(theta)*sin(phi) ; xi(3) = r*cos(theta)

!Perform the transformations


xi_B = RPSR*1e3*convert_m*Bvector
!xi_B = Bvector

if (stheta .NE. 0.0_dp)then

    xi_B = matmul(Ry,xi_B)
    xi_B = matmul(Rz,xi_B) + xi

else

xi_B = xi + xi_B


endif



end subroutine transform_to_radiation_point



subroutine PostKeplerianDelays(r,t,ur,&
                               delta_E,delta_R,delta_S)
!Arguments
real(kind=dp) :: r,t,delta_E,delta_R,ur,delta_S

!Other
real(kind=dp) :: cosE, sinE
real(kind=dp) :: gBAR, pBAR, aBAR, tBAR
real(kind=dp) :: alpha, beta

!Calculate the eccentric anomaly
cosE = (1.0_dp - r/semi_major)/eccentricity
sinE = sign(sqrt(1.0_dp -cosE**2),ur)



!And get the Einstein delay.
!Assumes r0 = rp
delta_E = gam * (sinE) +1.50_dp * t / semi_major



!Get the Roemer delay
!Takes sin(I) = 1
alpha = semi_major*1.0_dp*sin(phi_init)
beta = sqrt(1-eccentricity**2) * semi_major*1.0_dp*cos(phi_init)

delta_R = alpha*(cosE - eccentricity) + beta*sinE






end subroutine PostKeplerianDelays




subroutine transform_to_comoving_frame(vector,u,x,vector_comoving,vector_type)

!Takes a covariant vector and a contravariant velocity

!Arguments
real(kind=dp), intent(IN), dimension(4) :: vector,u,x !vector, 4-velocity, position
real(kind=dp), intent(OUT), dimension(4) :: vector_comoving
character(len=*), intent(in) :: vector_type


!Other
real(kind=dp), dimension(4) :: vector_contra, vector_covar
real(kind=dp), dimension(4,4) :: metric,metricCONTRA !covariant and contravariant metric components
real(kind=dp),dimension(4) :: u_covar
real(kind=dp) :: r, theta,phi,u0,u1,u2,u3,u0_covar,u1_covar,u2_covar,u3_covar
real(kind=dp) :: delta, grr,gthth, N1, N2, N3, pmag, mag1, mag2
real(kind=dp), dimension(4,4) :: transformation_matrix_inverse,transformation_matrix
real(kind=dp) :: magnitude_before, magnitude_after
integer(kind=dp) :: i

!Calculate the metric components
call calculate_covariant_metric(x(2), x(3), metric)
call calculate_contravariant_metric(x(2),x(3), metricCONTRA)


!Transform the contravariant components of the MSP 4-velocity
u_covar = MATMUL(metric,u)

!Check the 4 velocity
!print *, '4 velocity is =', u
!print *, '4-velocity magnitude =', u(1)*u_covar(1) + u(2)*u_covar(2)+u(3)*u_covar(3) +u(4)*u_covar(4)

!Extract variables from the vector format
r = x(2) ; theta = x(3); phi = x(4)
u0 = u(1) ; u1=u(2) ; u2 = u(3) ; u3 = u(4)
u0_covar = u_covar(1) ; u1_covar=u_covar(2) ; u2_covar = u_covar(3) ; u3_covar = u_covar(4)
grr = metric(2,2) ; gthth = metric(3,3)



!Construct transformation matrix
delta = r**2 + a**2 - 2.0_dp*r
N1 = sqrt(-grr * (u0_covar * u0 + u3_covar*u3)*(1.0_dp + u2_covar*u2))
N2 = sqrt(gthth*(1.0_dp + u2_covar*u2))
N3 = sqrt(-(u0_covar*u0 + u3_covar*u3)*delta*sin(theta)**2)



transformation_matrix_inverse(1,1) = -u0_covar
transformation_matrix_inverse(1,2) = -u1_covar
transformation_matrix_inverse(1,3) = -u2_covar
transformation_matrix_inverse(1,4) = -u3_covar

transformation_matrix_inverse(2,1) = u1_covar*u0_covar/N1
transformation_matrix_inverse(2,2) = -grr*(u0_covar*u0 + u3_covar*u3)/N1
transformation_matrix_inverse(2,3) = 0.0_dp
transformation_matrix_inverse(2,4) = u1_covar*u3_covar/N1

transformation_matrix_inverse(3,1) = u2_covar*u0_covar/N2
transformation_matrix_inverse(3,2) = u2_covar*u1_covar/N2
transformation_matrix_inverse(3,3) = gthth*(1.0_dp + u2_covar*u2)/N2
transformation_matrix_inverse(3,4) = u2_covar*u3_covar/N2

transformation_matrix_inverse(4,1) = -delta*sin(theta)**2*u3/N3
transformation_matrix_inverse(4,2) = 0.0_dp
transformation_matrix_inverse(4,3) = 0.0_dp
transformation_matrix_inverse(4,4) = delta*sin(theta)**2*u0/N3



transformation_matrix(1,1) = u0
transformation_matrix(1,2) = u1
transformation_matrix(1,3) = u2
transformation_matrix(1,4) = u3

transformation_matrix(2,1) = (u1_covar*u0) / N1
transformation_matrix(2,2) = (-(u0_covar * u0 + u3_covar * u3)) / N1
transformation_matrix(2,3) = (0.0_dp)/N1
transformation_matrix(2,4) = (u1_covar * u3)/N1

transformation_matrix(3,1) = (u2_covar*u0) / N2
transformation_matrix(3,2) = (u2_covar*u1) / N2
transformation_matrix(3,3) = (1.0_dp + u2_covar*u2)/N2
transformation_matrix(3,4) = (u2_covar*u3)/N2

transformation_matrix(4,1) = (u3_covar) / N3
transformation_matrix(4,2) = (0.0_dp) / N3
transformation_matrix(4,3) = (0.0_dp)/N3
transformation_matrix(4,4) = (-u0_covar)/N3



if (vector_type .EQ. 'covariant') then

vector_contra = matmul(metricCONTRA,vector)
magnitude_before =vector(1)*vector_contra(1) + vector(2)*vector_contra(2) + vector(3)*vector_contra(3) + vector(4)*vector_contra(4)

vector_comoving = matmul(transformation_matrix, vector)



else if (vector_type .EQ. 'contravariant') then

vector_covar = matmul(metric,vector)
magnitude_before =vector(1)*vector_covar(1) + vector(2)*vector_covar(2) + vector(3)*vector_covar(3) + vector(4)*vector_covar(4)


vector_comoving = matmul(transformation_matrix_inverse, vector)

else

print *, 'Vector type ', vector_type, ' is not recognized'

endif





magnitude_after = -vector_comoving(1)**2 + vector_comoving(2)**2 + vector_comoving(3)**2 + vector_comoving(4)**2


!print *, magnitude_before
!print *, magnitude_after



end subroutine transform_to_comoving_frame



subroutine pitch_angle(n,k,pitch)

!Takes vector normal n, tangent ray k

!Arguments
real(kind=dp), intent(in), dimension(3) :: n,k
real(kind=dp),intent(out) :: pitch
!Other
real(kind=dp) :: magN, magK


magN = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
magK = sqrt(k(1)**2 + k(2)**2 + k(3)**2)


pitch = acos((n(1)*k(1) + n(2)*k(2) +n(3)*k(3))/(magN*magK) )



end subroutine pitch_angle




end module analysis
