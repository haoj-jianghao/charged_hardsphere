SUBROUTINE initialize(delta_gamma)
	
	USE constant
	USE menu_parameters
	USE main_variables
	USE state_variables
	
	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(OUT) :: delta_gamma(0:n)
	
	DOUBLE PRECISION :: s, s6, u1,u2,u3, alpha
!	REAL :: ff(0:n)

	INTEGER :: i, j, jj

	delta_gamma = 0.0
	j = ANINT(1.0/dr)
	jj = ANINT(lambda/dr)

! Calculate reduced r(i) 

	DO i=0,n
	  r(i) = i*dr
	END DO

! Calculate Mayer function f(i) using potential model

	SELECT CASE (potential)

	  CASE (1)	! Lennard-Jones
	
!		f(:14) = -1.0
!		DO i=15,n
!		  s = 1.0/r(i)
!		  s6 = s**6
!		  u = 4.0/tstar*s6*(s6-1.0)
!		  f(i) = EXP(-u)-1.0
!		END DO

	  CASE (2)  ! Square-Well

 !		f = 0.0
!		f(:j) = -1.0
!		DO i=j+1,jj
!		  f(i) = EXP(1.0/tstar)-1.0
!		END DO

	  CASE (3)  ! Hard-Sphere

 		f = 0.0      

		DO i = 0, j
 		f(1,i) = -1
 		f(2,i) = -1
 		f(3,i) = -1
 		END DO
 		
 		DO i = j+1, n
 !      DO i = 0, n
 		u1 = z1*z1/r(i)/tstar;
 		u2 = z2*z2/r(i)/tstar;
 		u3 = z1*z2/r(i)/tstar;  ! u -- -1.0*beta*potential
		f(1,i) = DEXP(-u1) - 1.0
		f(2,i) = DEXP(-u2) - 1.0
		f(3,i) = DEXP(-u3) - 1.0
		ENDDO
 	
	END SELECT

!	alpha = 1.9
!	bridge = 1.0-EXP(-alpha*r)

! Calculate the Fourier tranform factor k(j)

	dk = 2.0*pi/((n+1)*dr)
	DO j=0,n
	  k(j) = j
	END DO
	k = dk*k


	
END SUBROUTINE initialize