SUBROUTINE output(rc,nrcycl,eta,En,Pressure)

	USE constant
	USE menu_parameters, ONLY : closure, potential
	USE main_variables
	USE state_variables

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: rc, nrcycl(200)
	DOUBLE PRECISION, INTENT(IN) :: eta(rc)
	DOUBLE PRECISION, INTENT(OUT) :: En, Pressure
	
	INTEGER :: i, j, jj
	DOUBLE PRECISION :: h(1:3,0:n), g(1:3,0:n), coeff, cj, cjj
	
	!OPEN(101, file='gamma1.txt')
   ! OPEN(3,file='input.txt')
   OPEN (4,file='c.txt')
    
	WRITE(2,'(/,3X,"*********** SUMMARY OF PROCESS ***********")')
	
	WRITE(2,'(/,8X,"REFINEMENT",2X,"NEWT-RAPHN",/,9X," CYCLE  ",4X," CYCLE  ",9X,"ETA",/)')

	DO i=1,rc
	  WRITE(2,'(10X,I3,9X,I3,5X,G15.6)') i, nrcycl(i), eta(i)
	END DO

	WRITE(2,'(/,5X,"*********** END OF SUMMARY ***********")')

	WRITE(2,'(//,1X,"Values obtained for the correlation functions")')
	
	WRITE(2,'(/,5X,"For: rhostar = ",F6.4," at Tstar = ",F6.3)') rho, tstar

	WRITE(2,'(/,3X,"i",5X,"r(i)",9X,"gamma(i)",12X,"c(i)",15X,"h(i)",15X,"g(i)",/)')

	j = NINT(1.0/dr)
	jj = NINT(1.0/dr)
    g = 1.d0
    
	  SELECT CASE (closure)
	    CASE(1)
	!	  DO i=0,180
	!        c(i) = (1.0 + gamma1(i))*f(i)
	!		h(i) = gamma1(i) + c(i)
	!		g(i) = h(i) + 1.0
	!		WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(i), c(i), h(i), g(i)
	!	  END DO
		CASE(2)
		  DO j=1,3
		  DO i=0,n/2
		    c(j,i) = (f(j,i) + 1.0)*DEXP(gamma1(j,i))-gamma1(j,i)-1.0
		    
		 !   c(j,i) = EXP(f(j,i)+gamma1(j,i)) - gamma1(j,i) - 1.0
			h(j,i) = gamma1(j,i) + c(j,i)
			g(j,i) = h(j,i) + 1.0
			WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(j,i), c(j,i), h(j,i), g(j,i)
		  END DO
		    WRITE(2,*) j
		  END DO
    
	    DO i = 0 , n/2
	  	WRITE(4,*) gamma1(1,i), gamma1(2,i), gamma1(3,i)
	    END DO

		!   DO i = 0, 255
		!   write(3, *) gamma1(1,i), gamma1(2,i), gamma1(3,i)
		!   END DO
		     
		CASE(3)
	!	  DO i=0,j
	!		c(i) = (1.0 + gamma1(i))*f(i)
	!		h(i) = gamma1(i) + c(i)
	!		g(i) = h(i) + 1.0
	!		WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(i), c(i), h(i), g(i)
	!	  END DO
	!	  DO i=j+1,180
	!		c(i) = LOG(f(i)+1.0)
	!		h(i) = gamma1(i) + c(i)
	!		g(i) = h(i) + 1.0
	!		WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(i), c(i), h(i), g(i)
	!	  END DO
		CASE(4)
	!     i=0
	!	  c(0) = (1.0 + gamma1(0))*f(0)
	!	  h(0) = gamma1(0) + c(0)
	!	  g(0) = h(0) + 1.0
	!	  WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(i), c(i), h(i), g(i)
	!	  DO i=1,100
	!	    c(i) = (f(i)+1.0)*(1.0+(EXP(bridge(i)*gamma1(i))-1.0)/bridge(i))-gamma1(i)-1.0
	!		h(i) = gamma1(i) + c(i)
	!		g(i) = h(i) + 1.0
	!		WRITE(2,'(1X,I3,4X,F5.3,4(4X,G15.6))') i, r(i), gamma1(i), c(i), h(i), g(i)
	!	  END DO
	  END SELECT
	
	coeff = 4.0*pi
!	cj = c(j)
!	cjj = c(jj)

!	IF (potential/=1) THEN
!	  c(j) = (c(j+1) + (c(j+1)-c(j+2)) + cj)/2.0
!	  c(jj) = cjj/2.0
!	END IF
!	CALL FFT(n,coeff,dr,r,k,c,ftc)

!	c(j) = cj
!	c(jj) = cjj

!	CALL FFT(n,coeff,dr,r,k,gamma1,ftg)
	
	CALL ENERGY(g, En , Pressure)
	
	

!	WRITE(2,'(//,1X,"Fourier transform of the direct correlation function c(k), rho*c(k) and gamma(k)",/)')
!	WRITE(2,'(1X,I3,4X,F6.2,4X,3G15.6)') (i, k(i), ftc(i), rho*ftc(i), ftg(i), i=0,n)

END SUBROUTINE output  