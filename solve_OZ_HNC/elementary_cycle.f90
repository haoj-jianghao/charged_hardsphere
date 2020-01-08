SUBROUTINE elementary_cycle(gamma_prime)

	USE constant
	USE menu_parameters
	USE main_variables
	USE state_variables
	
	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(OUT) :: gamma_prime(1:3,0:n)

	DOUBLE PRECISION :: coeff, cj(1:3,0:n)
	INTEGER :: i, j, jj
	DOUBLE PRECISION :: deltak(0:n)
	DOUBLE PRECISION :: lambda_rho, alpha_l
	 
	j = NINT(1.0/dr)
	jj = NINT(lambda/dr)
    
    lambda_rho = 0.3
    alpha_l = 1.0/lambda_rho
    
    
! Calculate new c(i)

	SELECT CASE (closure)

	  CASE(1)  ! PY closure
	
!		c = (1.0 + gamma1)*f

	  CASE(2)  ! HNC closure
        
         
	    
	   DO i=0,n
	   c(1,i) = (f(1,i) + 1.0)*DEXP(gamma1(1,i))-gamma1(1,i)-1.0   
	   c(2,i) = (f(2,i) + 1.0)*DEXP(gamma1(2,i))-gamma1(2,i)-1.0    
	   c(3,i) = (f(3,i) + 1.0)*DEXP(gamma1(3,i))-gamma1(3,i)-1.0	
	   END DO
	   
	   c(1,j) = (c(1,j+1) + (c(1,j+1)-c(1,j+2)) + c(1,j))/2.0
	   c(2,j) = (c(2,j+1) + (c(2,j+1)-c(2,j+2)) + c(2,j))/2.0
	   c(3,j) = (c(3,j+1) + (c(3,j+1)-c(3,j+2)) + c(3,j))/2.0   	    
       cj = c
	  
	   DO i = 1, n
	   c(1,i) = c(1,i) + z1*z1/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	   c(2,i) = c(2,i) + z2*z2/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	   c(3,i) = c(3,i) + z1*z2/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	   END DO	   
		    

	   

	    
     CASE(3)  ! MSA closure

!		DO i=0,j
!		  c(i) = (1.0 + gamma1(i))*f(i)
!		END DO

!		DO i=j+1,n
!		  c(i) = LOG(f(i)+1.0)
!		END DO

	  CASE(4)  ! TC closure

!	    c(0) = (1.0 + gamma1(0))*f(0)
!		DO i=1,n
!		  c(i) = (f(i)+1.0)*(1.0+(EXP(bridge(i)*gamma1(i))-1.0)/bridge(i))-gamma1(i)-1.0
!		END DO

	END SELECT

! Calculate Fourier transform of c(i)

	  coeff = 4.0*pi

      CALL FFT(n,coeff,dr,r,k,c,ftc)           

	  	  
	  DO i = 1, n
	  ftc(1,i) = ftc(1,i) - 4*pi*z1*z1/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
	  ftc(2,i) = ftc(2,i) - 4*pi*z2*z2/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
      ftc(3,i) = ftc(3,i) - 4*pi*z1*z2/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
      END DO
            
      c = cj   
      
      ftc(1,0) = ftc(1,1)-(ftc(1,2)-ftc(1,1))
	  ftc(2,0) = ftc(2,1)-(ftc(2,2)-ftc(2,1))
	  ftc(3,0) = ftc(3,1)-(ftc(3,2)-ftc(3,1))


! Calculate Fourier transform of gamma_prime(i)
 
	SELECT CASE (model)

	  CASE (1)    ! Particle-particle OZ
	    DO i = 0,n
	    deltak(i) = (1.0 - rho1*ftc(1,i))*(1.0 - rho2*ftc(2,i)) - rho1*rho2*ftc(3,i)*ftc(3,i)	    
		ftg(1,i) = rho1*ftc(1,i)**2*(1.0 - rho2*ftc(2,i)) + rho2*ftc(3,i)**2*(1.0 + rho1*ftc(1,i))
		ftg(2,i) = rho2*ftc(2,i)**2*(1.0 - rho1*ftc(1,i)) + rho1*ftc(3,i)**2*(1.0 + rho2*ftc(2,i))
		ftg(1,i) = ftg(1,i) /deltak(i)
		ftg(2,i) = ftg(2,i) /deltak(i)
		ftg(3,i) = ftc(3,i) * ( 1.0/deltak(i) -1.0)		
		END DO
	 END SELECT

     DO i = 1, n
       ftg(1,i) = ftg(1,i) - 4*pi*z1*z1/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
       ftg(2,i) = ftg(2,i) - 4*pi*z2*z2/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
       ftg(3,i) = ftg(3,i) - 4*pi*z1*z2/tstar*alpha_l**2/k(i)**2/( k(i)**2 + alpha_l**2 )
     END DO

  ! Calculate gamma_prime(i)

	
	coeff = 1.0/(2.0*pi*pi)
	CALL FFT(n,coeff,dk,k,r,ftg,gamma_prime)
	
	DO i = 1, n
	gamma_prime(1,i) = gamma_prime(1,i) + z1*z1/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	gamma_prime(2,i) = gamma_prime(2,i) + z2*z2/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	gamma_prime(3,i) = gamma_prime(3,i) + z1*z2/tstar/r(i)*(1.0 - dexp(-1.0*r(i)*alpha_l))
	END DO
      
END SUBROUTINE elementary_cycle
