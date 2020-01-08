SUBROUTINE Jacobian(S,limit,Q,P,J_inv)

!	USE MSIMSL, ONLY : LINRG
	USE constant
	USE menu_parameters
	USE main_variables
	USE state_variables

	IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: S
	INTEGER, INTENT(IN) :: limit
	DOUBLE PRECISION, INTENT(IN) :: P(1:3,1:numnod,0:n), Q(1:3,1:numnod,0:n)
	DOUBLE PRECISION, INTENT(OUT) :: J_inv(1:numnod,1:numnod)

	DOUBLE PRECISION :: dgpdg(0:limit,0:limit), J(1:numnod,1:numnod)
	DOUBLE PRECISION :: zz(0:511), dd(0:70), m1,m2,m3, sum, delta
	DOUBLE PRECISION :: jacobian_m(1:numnod, 1:numnod), jacobian_m_inv(1:numnod, 1:numnod)
	

	INTEGER :: i, j1, mm, l, alpha, beta1, jj

	jj = NINT(1.0/dr)

	SELECT CASE (S)

	  CASE (1) ! Particle-particle OZ
	
		  DO i=0,(n+1)/2-1 
	 !     DO i=0,n/2-1
		  m1 = (rho1*ftc(1,i) - 1.0)*(ftc(2,i)*rho2 - 1.0) - ftc(3,i)**2*rho1*rho2
		  m2 = rho1*rho2*ftc(3,i)**2 - 2.0*ftc(1,i)*rho1*(ftc(2,i)*rho2 - 1.0)
		  m3 = rho1*(ftc(2,i)*rho2-1.0)*(ftc(1,i)**2*rho1*(ftc(2,i)*rho2 - 1.0) - ftc(3,i)**2*rho2*(ftc(1,i)*rho1 + 1.0) )
		  zz(i) = m2/m1 + m3/m1**2
		  END DO
		
	 CASE(2)
	    
	    DO i=0,(n+1)/2-1
	!    DO i=0,n/2-1
	    m1 = (rho1*ftc(1,i) - 1.0)*(ftc(2,i)*rho2 - 1.0) - ftc(3,i)**2*rho1*rho2
	    m2 = rho1*rho2*ftc(3,i)**2 - 2.0*ftc(2,i)*rho2*(ftc(1,i)*rho1 - 1.0)
	    m3 = rho2*(ftc(1,i)*rho1-1.0)*(ftc(2,i)**2*rho2*(ftc(1,i)*rho1 - 1.0) - ftc(3,i)**2*rho1*(ftc(2,i)*rho2 + 1.0) )
	    zz(i) = m2/m1 + m3/m1**2	    
	    END DO
	    
	  CASE(3)
	  
	     DO i=0,(n+1)/2-1
	     
	!     DO i=0,n/2-1
	     m1 = (ftc(1,i)*rho1-1.0)*(ftc(2,i)*rho2 - 1.0) - ftc(3,i)**2*rho1*rho2;
	     m2 = 2.0*ftc(3,i)**2*rho1*rho2
	     zz(i) = 1.0/m1 + m2/m1**2 - 1.0
	     END DO

	END SELECT


! Calculate d_gamma_prime_d_gamma(0,j)

	DO j1=0,limit
	  sum = 0.0
	  DO mm=0,(n+1)/2-1
!      DO mm=0,n/2-1
	    sum = sum + k(mm)*zz(mm)*SIN(k(mm)*r(j1))
	  END DO
	  SELECT CASE (closure)
	    CASE(1)  ! PY closure
!   	  dgpdg(0,j1) = 2.0*dr*r(j1)/pi*f(j1)*dk*sum
		CASE(2)  ! HNC closure
    	!  dgpdg(0,j1) = 2.0*dr*r(j1)/pi*(EXP(f(S,j1)+gamma1(S,j1)) -1.0)*dk*sum
    	   dgpdg(0,j1) = 2.0*dr*r(j1)/pi*((f(S,j1) + 1.0)*DEXP(gamma1(S,j1))-1.0)*dk*sum
    	CASE(3)  ! MSA closure
!		  IF (j1<=jj) THEN
!   	    dgpdg(0,j1) = 2.0*dr*r(j1)/pi*f(j1)*dk*sum
!		  ELSE
!			dgpdg(0,j1) = 0.0
!		  END IF
		CASE(4)
!   	  dgpdg(0,j1) = 2.0*dr*r(j1)/pi*((f(j1)+1.0)*EXP(bridge(j1)*gamma1(j1))-1.0)*dk*sum
	  END SELECT
	END DO

! Calculate the remaining d_gamma_prime_d_gamma(i,j)

	DO l=0,2*limit+1
	  sum = 0.0
	  DO mm=0,(n+1)/2-1
!      DO mm=0,n/2-1
	    sum = sum + zz(mm)*COS(k(mm)*r(l))
	  END DO
	  dd(l) = dk*sum
	END DO

	DO i=1,limit
	  DO j1=0,limit
	    l = ABS(i-j1)
		mm = i+j1
	    SELECT CASE (closure)
	      CASE(1)  ! PY closure
	!	    dgpdg(i,j1) = dr*r(j1)/(pi*r(i))*f(j1)*(dd(l)-dd(mm))
		  CASE(2)  ! HNC closure
		!    dgpdg(i,j1) = dr*r(j1)/(pi*r(i))*( exp(f(S,j1)+gamma1(S,j1))-1.0 )*(dd(l)-dd(mm))
		     dgpdg(i,j1) = dr*r(j1)/(pi*r(i))*( (f(S,j1) + 1.0)*DEXP(gamma1(S,j1))-1.0 )*(dd(l)-dd(mm))		
		  CASE(3)  ! MSA closure
	!	    IF (j1<=jj) THEN
	!	      dgpdg(i,j1) = dr*r(j1)/(pi*r(i))*f(j1)*(dd(l)-dd(mm))
	!	    ELSE
	!		  dgpdg(i,j1) = 0.0
	!	    END IF
		  CASE(4)  ! TC closure
	!	    dgpdg(i,j1) = dr*r(j1)/(pi*r(i))*&
	!				&((f(j1)+1.0)*EXP(bridge(j1)*gamma1(j1))-1.0)*(dd(l)-dd(mm))
	    END SELECT
	  END DO
	END DO

! Jacobian matrix
 
	DO alpha=1,numnod
	  DO beta1=1,numnod
	    delta = 0.0
		IF(alpha == beta1) delta = 1.0
		sum = 0.0
		DO i=1,limit
		  DO j1=1,limit
		    sum = sum + Q(S,alpha,i)*dgpdg(i,j1)*P(S,beta1,j1)
		  END DO
		END DO
		J(alpha,beta1) = delta - sum
	  END DO
	END DO

   
   ! Jacobian_m =0; Jacobian_m_inv =0;


   CALL DLINRG(numnod,J,numnod,J_inv,numnod)



END SUBROUTINE Jacobian