MODULE constant

	IMPLICIT NONE

	DOUBLE PRECISION :: pi = 3.1415927, lambda = 1.69925

END MODULE constant

MODULE menu_parameters

	IMPLICIT NONE

	INTEGER :: model, potential, closure, iteration

END MODULE menu_parameters

MODULE state_variables

	IMPLICIT NONE

	DOUBLE PRECISION :: rho, tstar
	DOUBLE PRECISION :: rho1, rho2

END MODULE state_variables

MODULE main_variables

	IMPLICIT NONE

	INTEGER :: numnod, n
	PARAMETER (numnod=15,n=1023)
	DOUBLE PRECISION ::	dr, dk, r(0:n), k(0:n)
	DOUBLE PRECISION :: f(1:3,0:n), gamma1(1:3,0:n), c(1:3,0:n), bridge(1:3,0:n)   ! 2 component ionic system nofcomp = 3
	DOUBLE PRECISION :: ftc(1:3,0:n), ftg(1:3,0:n)
	DOUBLE PRECISION :: z1 = 1.0, z2 = -1.0
	DOUBLE PRECISION :: dielectric = 78.15
	DOUBLE PRECISION :: charge = 1.609

END MODULE main_variables

!-------------------------------------------------------------------

PROGRAM Integral_equation

! For Compaq Visual Fortran Compiler, to use IMSL library, use the following statement in 
! all subroutines using IMSL library:
    

!	USE MSIMSL 

! or if only one IMSL subroutine is used, for example LINDS:

!    USE MSIMSL, ONLY : LINDS
    
! For Intel Fortran Compiler, to use IMSL library, use the following statements in the main program only:
    
    INCLUDE 'link_fnl_static.h'
    !DEC$ OBJCOMMENT LIB:'libiomp5md.lib' 

!--------------------------------------------------------------------------------

	USE constant
	USE menu_parameters
	USE state_variables
	USE main_variables
	
	IMPLICIT NONE
	
	INTEGER :: i, rc, nrc, limit, nrcycl(200), j
	PARAMETER (limit=30)

	DOUBLE PRECISION :: P(1:3,1:numnod,0:n), Q(1:3,1:numnod,0:n), a(1:3,1:numnod)
	DOUBLE PRECISION :: delta_gamma(1:3,0:n), dgamma_prime(1:3,0:n), d(1:3,1:numnod), gamma_prime1(1:3,0:n)
	DOUBLE PRECISION :: a_prime(1:3,1:numnod), J_inv(1:numnod,1:numnod)
	DOUBLE PRECISION :: eta(200)
	DOUBLE PRECISION :: hel_A(0:500),pressure, hel, en
	DOUBLE PRECISION :: T_star(0:500), delta_T1, delta_T2, dr1

	LOGICAL :: converge

	INTEGER :: nn, maxfcn, excess, kk
    integer :: clck_counts_beg, clck_counts_end, clck_rate
    
    
    call system_clock ( clck_counts_beg, clck_rate )

	OPEN(3,file='input.txt')
	OPEN(2,file='results.txt')
	OPEN(10, file='Hel.txt')

	!CALL menu
	model = 1
	potential = 3
	closure = 2
	iteration = 2
	excess = 0
		
	PRINT *, 'Input values of rho* and T*: '
!	READ *, rho, tstar
    
    rho = 0.35   !packing fraction < 0.74
 !   IF (excess==0) THEN
    tstar = 0.59441
 !   END IF
        
	dr = 0.009
    dr1 = dr
    
	rho = rho * 6.0 /pi         ! reduced density

	rho1 = rho/2.0;
	rho2 = rho/2.0;
	
	IF (iteration == 2) THEN
	  gamma1 = 0.0
	  DO i=0,n/2
	    READ(3,*) gamma1(1,i), gamma1(2,i), gamma1(3,i)
	  END DO
	END IF
	
!	OPEN (UNIT=3, FILE="input.txt", STATUS="OLD")  ! Current directoryCLOSE (UNIT=5, STATUS="DELETE")
 !   CLOSE (UNIT=3, STATUS="DELETE")
    
    IF (excess==0) THEN
    CALL initialize(delta_gamma)
	CALL PQ_generator(n,numnod,limit,P,Q)
	CALL decompose_gamma(numnod,limit,n,Q,P,gamma1,a,delta_gamma)

! Begin main loop

	 DO rc=1,200

	   DO nrc=1,200
	    
	    nrcycl(rc)=nrc

		CALL new_gamma(n,numnod,P,a,delta_gamma,gamma1)
		CALL elementary_cycle(gamma_prime1)
		CALL decompose_gamma(numnod,limit,n,Q,P,gamma_prime1,a_prime,dgamma_prime)

! Test for Newton-Raphson cycle convergence
		CALL NR_test(numnod,a,a_prime,d,converge)
		IF (converge) EXIT
		CALL Jacobian(1,limit,Q,P,J_inv)
		CALL new_a(1,numnod,J_inv,d,a)
		CALL Jacobian(2,limit,Q,P,J_inv)
		CALL new_a(2,numnod,J_inv,d,a)
		CALL Jacobian(3,limit,Q,P,J_inv)
		CALL new_a(3,numnod,J_inv,d,a)						
	  END DO

! Test for refinement cycle convergence

	  CALL RC_test(n,rc,dr,gamma1,gamma_prime1,eta,converge)
	  IF(converge .AND. nrc<=200) EXIT
          
      DO j=1,3
      DO i=0,n
	    delta_gamma(j,i) = dgamma_prime(j,i)
	  END DO
	  END DO
	END DO
	
	CALL output(rc,nrcycl,eta, En, Pressure)
	write(*,*) En, Pressure
	
	END IF
	  
	!----------------calculate Hel energy by thermodynamic integration-------
	delta_T1 = 0.1  
	delta_T2 = 1
	!delta_T3 = 0.001
	IF (excess ==1) THEN
	
	 DO kk = 1,70
	   IF (kk<=30) THEN
	   T_star(kk) = tstar + delta_T1*(kk-1)
	   ELSE
	   T_star(kk) = T_star(kk-1) + delta_T2
	   END IF
	   
	!   write(4,*) T_star(kk)
	 END DO
	 Hel_A = 0;
	  
	 DO kk = 1 , 70
	 
	 tstar = T_star(kk)
	 
	 if (tstar > 1 .and. tstar <= 5) then
	   dr = dr1 * 1.0
	 end if
	 if (tstar >5  .and. tstar <=15) then
	   dr = dr1*1.0
	 end if
	 if (tstar >15 ) then
	   dr = dr1*1.0
	 end if 
	  
	 	
	 CALL initialize(delta_gamma)
	 CALL PQ_generator(n,numnod,limit,P,Q)
	 CALL decompose_gamma(numnod,limit,n,Q,P,gamma1,a,delta_gamma)

! Begin main loop

	  DO rc=1,200

	    DO nrc=1,200
	    
	    nrcycl(rc)=nrc
		CALL new_gamma(n,numnod,P,a,delta_gamma,gamma1)
		CALL elementary_cycle(gamma_prime1)
		CALL decompose_gamma(numnod,limit,n,Q,P,gamma_prime1,a_prime,dgamma_prime)

! Test for Newton-Raphson cycle convergence
		CALL NR_test(numnod,a,a_prime,d,converge)
		IF (converge) EXIT
		CALL Jacobian(1,limit,Q,P,J_inv)
		CALL new_a(1,numnod,J_inv,d,a)
		CALL Jacobian(2,limit,Q,P,J_inv)
		CALL new_a(2,numnod,J_inv,d,a)
		CALL Jacobian(3,limit,Q,P,J_inv)
		CALL new_a(3,numnod,J_inv,d,a)						
	  END DO

! Test for refinement cycle convergence

	  CALL RC_test(n,rc,dr,gamma1,gamma_prime1,eta,converge)
	  IF(converge .AND. nrc<=200) EXIT          
      DO j=1,3
         DO i=0,n
		   delta_gamma(j,i) = dgamma_prime(j,i)
		 END DO
	   END DO
	
	END DO
		
	CALL output(rc,nrcycl,eta, En, Pressure)
	Hel_A(kk) = En
	write(*,*) 'number of T T=',kk, tstar 
	write(10,*) Hel_A(KK)
    
    END DO
    
    Hel = 0;
    DO i = 1, 69
       IF (i<30) THEN
       Hel = Hel + (Hel_A(i)/T_star(i) + Hel_A(i+1)/T_star(i+1))*(delta_T1)/2.0  ! Trapezio rule
       ELSE
       Hel = Hel + (Hel_A(i)/T_star(i) + Hel_A(i+1)/T_star(i+1))*(delta_T2)/2.0
       END IF
    END DO

    write(*,*) 'Hel =', HEL
    
    ENDIF
    
    call system_clock ( clck_counts_end, clck_rate )
    write (*, *) 'clock time =', (clck_counts_end - clck_counts_beg) / real (clck_rate)

END PROGRAM Integral_equation
	
	 