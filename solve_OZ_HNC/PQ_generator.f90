SUBROUTINE PQ_generator(n,numnod,limit,P,Q)

!	USE MSIMSL, ONLY : LINDS
	USE constant
	USE menu_parameters, ONLY : iteration

	IMPLICIT NONE

	INTEGER, INTENT(IN)	:: n, numnod, limit
	DOUBLE PRECISION, INTENT(OUT) :: P(1:3,1:numnod,0:n),Q(1:3,1:numnod,0:n)
	
	DOUBLE PRECISION :: rr(1:3,1:numnod,1:numnod),bb(1:3,1:numnod,1:numnod)
	DOUBLE PRECISION :: ia(1:3,0:numnod), sum
	DOUBLE PRECISION :: R_R(1:numnod,1:numnod), B_B(1:numnod,1:numnod)
	
	INTEGER	:: alpha, beta1, i, j, k
        
    ia = 0; R_R = 0 ; B_B =0;  rr = 0; bb =0;
	DO alpha=1,numnod
	  ia(1,alpha) = limit/numnod*alpha 
	  ia(2,alpha) = limit/numnod*alpha 
	  ia(3,alpha) = limit/numnod*alpha  
	END DO
		 
! Initialize base functions and conjugates to zero

	Q = 0.0
	P = 0.0

	SELECT CASE (iteration)

!	  CASE(1)

! Calculate P(alpha,i)

!		DO alpha=1,numnod
!		  DO i=0,limit
!			P(alpha,i) = SIN(pi/(n+1)*i*alpha)
!		  END DO
!		END DO
!		P = P/SQRT(REAL(n+1))
!		Q = P

	  CASE(2)

! Initialize R & B matrices to zero

		rr = 0.0
		bb = 0.0

! Calculate P(1,i)
		
		DO j = 1,3
		i = 0
		DO WHILE(i <= ia(j,1))
		  P(j,1,i) = (ia(j,1)-i)/ia(j,1)
		  i = i + 1
		END DO
		END DO

! Calculate P(alpha,i) for alpha>=2
        DO i = 1, 3
		DO alpha=2,numnod
		  CALL P_generator(i,n,numnod,alpha,ia(i,alpha-2),ia(i,alpha-1),ia(i,alpha),limit,P)
		END DO
		END DO

! Calculate R matrix
       DO j = 1, 3
		DO alpha=1,numnod
		  DO beta1=1,numnod
			sum = 0.0
			DO i=0,limit
			  sum = sum + P(j,alpha,i)*P(j,beta1,i)
			END DO
			rr(j,alpha,beta1) = sum
		  END DO
		END DO
	  END DO

! Calculate B matrix
     DO k = 1, 3
       DO i = 1,numnod
         DO j = 1, numnod
            R_R(i,j) = rr(k,i,j)
         ENDDO
       ENDDO
	   CALL DLINDS(numnod,R_R,numnod,B_B,numnod)
	   DO i = 1,numnod
         DO j = 1, numnod
            bb(k,i,j) = B_B(i,j)
         ENDDO
       ENDDO
       R_R = 0; B_B=0
     END DO

! Calculate base function conjugates Q(alpha,i)
       DO j = 1, 3
		DO alpha=1,numnod
		  DO i=0,limit
			sum = 0.0
			DO beta1=1,numnod
			  sum = sum + bb(j,alpha,beta1)*P(j,beta1,i)
			END DO
			Q(j,alpha,i) = sum
		  END DO
		END DO
	 END DO

	END SELECT

END SUBROUTINE PQ_generator


SUBROUTINE P_generator(k,n,numnod,alpha,lo,med,hi,limit,P)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, numnod, alpha, limit, k
	DOUBLE PRECISION, INTENT(IN) :: lo, med, hi
	DOUBLE PRECISION, INTENT(INOUT) :: P(1:3,numnod,0:n)

	INTEGER :: i

	DO i=0,limit
	  IF((i<=lo).OR.(i>hi)) P(k,alpha,i) = 0.0
	  IF((i>lo).AND.(i<=med)) P(k,alpha,i) = (i-lo)/(med-lo)
	  IF((i>med).AND.(i<=hi)) P(k,alpha,i) = (hi-i)/(hi-med)
	END DO

END SUBROUTINE P_generator


