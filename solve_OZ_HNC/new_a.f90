SUBROUTINE new_a(S,numnod,J_inv,d,a)

	IMPLICIT NONE
    
    INTEGER, INTENT (IN) :: S
	INTEGER, INTENT(IN) :: numnod
	DOUBLE PRECISION, INTENT(IN) :: J_inv(1:numnod,1:numnod), d(1:3,1:numnod)
	DOUBLE PRECISION, INTENT(INOUT) :: a(1:3,1:numnod)

	INTEGER :: alpha, beta, k
	DOUBLE PRECISION :: sum

	DO alpha=1,numnod
	  sum = 0.0
	  DO beta=1,numnod
	    sum = sum + J_inv(alpha,beta)*d(S,beta)
	  END DO
	  a(S,alpha) = a(S,alpha) - sum
	END DO


END SUBROUTINE new_a
