SUBROUTINE decompose_gamma(numnod,limit,n,Q,P,gamma,a,delta_gamma)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: numnod, limit, n
	DOUBLE PRECISION, INTENT(IN) :: Q(1:3,1:numnod,0:n), P(1:3,1:numnod,0:n), gamma(1:3,0:n)
	DOUBLE PRECISION, INTENT(OUT) :: a(1:3,1:numnod), delta_gamma(1:3,0:n)

	DOUBLE PRECISION :: sum1, sum2, sum3
	INTEGER :: alpha, i

! Calculate a(alpha)
	
	a = 0.0
	DO alpha=1,numnod
	  DO i=0,limit
	    a(1,alpha) = a(1,alpha) + Q(1,alpha,i)*gamma(1,i)
	    a(2,alpha) = a(2,alpha) + Q(2,alpha,i)*gamma(2,i)
	    a(3,alpha) = a(3,alpha) + Q(3,alpha,i)*gamma(3,i)
	  END DO
	END DO

! Calculate delta_gamma(i)

	DO i=0,n
	  sum1 = 0.0; sum2 = 0; sum3 = 0
	  DO alpha=1,numnod
	    sum1 = sum1 + a(1,alpha)*P(1,alpha,i)
	    sum2 = sum2 + a(2,alpha)*P(2,alpha,i)
	    sum3 = sum3 + a(3,alpha)*P(3,alpha,i)
	  END DO
	  delta_gamma(1,i) = gamma(1,i) - sum1
	  delta_gamma(2,i) = gamma(2,i) - sum2
	  delta_gamma(3,i) = gamma(3,i) - sum3
	END DO

END SUBROUTINE decompose_gamma