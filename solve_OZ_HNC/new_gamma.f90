SUBROUTINE new_gamma(n,numnod,P,a,delta_gamma,gamma)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, numnod
	DOUBLE PRECISION, INTENT(IN) :: P(1:3,1:numnod,0:n), a(1:3,1:numnod), delta_gamma(1:3,0:n)
	DOUBLE PRECISION, INTENT(OUT) :: gamma(1:3,0:n)

	INTEGER :: i, alpha
	DOUBLE PRECISION :: sum1, sum2, sum3

	DO i=0,n
	  sum1 = 0.0;sum2 =0.0; sum3 =0.0
	  DO alpha=1,numnod
	    sum1 = sum1 + a(1,alpha)*P(1,alpha,i)
	    sum2 = sum2 + a(2,alpha)*P(2,alpha,i)
	    sum3 = sum3 + a(3,alpha)*P(3,alpha,i)
	  END DO
	  gamma(1,i) = sum1 + delta_gamma(1,i)
	  gamma(2,i) = sum2 + delta_gamma(2,i)
	  gamma(3,i) = sum3 + delta_gamma(3,i)
	END DO

END SUBROUTINE new_gamma
