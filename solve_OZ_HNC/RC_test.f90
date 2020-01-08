SUBROUTINE RC_test(n,rc,dr,gamma,gamma_prime,eta,converge)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, rc
	DOUBLE PRECISION, INTENT(IN) :: dr, gamma(1:3,0:n), gamma_prime(1:3,0:n)
	DOUBLE PRECISION, INTENT(OUT) :: eta(rc)
	LOGICAL, INTENT(OUT) :: converge

	INTEGER :: i, j
	DOUBLE PRECISION :: sum, diff
	
! Test for convergence of a refinement cycle

	converge = .false.
	sum = 0.0
	DO j=1,3
	DO i=0,n
	  diff = gamma_prime(j,i)-gamma(j,i)
	  sum = sum + diff*diff
	END DO
	END DO
	
	eta(rc) = SQRT(dr*sum)

	IF(eta(rc) < 1d-5) converge = .true.

END SUBROUTINE RC_test 