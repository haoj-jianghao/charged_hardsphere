SUBROUTINE NR_test(numnod,a,a_prime,d,converge)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: numnod
	DOUBLE PRECISION, INTENT(IN) :: a(1:3,1:numnod), a_prime(1:3,1:numnod)
	DOUBLE PRECISION, INTENT(OUT) :: d(1:3,1:numnod)
	LOGICAL, INTENT(OUT) :: converge

	INTEGER :: alpha,i
	DOUBLE PRECISION :: sum, dtest

! Test for Newton-Raphson cycle convergence

	dtest = 1.0
	converge = .false.
!	d = a - a_prime
	sum	= 0.0
	
	DO i = 1, 3
	DO alpha=1,numnod
	  d(i,alpha) = a(i,alpha) - a_prime(i,alpha)
	  sum = sum + d(i,alpha)*d(i,alpha)
	END DO
	END DO

	IF(sum < dtest) GOTO 522

!	PRINT '(1X,"d(alpha) limit surpassed")'
    WRITE (*,*) 'limit surpassed' , sum
    
 !   DO i = 1, 3
!	DO alpha=1,numnod
!	  d(i,alpha) = SQRT(dtest)*d(i,alpha)/SQRT(sum)
!	END DO
!	END DO

522	IF(SQRT(sum) < 1.0e-5) converge = .true.

END SUBROUTINE NR_test


