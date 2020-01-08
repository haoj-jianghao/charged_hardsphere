SUBROUTINE FFT(n,coeff,dx,x,k,c,ftc)

	USE constant
	
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n
	DOUBLE PRECISION, INTENT(IN) :: coeff, dx, x(0:n), c(1:3,0:n), k(0:n)
	DOUBLE PRECISION, INTENT(OUT) :: ftc(1:3,0:n)

	DOUBLE PRECISION :: sum1, sum2, sum3
	INTEGER :: i, j, fn

! Calculate ftc(0)

	fn = (n+1)/2-1
 
 !   fn = n/2 - 1
 
!	sum = 0.0
!	DO i=1,fn
!	  sum = sum + x(i)*x(i)*c(i)
!	END DO
!	ftc(0) = coeff*dx*sum

! Calculate the remaining ftc(i)



	DO j=1,n
	  sum1 = 0.0; sum2 =0; sum3 =0;
	  DO i=1,fn
        sum1 = sum1 + x(i)*SIN(k(j)*x(i))*c(1,i)
 	    sum2 = sum2 + x(i)*SIN(k(j)*x(i))*c(2,i)
	    sum3 = sum3 + x(i)*SIN(k(j)*x(i))*c(3,i)
	  END DO
	  ftc(1,j) = coeff*dx*sum1/k(j)
	  ftc(2,j) = coeff*dx*sum2/k(j)
	  ftc(3,j) = coeff*dx*sum3/k(j)
	END DO	   

!	ftc(1,0) = ftc(1,1)-(ftc(1,2)-ftc(1,1))
!	ftc(2,0) = ftc(2,1)-(ftc(2,2)-ftc(2,1))
!	ftc(3,0) = ftc(3,1)-(ftc(3,2)-ftc(3,1))

END SUBROUTINE FFT


