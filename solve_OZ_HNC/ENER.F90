SUBROUTINE ENERGY(g, En, Pressure)

USE main_variables
USE state_variables
USE constant

IMPLICIT NONE

DOUBLE PRECISION, INTENT (INOUT) :: g(1:3,0:n)
DOUBLE PRECISION, INTENT (OUT) :: En, Pressure
DOUBLE PRECISION :: Debye_l
DOUBLE PRECISION :: gab(0:n), gab1(0:n), sum, sum1, gabh(1:3,0:n), M(3), h(1:3,0:n), osm
INTEGER :: i,j, jj

jj = anint(1.0/dr)
DO j = 1, 3
g(j,jj) = 2.0*g(j,jj+1) -  g(j,jj+2) 
END DO

Debye_l = g(1,180)*g(3,180)

DO i = jj, n
gab(i) = ( g(1,i) )*r(i)
ENDDO

DO i = jj, n
gab1(i) =  g(3,i) * r(i)
END DO


sum = 0
DO i = jj, n-1
sum = sum + (gab(i) + gab(i+1))*dr/2.0
END DO

sum1 = 0
DO i = jj, n-1
sum1 = sum1 + (gab1(i) + gab1(i+1))*dr/2.0
END DO


En = (sum-sum1)*rho/tstar*pi
Pressure = 1.0 + En/3.0 + 2*rho*pi/6.0*(g(1,jj) + g(3,jj))
write(*,*) pressure

! Check the 0th moment
DO j = 1, 3
DO i = 0, n
h(j,i) = g(j,i) - 1.d0
END DO
END DO


 DO j = 1, 3
   DO i = 0, n
    gabh(j,i) = r(i)**2.0*h(j,i)
   END DO
 END DO
    
  M = 0
  DO j = 1, 3
  DO i = 0, n-1
    M(j) = M(j) + gabh(j, i) + gabh(j,i+1)
  END DO
  END DO
  
  M = M*dr/2.0*4.0*pi
  
  sum = 1 + rho1*(M(1) - M(3))
  write(*,*) sum
  sum = 1 + rho2*(M(2) - M(3))
  write(*,*) sum
  
    
    

  

!write(*,*) 'Internal energy = ', En
!write(*,*) 'Pressure = ',1.0 + En/3.0 + 2*rho*pi/6.0*(g(1,jj) + g(3,jj))

END SUBROUTINE ENERGY







