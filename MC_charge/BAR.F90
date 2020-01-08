SUBROUTINE BAR_sample(i,en, hist4, hist5, hist6)

USE system
USE conf

IMPLICIT NONE

DOUBLE PRECISION, INTENT (IN) :: en
integer, intent (inout) :: i
DOUBLE PRECISION, INTENT(OUT) :: hist4(6000), hist5(6000),hist6(6000)
DOUBLE PRECISION :: VR, VK
DOUBLE PRECISION :: para_l


i =  i + 1
hist4(i) = en   !energy in system 1

Call  RWALD ( VR, Xh, Yh, Zh)   
Call  KWALD ( VK, Xh, Yh, Zh)

hist5(i) = VR + VK  ! energy in system 0

para_l = 0.5

Call  RWALD ( VR, Xi, Yi, Zi)   
Call  KWALD ( VK, Xi, Yi, Zi)

hist6(i) = (VR+VK)*para_l



END SUBROUTINE BAR_sample


SUBROUTINE BAR_C(K,hist4,hist5,hist6,delta_F)

USE system
USE conf
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: hist4(6000), hist5(6000), hist6(6000)
DOUBLE PRECISION, INTENT(OUT) :: delta_F
integer, INTENT(IN) :: K
DOUBLE PRECISION:: delta_F1, delta_F2
integer :: i
DOUBLE PRECISION :: C1, C2, sum1, sum2, error

!initial guess for C

C1 = -80.0
error = 1

DO i=1,K
write(22,*) hist5(i)
END DO

write(22,*) ' Hist5 --- Hist4'

DO i=1,K
write(22,*) hist4(i)
END DO

write(22,*) ' Hist5 --- Hist4'
DO i=1,K
write(22,*) hist6(i)
END DO


do while (error>1.d-4)

sum1 = 0; sum2 =0;
DO i = 1, K
sum1 = sum1 + 1.d0/(1.d0 + dexp(beta*(-1.d0*hist6(i) + C1)))
sum2 = sum2 + 1.d0/(1.d0 + dexp(beta*(0.5*hist5(i) - C1)))
ENDDO

delta_F1 = dlog(sum1/sum2)/beta + C1

error = abs(delta_F1 - C1)

C1 = delta_F1

end do

error= 1.0
C2 = -50.0
do while (error>1.d-4)

sum1 = 0; sum2 =0;
DO i = 1, K
sum1 = sum1 + 1.d0/(1.d0 + dexp(beta*(-0.50*hist4(i) + C2)))
sum2 = sum2 + 1.d0/(1.d0 + dexp(beta*(hist6(i) - C2)))
ENDDO

delta_F2 = dlog(sum1/sum2)/beta + C2

error = abs(delta_F2 - C2)

C2 = delta_F2

end do

delta_F = delta_F1 + delta_F2
write(22,*) 'Hel energy=', delta_F*beta/Npart
write(*,*) sum1, sum2

end subroutine BAR_C


