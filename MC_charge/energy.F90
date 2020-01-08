      SUBROUTINE ENER(En,r2,i,j)
!c     calculate energy
!c
!c     En : (output) energy
!c     Vir: (output) virial
!c     R2 : (input) distance squared between two particles

       USE conf
       USE system
       USE parameters
       USE potential

      IMPLICIT NONE
      DOUBLE PRECISION R2, r2i, r6i, Vir,Enhs,xi,yi,zi
	  DOUBLE PRECISION KAPPA
	  REAL VR, VK, En


	
	  INTEGER i,j,n
 
         IF (R2 .LT. sig2) THEN
         En = 10**6
         Vir = 0
	
	   ELSE 
          En = 0
          Vir = 0
	
	   END IF

      RETURN
      END