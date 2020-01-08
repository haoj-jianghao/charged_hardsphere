!**==toterg.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
      SUBROUTINE TOTERG(Ener, Vir)
!c
!c     calculates total energy
!c
!c  Ener (output) : total energy
!c  Vir  (output) : total virial
!c
      USE parameters
      USE conf
      USE potential
      USE system
 
      IMPLICIT NONE

 
      
      INTEGER i, jb
      DOUBLE PRECISION Ener, eni, viri, Vir
 
      Ener = 0
      Vir = 0

      CALL ENERI(X, Y, Z, eni, viri)

      Ener = eni

      
      RETURN
      END