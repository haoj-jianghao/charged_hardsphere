 
      SUBROUTINE ADJUST(Attemp, Nacc1, Nacc2, Nacc3, Dr1, Dr2, Dr3)
!c
!c     adjusts maximum displacement such that 50% of the
!c     movels will be accepted
!c
!c  Attemp (input) number of attemps that have been performed to displace a particle
!c  Nacc   (input) number of successful attemps to displace a particle
!c  Dr     (output) new maximum displacement
!c
      USE system
      USE conf
      
      
      IMPLICIT NONE
     
      INTEGER, INTENT(INOUT):: Attemp, Nacc1, Nacc2, Nacc3
      DOUBLE PRECISION, INTENT(INOUT) :: dr1, dr2, dr3
      INTEGER :: attempp, naccp1, naccp2, naccp3
    
      DOUBLE PRECISION :: dro1,dro2,dro3, frac1, frac2, frac3
     
      SAVE naccp1, naccp2, naccp3, attempp
 
      IF (Attemp.EQ.0.OR.attempp.GE.Attemp) THEN
         naccp1 = Nacc1
         naccp2 = Nacc2
         naccp3 = Nacc3
         attempp = Attemp
      ELSE
         frac1 = DBLE(Nacc1-naccp1)/DBLE(Attemp-attempp)
         frac2 = DBLE(Nacc2-naccp2)/DBLE(Attemp-attempp)
         frac3 = DBLE(Nacc3-naccp3)/DBLE(Attemp-attempp)
         dro1 = Dr1
         dro2 = Dr2
         dro3 = Dr3
         Dr1 = Dr1*ABS(frac1/0.5D0)
         Dr2 = Dr2*ABS(frac2/0.5D0)
         Dr3 = Dr3*ABS(frac3/0.5D0)
         
!c        ---limit the change:
         IF (Dr1/dro1.GT.1.5D0) Dr1 = dro1*1.5D0
         IF (Dr1/dro1.LT.0.5D0) Dr1 = dro1*0.5D0
         IF (Dr1.GT.HBOX/2.D0) Dr1 = HBOX/2.D0
         
         IF (Dr2/dro2.GT.1.5D0) Dr2 = dro2*1.5D0
         IF (Dr2/dro2.LT.0.5D0) Dr2 = dro2*0.5D0
         IF (Dr2.GT.HBOX/2.D0) Dr2 = HBOX/2.D0
         
         IF (Dr3/dro3.GT.1.5D0) Dr3 = dro3*1.5D0
         IF (Dr3/dro3.LT.0.5D0) Dr3 = dro3*0.5D0
         IF (Dr3.GT.HBOX/2.D0) Dr3 = HBOX/2.D0
         
         WRITE (6, 99001) Dr1, dro1, frac1, Attemp - attempp, Nacc1 - naccp1
         WRITE (6, 99001) Dr2, dro2, frac2, Attemp - attempp, Nacc2 - naccp2
         WRITE (6, 99001) Dr3, dro3, frac3, Attemp - attempp, Nacc3 - naccp3
!c        ---store nacc and attemp for next use
         naccp1 = Nacc1
         naccp2 = Nacc2
         naccp3 = Nacc3
         attempp = Attemp
      END IF
      RETURN
99001 FORMAT (' Max. displ. set to : ', f6.3, ' (old : ', f6.3, ')', /, &
     &        ' Frac. acc.: ', f5.4, ' attempts: ', i7, ' succes: ', i7)
     
      END SUBROUTINE ADJUST
 