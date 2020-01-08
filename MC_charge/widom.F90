SUBROUTINE WIDOM (En,wtest1,wtest2,wtest3,inser)

USE system
USE conf
USE K_VEC
USE potential

IMPLICIT NONE

DOUBLE PRECISION, INTENT(INOUT)::wtest1(6),wtest2(6),wtest3
REAL, INTENT(IN) :: EN
INTEGER, INTENT (INOUT) :: inser

DOUBLE PRECISION :: scale_lambda(6)
DOUBLE PRECISION :: xn,yn,zn
DOUBLE PRECISION :: En_test, En_elec(6)
integer :: particle
integer :: Iseed
integer :: KMAX, KSQMAX
DOUBLE PRECISION ::  RANF
DOUBLE PRECISION :: Q_charge, TWOPI, RSQPI,r2, QI,QJ
PARAMETER ( KMAX = 7, KSQMAX = 125 )
DOUBLE PRECISION :: VR(6), VK(6), VD,VS
REAL ::      RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, RKX,RKY,RKZ
REAL ::      RIJSQ, RIJ, KRIJ, ERFC, VIJ, FACTOR
INTEGER ::   I,J, K, KX, KY, KZ, TOTK, KSQ, N

COMPLEX     EIKX(1:Npart+1, 0:KMAX)
COMPLEX     EIKY(1:Npart+1, -KMAX:KMAX)
COMPLEX     EIKZ(1:Npart+1, -KMAX:KMAX)
COMPLEX     EIKR(Npart+1), SUM

! parameters for ewald summation

TWOPI = 6.2831853
RSQPI = 0.5641896
Iseed = 32657

N = Npart
DO i = 1,6
scale_lambda(i) = 1.0  !(i-1)*1.d0/5.d0
END DO


! Only 1 ghost particle was inserted into the box.
Q_charge = 1.0
Inser = Inser + 1

xn = RANF(Iseed) * BOX
yn = RANF(Iseed) * BOX
zn = RANF(Iseed) * BOX

VR = 0.0
VK = 0.0




En_test = 0.0

!! calculate the interaction between ghost particle with real particles    
!calculate the interaction in Rspace

DO i = 1 , Npart        
        
              RXIJ = xn - X(I)
              RYIJ = yn - Y(I)
              RZIJ = zn - Z(I)
              
	          IF (ABS(RXIJ).GT.HBOX) THEN
                 RXIJ = RXIJ - BOX*ANINT(RXIJ/BOX)                             
              END IF                              
              IF (ABS(RYIJ).GT.HBOX) THEN
                 RYIJ = RYIJ - BOX*ANINT(RYIJ/BOX)     
              END IF            
  		      IF (ABS(RZIJ).GT.HBOX) THEN
                 RZIJ = RZIJ - BOX*ANINT(RZIJ/BOX)      
              END IF 
              r2 = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
              
              IF (R2 .LT. sig2) THEN
              En_test = 10**6 + En_test
	          ELSE 
              En_test = 0 + En_test
              ENDIF
 ENDDO

IF (En_test<1.d0) THEN
!===============================================================================

DO K = 1, 6 
     VR(k) = 0;
     DO I = 1, N
           RXI = X(I)
           RYI = Y(I)
           RZI = Z(I)
           QI  = Q_scale(k,i)         
           DO J = I + 1, N + 1
              IF (j>N) THEN
              RXIJ = RXI - Xn
              RYIJ = RYI - Yn
              RZIJ = RZI - Zn
              ELSE
              RXIJ = RXI - X(J)
              RYIJ = RYI - Y(J)
              RZIJ = RZI - Z(J)
              END IF
              RXIJ = RXIJ - BOX*ANINT ( RXIJ /BOX )
              RYIJ = RYIJ - BOX*ANINT ( RYIJ /BOX )
              RZIJ = RZIJ - BOX*ANINT ( RZIJ /BOX )
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
              RIJ   = SQRT ( RIJSQ )
              KRIJ  = KAPPA * RIJ              
              VIJ   = QI * Q_scale(k,j) * ERFC ( KRIJ ) / RIJ
              VR(k)    = VR(k) + VIJ
          ENDDO
     ENDDO
  !===============================================================================
  !calculate the interaction in Kspace       
       DO  I = 1, N + 1
           EIKX(I, 0) = (1.0, 0.0)
           EIKY(I, 0) = (1.0, 0.0)
           EIKZ(I, 0) = (1.0, 0.0)
           IF(I>N) THEN
           EIKX(I, 1) = CMPLX ( COS ( TWOPI * Xn /BOX ) ,  SIN ( TWOPI * Xn /BOX ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPI * Yn /BOX ) ,  SIN ( TWOPI * Yn /BOX ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPI * Zn /BOX ) ,  SIN ( TWOPI * Zn /BOX ) )
           ELSE
           EIKX(I, 1) = CMPLX ( COS ( TWOPI * X(I) /BOX ) ,  SIN ( TWOPI * X(I) /BOX ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPI * Y(I) /BOX ) ,  SIN ( TWOPI * Y(I) /BOX ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPI * Z(I) /BOX ) ,  SIN ( TWOPI * Z(I) /BOX ) )
           ENDIF
           EIKY(I, -1) = CONJG ( EIKY(I, 1) )
           EIKZ(I, -1) = CONJG ( EIKZ(I, 1) )
      ENDDO

     DO KX = 2, KMAX
           DO I = 1, N + 1 
              EIKX(I, KX) = EIKX(I, KX-1) * EIKX(I, 1)
           ENDDO
     ENDDO

     DO KY = 2, KMAX
           DO I = 1, N + 1
              EIKY(I,  KY) = EIKY(I, KY-1) * EIKY(I, 1)
              EIKY(I, -KY) = CONJG ( EIKY(I, KY) )
           ENDDO
     ENDDO

      DO KZ = 2, KMAX
           DO I = 1, N + 1 
              EIKZ(I,  KZ) = EIKZ(I, KZ-1) * EIKZ(I, 1)
              EIKZ(I, -KZ) = CONJG ( EIKZ(I, KZ) )
           ENDDO
      ENDDO

!C    ** SUM OVER ALL VECTORS **
        VD   = 0.0
        TOTK = 0

        DO KX = 0, KMAX
           IF ( KX .EQ. 0 ) THEN
              FACTOR = 1.0
           ELSE
              FACTOR = 2.0
           ENDIF
           DO  KY = -KMAX, KMAX
              DO  KZ = -KMAX, KMAX
                 KSQ = KX * KX + KY * KY + KZ * KZ
                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN
                    TOTK = TOTK + 1
                    SUM  = (0.0, 0.0)
                    DO I = 1, N + 1
                       EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ)
                       SUM     = SUM + Q_scale(K,I) * EIKR(I)
                   ENDDO
                    VD = VD + FACTOR * KVEC(TOTK) * CONJG ( SUM ) * SUM
                 ENDIF
            ENDDO
          ENDDO
      ENDDO
      VS = 0.0
        DO I = 1, N+1
           VS = VS + Q_scale(k,I) * Q_scale(k,I)
        ENDDO
        VS = RSQPI * KAPPA * VS
        VD = VD/BOX/BOX/BOX
               
        VK(k) = VD  - VS 
        
        EN_elec(K) = VK(K) + VR(K) - EN    
 ENDDO
        write(22,*) wtest1(6)/wtest2(6)
 ENDIF     


   !   DO k = 1, 6
   !   EN_elec(k) = En_elec(k) 
   !   ENDDO
      
      wtest3 = wtest3 + dexp(-beta*En_test)

      DO i = 1,6
  !    wtest1(i) = wtest1(i) + beta*EN_elec(i)*dexp(-beta*(En_test + scale_lambda(i)*EN_elec(i)))
      wtest2(i) = wtest2(i) + dexp(-beta*(En_test+scale_lambda(i)*EN_elec(i)))
      END DO
      
      
    !  IF (abs(wtest1(6))>1000) THEN
    !  write(*,*) '123'
    !  END IF
 !    wtest = wtest + dexp(-beta*en_test)
      
 END SUBROUTINE WIDOM
      
      
      
      
      
 SUBROUTINE chemical_potential (wtest1,wtest2,wtest3,inser, chem_p)
 
 USE system
 
 IMPLICIT NONE
 
 DOUBLE PRECISION, INTENT (INOUT) :: wtest1(6),wtest2(6),wtest3
 INTEGER, INTENT (IN) :: inser
 DOUBLE PRECISION, INTENT (OUT) :: chem_p
 DOUBLE PRECISION :: scale_lambda(6),wtest(6)
 DOUBLE PRECISION :: SUM
 INTEGER :: i
 
 

 wtest1 = wtest1/inser
 wtest2 = wtest2/inser
 wtest3 = -dlog(wtest3/inser)
 
 wtest = wtest1/wtest2
 SUM = 0
 
 DO i = 1, 5
  SUM = SUM + (wtest(i) + wtest(i+1))*1.0/5.0/2.0
 END DO
 
 
 !chem_p = (sum + wtest3)/beta
 chem_p =  wtest2(6)
 write(22,*) chem_p
 END SUBROUTINE chemical_potential
   
   
   
 
 
 
 

     
