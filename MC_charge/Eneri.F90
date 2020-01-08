   
      SUBROUTINE ENERI(XX, YY, ZZ, En, Vir)
      
      
      USE conf
      USE system
      USE parameters
      USE potential  
      USE K_VEC  

      IMPLICIT NONE 

      
!      USE conf
    
      DOUBLE PRECISION,INTENT(IN) ::  XX(NPART), YY(NPART), ZZ(NPART)
      
      DOUBLE PRECISION, INTENT(OUT) :: En,Vir

      DOUBLE PRECISION :: dx, dy, dz, r2
      
      DOUBLE PRECISION :: VR, Vk,En1,enii, virij, enij

      INTEGER I, j
  
      En  = 0
      Vir = 0
	  En1 = 0
	  VR = 0; VK = 0;
      DO I = 1, NPART-1
       DO j = I + 1 , NPART
            dx = XX(i) - XX(j)
            dy = YY(i) - YY(j)
            dz = ZZ(i) - ZZ(j)
            IF (dx.GT.HBOX) THEN
               dx = dx - BOX
            ELSE
               IF (dx.LT.-HBOX) dx = dx + BOX
            END IF
           IF (dy.GT.HBOX) THEN
               dy = dy - BOX
            ELSE
               IF (dy.LT.-HBOX) dy = dy + BOX
            END IF
            IF (dz.GT.HBOX) THEN
              dz = dz - BOX 
            ELSE
               IF (dz.LT.-HBOX) dz = dz + BOX
            END IF
            r2 = dx*dx + dy*dy + dz*dz
  
            if (r2 .LE. 1.0) then
               en = en + 10**6.0
               r2 = 1.d0
            else
                en = 0 + en              
            end if
            
            en = en + Q(I)*Q(J)/dsqrt(r2)
	 
      END DO
      END DO
      
      
      KAPPA= 7.d0/BOX     
	  
!	  Call  SETUP ( KAPPA )
     
 !     Call  RWALD ( VR, XX, YY, ZZ)
     
  !    Call  KWALD ( VK, XX, YY, ZZ)
 
   !   call simple_ener(xx, yy, zz, en)
      
   !   En = En + VR + VK
  

!C	  En=(En+En1)
!C      write (10,*) En, En1 

 !    RETURN
      END SUBROUTINE ENERI