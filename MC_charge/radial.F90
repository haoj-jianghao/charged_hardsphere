
      SUBROUTINE radial(Iout,hist1,hist2,hist3, hist11, hist12)
      use parameters
      use conf
      use potential
      use system
     
      IMPLICIT NONE

      INTEGER Iout, i,j,k,l,n,m,p
      INTEGER maxbin
      PARAMETER (maxbin=6000)
      REAL*8 gr1(maxbin),gr11(maxbin),gr12(maxbin),gr22(maxbin)
      INTEGER bin
      INTEGER hist1(maxbin),hist2(maxbin),hist3(maxbin), hist11(maxbin), hist12(maxbin)
      REAL*8 xj,yj,zj,xk,yk,zk,boxf,rhof,dr
      REAL*8 dx,dy,dz
      REAL*8 rjk,width,r
      REAL*8 rlower,rupper
      REAL*8 factor
      REAL*8 volume,expect,tmpSum2
      integer a1(4),a2(4),a11(4)
      INTEGER near; near=0;
 
    
      width=0.025
          
          DO i = 1, 4
          a1(i) = hist1(int(1.0/width)+i);
          END DO
          
          DO i = 1,4
          a2(i) = hist2(int(1.0/width)+i)
          END DO
          
          DO i = 1,4
          a11(i) = hist11(int(1.0/width)+i)
          END DO
          
	 
	     	
	 DO j = 1, NPART-1
               
                  xj = X(j)
                  yj = Y(j)
                  zj = Z(j)


!C	  periodic boundary conditions

	       IF (xj.GT.HBOX) THEN
              xj = xj - BOX	        
            ELSE
               IF (xj.LT.-HBOX) xj = xj + BOX
            END IF
            IF (yj.GT.HBOX) THEN
               yj = yj - BOX
            ELSE
               IF (yj.LT.-HBOX) yj = yj + BOX
            END IF
            IF (zj.GT.HBOX) THEN
               zj = zj - BOX
            ELSE
               IF (zj.LT.-HBOX) zj = zj + BOX
            END IF
                                                
                        DO k = j+1, NPART
                               
							 
	                     xk= x(k)
						 yk= y(k) 
                         zk= z(k)
							 
!C	  periodic boundary conditions

							   IF (xk.GT.HBOX) THEN
                                    xk = xk - BOX
	                           ELSE
                                    IF (xk.LT.-HBOX) xk = xk + BOX
                               END IF

                               IF (yk.GT.HBOX) THEN
                                     yk = yk - BOX
                               ELSE
                                    IF (yk.LT.-HBOX) yk = yk + BOX
                               END IF

                               IF (zk.GT.HBOX) THEN
                                      zk = zk - BOX
                               ELSE
                                   IF (zk.LT.-HBOX) zk = zk + BOX
!
                               END IF

                           dx = xk - xj
                           dy = yk - yj
                           dz = zk - zj

!C                   Minimum image convention

                             IF (ABS(dx).GT.HBOX) THEN
                                 dx = dx - BOX*ANINT(dx/BOX)
                                             
                               END IF
                              
                             IF (ABS(dy).GT.HBOX) THEN
                                  dy = dy - BOX*ANINT(dy/BOX)     
                                END IF
                              IF (ABS(dz).GT.HBOX) THEN
                                   dz = dz - BOX*ANINT(dz/BOX)      
                                END IF
  
	                   
                    rjk = SQRT (dx*dx + dy*dy + dz*dz)         
                                            
					bin = INT((rjk)/width) + 1        
				   	 
	               IF (bin .LE. maxbin .and. j.LE.Nelec .and. k.LE.Nelec) THEN
	                     IF (MOD(j,2)==0 .AND. MOD(k,2)==0) THEN            
		                 hist1(bin) = hist1(bin) + 2                                              
	                     END IF
	                     IF (MOD(j,2)==0 .AND. MOD(k,2)==1) THEN	                       
		                 hist2(bin)= hist2(bin) + 2                                               
	                     END IF	                    	                     
	                     IF (MOD(j,2)==1 .AND. MOD(k,2)==1) THEN
		                 hist3(bin) = hist3(bin) + 2                                               
	                     END IF
                   END IF
                  
                  IF (Nelec < Npart) THEN
                  IF (bin .LE. maxbin .and. j>Nelec .and. k>Nelec) THEN
                        hist11(bin) = hist11(bin) + 2         !neutral - neutral
                  END IF                  
                  IF (bin .LE. maxbin .and. j.LE.Nelec .and. k>Nelec .and. MOD(j,2)==0) THEN
                        hist12(bin) = hist12(bin) + 1
                  END IF
                  ELSE
                        hist11 = 0; hist12 = 0;
                  END IF
                   
	                           
              END DO
  
          END DO                                
            
                  
        write(23,*) hist1(int(1.0/width) + 1 ) - a1(1), hist1(int(1.0/width) + 2 ) - a1(2), hist1(int(1.0/width) + 3 ) - a1(3), hist1(int(1.0/width) + 4 ) - a1(4)
        write(24,*) hist2(int(1.0/width) + 1 ) - a2(1), hist2(int(1.0/width) + 2 ) - a2(2), hist2(int(1.0/width) + 3 ) - a2(3), hist2(int(1.0/width) + 4 ) - a2(4)   
        write(26,*) hist11(int(1.0/width) + 1 ) - a11(1), hist11(int(1.0/width) + 2 ) - a11(2), hist11(int(1.0/width) + 3 ) - a11(3), hist11(int(1.0/width) + 4 ) - a11(4)
      
      END SUBROUTINE radial