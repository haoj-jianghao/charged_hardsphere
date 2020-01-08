      SUBROUTINE SAMPLE(I, En, Vir)
!c
!c      write quantities (pressure and energy) to file
!c
!c
!c  Ener (input) : total energy
!c  Vir  (input) : total virial
!c
!c
      USE parameters
      USE conf
      USE system
      USE potential
      
      IMPLICIT NONE

      INTEGER I,J
      DOUBLE PRECISION enp, press, CORP, vol, rho
      DOUBLE PRECISION En, Vir
      
      IF (NPART.NE.0) THEN
         enp = En/DBLE(NPART)
         vol = BOX**3
         press = (NPART/vol)/BETA + Vir/(3*vol)
         rho = NPART/vol
         IF (TAILCO) press = press          
         
	   ELSE
         enp = 0
         press = 0
      END IF
      WRITE (66, *) I, enp, press
      RETURN
      END



!**==store.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
      SUBROUTINE STORE(Iout, Dr)
!c
!c     writes configuration to disk
!c
!c  Iout (input) file number
!c  Dr   (input) maximum displacement
!c
!c
      USE parameters
      USE conf
      USE system
      
      IMPLICIT NONE

      INTEGER Iout, i
      DOUBLE PRECISION Dr
 
      WRITE (Iout, *) BOX, HBOX
      WRITE (Iout, *) NPART
      WRITE (Iout, *) Dr
      DO i = 1, NPART
         WRITE (Iout, *) X(i), Y(i), Z(i)
      END DO
      REWIND (Iout)
      RETURN
      END


!**==rantest.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
 
      SUBROUTINE RANTEST(Iseed)
!c
!c     test and initialize the random number generator
!c
      IMPLICIT NONE
      INTEGER Iseed, i
      DOUBLE PRECISION RANF
 
      CALL RANSET(Iseed)
      PRINT *, ' ******** test random numbers ***********'
      DO i = 1, 5
         PRINT *, ' i,ranf() ', i, RANF(Iseed)
      END DO
      RETURN
      END
      
      SUBROUTINE LATTICE
!c
!c     place `npart' particles on a simple cubic
!c     lattice with density 'rho'

      USE parameters
      USE conf
      USE system

      IMPLICIT NONE

      INTEGER i, j, k, itel, n
      DOUBLE PRECISION dx, dy, dz, del
 
      n = INT(NPART**(1./3.)) !+ 1
      IF (n.EQ.0) n = 1
      del = BOX/DBLE(n)
      itel = 0
      dx = -del
      DO i = 1, n
         dx = dx + del
         dy = -del
         DO j = 1, n
            dy = dy + del
            dz = -del
            DO k = 1, n
               dz = dz + del
               IF (itel.LT.NPART) THEN
                  itel = itel + 1
                  X(itel) = dx
                  Y(itel) = dy
                  Z(itel) = dz
                  Xh(itel) = dx
                  Yh(itel) = dy
                  Zh(itel) = dz
                  Xi(itel) = dx
                  Yi(itel) = dy
                  Zi(itel) = dz
			   If (MOD(itel,2)==0) then
		           Q(itel)=1
		       Else 
		           Q(itel)=-1
		       End if
       
			   write (3,10) x(itel),y(itel),z(itel),Q(itel)
10		   format(2X,3(2X,F10.6))


               END IF
            END DO
         END DO
      END DO
      
      CALL FCC(BOX,8, NPART)
      ratio_ek = 0.5           !Nelec/Npart without neutral --- ratio ==1
      Nelec = int(Npart * ratio_ek)

        DO  i = 1,NPART
         	   If (MOD (i,2) ==0 .and. i<=Nelec) then
		           Q(i)=1
    	       ELSEif (MOD (i,2) ==1 .and. i<=Nelec) then
		           Q(i)=-1
		       ELSEIf (i>Nelec) then
		           Q(i) = 0;
		       End if
		END DO



      WRITE (6, 99001) itel
      RETURN
99001 FORMAT (' Initialisation on lattice: ', /, i10, & 
     &        ' particles placed on a lattice')
      END
      
      
       SUBROUTINE FCC(box_side,ncells_per_side,npart)

        USE parameters
        USE conf
        USE system
        
        INTEGER :: ncells_per_side
        DOUBLE PRECISION :: box_side
  !     DOUBLE PRECISION, INTENT(OUT) :: x(npart), y(npart), z(npart)

        DOUBLE PRECISION :: cell_side
        INTEGER :: ix, iy, iz, iref, m

! Calculate the side of the unit cell

        cell_side  = box_side/DBLE(ncells_per_side)

! Build the unit cell

! Point at the origin

        x(1) =  0.0;                y(1) =  0.0;                z(1) =  0.0

! Point at the face center on the XY plane

        x(2) =  0.5*cell_side;      y(2) =  0.5*cell_side;      z(2) =  0.0

! Point at the face center on the YZ plane

        x(3) =  0.0;                y(3) =  0.5*cell_side;      z(3) =  0.5*cell_side

! Point at the face center on the XZ plane

        x(4) =  0.5*cell_side;      y(4) =  0.0;                z(4) =  0.5*cell_side

! Construct the lattice from the four points above

        m = 0
        DO iz = 1, ncells_per_side
           DO iy = 1, ncells_per_side
              DO ix = 1, ncells_per_side
                 DO iref = 1, 4
                    x(iref+m) = x(iref) + cell_side*DBLE(ix-1)
                    y(iref+m) = y(iref) + cell_side*DBLE(iy-1)
                    z(iref+m) = z(iref) + cell_side*DBLE(iz-1)
                 END DO
                 m = m + 4
              END DO
           END DO
        END DO
        
        xi = x; yi = y; zi = z;
        xh = x; yh = y; zh =z;

        RETURN

    END




!DO j = 1, 6
!DO i=1,Npart
! q_scale(j,i) = Q(i)
!ENDDO
!ENDDO

!DO i = 1,6
!scale_lambda(i) = (i-1)*1.d0/5.d0
!END DO

!scale the counter ions
!DO j = 1, 6
!DO i = 1, Npart
!IF (MOD(i,2)==1) THEN
!q_scale(j,i) = q_scale(j,i)!*(1.d0 + 2.d0*scale_lambda(j)/Npart)
!END IF
!END DO
!q_scale(j,Npart+1) = 1.0  !insert a positive charge
!END DO
