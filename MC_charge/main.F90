      MODULE system
      DOUBLE PRECISION :: box,temp,beta,hbox
      END MODULE SYSTEM
      
      MODULE conf
      DOUBLE PRECISION, ALLOCATABLE :: x(:),y(:),z(:),q(:),q_scale(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: xh(:),yh(:),zh(:)
      DOUBLE PRECISION, ALLOCATABLE :: xi(:),yi(:),zi(:)
      integer:: npart, Nelec
      double precision :: ratio_ek
      END MODULE conf
      
      MODULE parameters
      integer:: npmax
      parameter (npmax=10000)
      END MODULE parameters
      
      MODULE potential
      double precision:: eps4,sig2,mass,ecut,rc2,pi,eps48,rc
      logical :: tailco,shift
      END MODULE potential
      
      MODULE K_vec
      DOUBLE PRECISION :: KVEC(2000), KAPPA
      END MODULE K_vec
      
      
      
      PROGRAM MC_NVT
      
       USE conf
       USE system
       USE parameters
       USE potential
       USE K_vec
 
       IMPLICIT NONE
      
      INTEGER :: iseed, equil, prod, nsamp, ii, icycl, ndispl, attempt, nacc1, nacc2, nacc3, ncycl, nmoves, imove, index
      DOUBLE PRECISION :: en, ent, vir, virt, en_i
      DOUBLE PRECISION :: delta_F
      DOUBLE PRECISION :: dr1,dr2,dr3,wtest1(6),wtest2(6),wtest, chem_p, wtest3
	  INTEGER :: maxbin,n
      PARAMETER (maxbin=6000)
      REAL*8 ::  gr1(maxbin),gr2(maxbin),gr3(maxbin),gr11(maxbin),gr12(maxbin)
      INTEGER ::  bin
      INTEGER :: hist1(maxbin),hist2(maxbin),hist3(maxbin), hist11(maxbin), hist12(maxbin)   
      INTEGER :: insertions
      DOUBLE PRECISION :: hist4(6000), hist5(6000),hist6(6000)
	  REAL*8 :: rhof
      REAL*8 :: rjk,width,r
      REAL*8 :: rlower,rupper
      REAL*8 :: factor, factor1, factor2
      REAL*8 :: volume,expect , expect1 , expect2
      integer :: clck_counts_beg, clck_counts_end, clck_rate
      call system_clock ( clck_counts_beg, clck_rate )


 
      WRITE (*, *) '**************** MC_NVT ***************'
!     ---initialize sysem
      CALL READDAT(equil, prod, nsamp, ndispl, dr1, iseed)
      nmoves = ndispl
      dr2 = dr1*1.2
      dr3 = dr1*1.2
!     ---total energy of the system
  
      KAPPA = 7/BOX; wtest1 = 0;wtest2 = 0; wtest3 = 0; insertions = 0; index = 0

      CALL SETUP(KAPPA);
      
      CALL TOTERG(en, vir)
      WRITE (*, 99001) en, vir
      en_i = en*0.5  ! intermedia lambda = 0.5 only one stage

!     ---start MC-cycle
      DO ii = 1, 2
!        --- ii=1 equilibration
!        --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (*, *) ' Start equilibration '
         ELSE
            IF (ncycl.NE.0) WRITE (*, *) ' Start production '
            ncycl = prod
         END IF
         attempt = 0
         nacc1 = 0
         nacc2 = 0
         nacc3 = 0
!        ---intialize the subroutine that adjust the maximum displacement
         CALL ADJUST(attempt, nacc1,nacc2,nacc3, dr1, dr2,dr3)
         DO icycl = 1, ncycl
            DO imove = 1, nmoves
!              ---attempt to displace a particle
               CALL MCMOVE(en, en_i, vir, attempt, nacc1,nacc2,nacc3, dr1, dr2, dr3, iseed)             
            END DO
            
     !       WRITE(*,*) 'cycle =   ', icycl
     !       WRITE(*,*) 'en =   ', en
     !       IF (abs(en)>=1000) THEN
     !       write(*,*) '1'
     !       ENDIF                 
              
            IF (ii.EQ.2) THEN
!              ---sample averages
               IF (MOD(icycl,nsamp).EQ.0) THEN
                  CALL SAMPLE(iCycl, en, vir)
              !   CALL WIDOM(en,wtest1,wtest2,wtest3,insertions)   
                  CALL BAR_sample(index,en,hist4, hist5, hist6)                
               END IF
         
               CALL  radial(icycl,hist1,hist2,hist3, hist11, hist12)
	         END IF
            IF (MOD(icycl,ncycl/10).EQ.0) THEN
               WRITE (*, *) '======>> Done ', icycl, ' out of ', ncycl
               WRITE (*,*) 'Energy =', en
!              ---write intermediate configuration to file
               CALL STORE(icycl, dr1)
!              ---adjust maximum displacements
               CALL ADJUST(attempt, nacc1,nacc2, nacc3, dr1, dr2, dr3)
   
            END IF
         END DO
         IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (*, 99003) attempt, nacc1, &
     &                               100.*FLOAT(nacc1)/FLOAT(attempt)
!           ---test total energy
            CALL TOTERG(ent, virt)
            IF (ABS(ent-en).GT.1.D-6) THEN
               WRITE (*, *) ' ######### PROBLEMS ENERGY ################ '
                         
            END IF
            IF (ABS(virt-vir).GT.1.D-6) THEN
               WRITE (*, *) ' ######### PROBLEMS VIRIAL ################ '
     &                    
            END IF
            WRITE (*, 99002) ent, en, ent - en, virt, vir, virt - vir
         END IF
      END DO
      CALL STORE(21, dr1)
      CALL BAR_C(index,hist4,hist5,hist6,delta_F)
         
      OPEN (11, FILE='11.dat',status='old')
	
!	 READ (11, *) boxf

!    READ (11, *) NPART
      !   CALL chemical_potential (wtest1,wtest2,wtest3,insertions, chem_p)
   
         pi=3.1415926
	     width=0.025
         volume = box * box * box
         factor=(4.0d0/3.0d0) * pi * DBLE(NPART) **2 / volume*(ratio_ek/2.0)**2
         factor1=(4.0d0/3.0d0) * pi * DBLE(NPART) **2 / volume*(1-ratio_ek)**2
         factor2=(4.0d0/3.0d0) * pi * DBLE(NPART) **2 / volume*(1-ratio_ek)*(ratio_ek)
         
         DO bin = 1, (HBOX/width)
           
          rlower = DBLE(bin-1) * width
	      rupper = rlower + width
	
          expect = factor * (rupper**3 - rlower**3)
          expect1 = factor1 * (rupper**3 - rlower**3)
          expect2 = factor2 * (rupper**3 - rlower**3)
          
          gr1(bin) = DBLE(hist1(bin)) /expect/ncycl
	      gr2(bin) = DBLE(hist2(bin)) /expect/ncycl
	      gr3(bin) = DBLE(hist3(bin)) /expect/ncycl
	      IF (ratio_ek .NE. 1.0) THEN
	      gr11(bin) = DBLE(hist11(bin)) /expect1/ncycl
	      gr12(bin) = DBLE(hist12(bin)) /expect2/ncycl
	      END IF
            
	      r=width*bin - width/2.0
	      
          WRITE (1, *) bin, r, gr1(bin),gr2(bin),gr3(bin),gr11(bin), gr12(bin)
   10     FORMAT (2X,(i10, 4(2X, F8.5)) )
       END DO
      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)  (clck_counts_end - clck_counts_beg) / real (clck_rate)

	STOP

99001 FORMAT (' Total energy initial configuration: ', f12.5, /, &
     &        ' Total virial initial configuration: ', f12.5)
99002 FORMAT (' Total energy end of simulation    : ', f12.5, /, &
     &        '       running energy              : ', f12.5, /, &
     &        '       difference                  :  ', e12.5, /, &
     &        ' Total virial end of simulation    : ', f12.5, /, &
     &        '       running virial              : ', f12.5, /, &
     &        '       difference                  :  ', e12.5) 
99003 FORMAT (' Number of att. to displ. a part.  : ', i10, /, &
     &        ' success: ', i10, '(= ', f5.2, '%)')
      END
