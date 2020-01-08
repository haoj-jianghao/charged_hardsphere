      SUBROUTINE MCMOVE(En, En_i, Vir, Attempt, Nacc1, Nacc2, Nacc3, Dr1, Dr2, Dr3, Iseed)


      USE parameters
      USE conf
      USE system

      IMPLICIT NONE

      DOUBLE PRECISION ::  RANF
      DOUBLE PRECISION :: enn, eno, En, Vir, viro, virn , enohs, ennhs,enoi,enni, En_i
      DOUBLE PRECISION ::  Dr1, Dr2, Dr3
      INTEGER :: o, Attempt, Nacc1, Nacc2, Nacc3, jb, Iseed
      DOUBLE PRECISION :: lambda
      
!      USE conf
      
      DOUBLE PRECISION ::  x1(Npart), y1(Npart), z1(Npart), xn,yn,zn
      DOUBLE PRECISION ::  x2(Npart), y2(Npart), z2(Npart), xnh,ynh,znh
      DOUBLE PRECISION ::  x3(Npart), y3(Npart), z3(Npart), xni,yni,zni
      
      X1 = X
      Y1 = Y
      Z1 = Z
      
      X2 = Xh
      Y2 = Yh
      Z2 = Zh
      
      X3 = Xi
      Y3 = Yi
      Z3 = Zi
      
      
      !------------------------------------------------
      Attempt = Attempt + 1
      jb = 1
!     ---select a particle at random
      o = INT(NPART*RANF(Iseed)) + 1
!     ---calculate energy old configuration
      CALL ENERI(X, Y, Z, eno, viro)
		

!     ---give particle a random displacement
      x1(o) = X(o) + (RANF(Iseed)-0.5D0)*Dr1
      y1(o) = Y(o) + (RANF(Iseed)-0.5D0)*Dr1
      z1(o) = Z(o) + (RANF(Iseed)-0.5D0)*Dr1
      
      xn = X1(o) 
      yn = Y1(o) 
      zn = Z1(o) 
!     ---calculate energy new configuration:
     
      CALL ENERI(x1, y1, z1, enn, virn)

!	write(12,*) enn, eno

!     ---acceptance test
      IF (RANF(Iseed).LT.EXP(-BETA*(enn-eno))) THEN
!      write (10,*) X(o),xn
!        --accepted
         Nacc1 = Nacc1 + 1
!	write (12,*) En
         En = En + (enn-eno)
!	write (13,*) En
         Vir = Vir + (virn-viro)
!        ---put particle in simulation box
         IF (xn.LT.0) xn = xn + BOX
         IF (xn.GT.BOX) xn = xn - BOX
         IF (yn.LT.0) yn = yn + BOX
         IF (yn.GT.BOX) yn = yn - BOX
         IF (zn.LT.0) zn = zn + BOX
         IF (zn.GT.BOX) zn = zn - BOX
         X(o) = xn
         Y(o) = yn
         Z(o) = zn
      ELSE
 !       write(*,*) "rejected"
      END IF
     
     
     !------------------------------------------------
     ! a parallel monte carlo of hard sphere for BAR
      
      o = INT(NPART*RANF(Iseed)) + 1
      CALL ENERI_hs(Xh, Yh, Zh, enohs)
      x2(o) = Xh(o) + (RANF(Iseed)-0.5D0)*Dr2
      y2(o) = Yh(o) + (RANF(Iseed)-0.5D0)*Dr2
      z2(o) = Zh(o) + (RANF(Iseed)-0.5D0)*Dr2
      xnh = X2(o) 
      ynh = Y2(o) 
      znh = Z2(o) 
      CALL ENERI_hs(x2, y2, z2, ennhs)
      IF (RANF(Iseed).LT.EXP(-BETA*(ennhs-enohs))) THEN
      Nacc2 = Nacc2 + 1
         IF (xnh.LT.0) xnh = xnh + BOX
         IF (xnh.GT.BOX) xnh = xnh - BOX
         IF (ynh.LT.0) ynh = ynh + BOX
         IF (ynh.GT.BOX) ynh = ynh - BOX
         IF (znh.LT.0) znh = znh + BOX
         IF (znh.GT.BOX) znh = znh - BOX
         Xh(o) = xnh
         Yh(o) = ynh
         Zh(o) = znh
      END IF
      
      !------------------------------------------------
      ! a parallel monte carlo of intermdeida potential for BAR U_inter = U0 + lambda*(U1-U0)
      o = INT(NPART*RANF(Iseed)) + 1
      lambda = 0.5
      CALL ENERI_inter(Xi, Yi, Zi, enoi)
      x3(o) = Xi(o) + (RANF(Iseed)-0.5D0)*Dr3
      y3(o) = Yi(o) + (RANF(Iseed)-0.5D0)*Dr3
      z3(o) = Zi(o) + (RANF(Iseed)-0.5D0)*Dr3
      xni = X3(o) 
      yni = Y3(o) 
      zni = Z3(o)
      CALL ENERI_inter(x3, y3, z3, enni)
      IF (RANF(Iseed).LT.EXP(-BETA*(enni-enoi))) THEN
      Nacc3 = Nacc3 + 1
      En_i = En_i + (enni-enoi)
         IF (xni.LT.0) xni = xni + BOX
         IF (xni.GT.BOX) xni = xni - BOX
         IF (yni.LT.0) yni = yni + BOX
         IF (yni.GT.BOX) yni = yni - BOX
         IF (zni.LT.0) zni = zni + BOX
         IF (zni.GT.BOX) zni = zni - BOX
         Xi(o) = xni
         Yi(o) = yni
         Zi(o) = zni
      END IF
      

      END SUBROUTINE MCMOVE