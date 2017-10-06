      SUBROUTINE FDYNFRI(XGDOT,ID,FM,FD)
*
*
*       Jpetts - obtain dynamical friction force
*       Intent in - XGDOT = VG, ID = ID of star, FM = Cummulative force,
*                   FD = Cummulative derivative of force
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8  XI(3),XGDOT(3),FM(3),FD(3)
      REAL*8  V2, VVD, H3
*
*       Calculate force and derivative
*      
      V2 = XGDOT(1)**2.0 + XGDOT(2)**2.0 + XGDOT(3)**2.0
      H3 = -DYNFCOEF/(V2*SQRT(V2))
*
      DO 30 K = 1,3
        FM(K) = -H3*XGDOT(K)        
   30 CONTINUE
      
      VVD = (FM(1)**2 + FM(2)**2 + FM(3)**2)**0.5

      DO 40 K = 1,3
         FD(K) = -H3*(FM(K) - 3.0*VVD*XGDOT(K)/(V2**0.5))
   40 CONTINUE
*

*      Store delta energy removed in DEDYNFRI(ID)      
       DEDYNFRI(ID) = 0.5D0*BODY(ID)*(FM(1)**2 + FM(2)**2 + FM(3)**2)
      RETURN
*
      END

