      SUBROUTINE FDYNFRI(XI,XIDOT,FM,FD)
*
*
*       Jpetts - obtain dynamical friction force
*       Intent in - FM = Acceleration, FD = Jerk
*                   
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8  XI(3),FM(3),FD(3), XGDD(3), FS(3), ACC(3), ACCD(3)
      REAL*8  V2, VVD, H3, VVDD, VDVD, XGDDD(3), XGD(3), XG(3)
      REAL*8  XIDOT(3)
*
*     Get galactocentric position and velocity
      DO 5 K = 1,3
        XG(K) = RG(K) + XI(K)
        XGD(K) = VG(K) + XIDOT(K)
    5 CONTINUE

      V2 = XGD(1)**2.0 + XGD(2)**2.0 + XGD(3)**2.0
      H3 = DYNFCOEF/(V2*SQRT(V2))
*
*     Calculate VDOT contribution from the potential
      DO 6 K = 1,3
        XGDD(K) = 0.d0
        ACC(K) = 0.d0
        ACCD(K) = 0.d0
    6 CONTINUE
      CALL FNUC(XG,XGD,ACC,ACCD)
      DO 61 K = 1,3
        XGDD(K) = XGDD(K) + ACC(K)
   61 CONTINUE
      CALL FBULGE(XG,XGD,ACC,ACCD)
      DO 62 K = 1,3
        XGDD(K) = XGDD(K) + ACC(K)
   62 CONTINUE
*
*
*     Calculate the acceleration and add it to XGDD
      DO 30 K = 1,3
        FM(K) = -H3*XGD(K)
        XGDD(K) = XGDD(K) + FM(K) !XGDD due to potential + df
   30 CONTINUE

      VVD = (XGD(1)*XGDD(1)+XGD(2)*XGDD(2)+XGD(3)*XGDD(3))

*     Calculate the jerk
      DO 40 K = 1,3
        FD(K) = -H3*(XGDD(K) - 3.0*VVD*XGD(K)/V2)
   40 CONTINUE
*
      RETURN
*
      END

