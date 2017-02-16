      SUBROUTINE LAGR2(C)
*
*
*       Mass distribution for two mass groups.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 R2,ZMI(NMAX)
      COMMON/WORK1/ R2(NMAX)
      PARAMETER (LX=11)
      REAL*8 C(3),FLAGR(LX),RLAGR(LX)
*     DATA FLAGR/-1.9,-1.7,-1.5,-1.3,-1.1,-.9,-.7,-.5,-.3,-.1/
*     DATA FLAGR/0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,
*    &           0.75,0.9/
      DATA FLAGR/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.625,0.75,0.99/
*
*
*       Set square radii of single particles & c.m. bodies (NAME <= NZERO/5).
      NAM1 = NZERO/5
*       Apply the mass test in same proportion for two Plummer spheres.
      IF (KZ(5).EQ.2) THEN
          NAM1 = (NZERO - N1)/5
          NAM2 = N1 + (NZERO - N1)/5
      END IF
      ITER = 0
      NBLACK = 0
      DO 5 I = 1,N
          IF (KSTAR(I).EQ.14) NBLACK = NBLACK + 1
    5 CONTINUE
*
    1 NP = 0
      ZM1 = 0.0
      DO 10 I = IFIRST,NTOT
          IF (KZ(5).EQ.1) THEN
              IF (NAME(I).GT.NAM1.AND.NAME(I).LE.NZERO) GO TO 10
          ELSE IF (KZ(5).EQ.2) THEN
              IF ((NAME(I).GT.NAM1.AND.NAME(I).LE.N1).OR.
     &            (NAME(I).GT.NAM2.AND.NAME(I).LE.NZERO)) GO TO 10
          END IF
          ZM1 = ZM1 + BODY(I)
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
          ILIST(NP) = I
          ZMI(NP) = BODY(I)
   10 CONTINUE
*
*       Improve choice of NAM1 using relative deviation (max 10 tries).
      DM = (ZM1 - 0.5*ZMASS)/ZMASS
      IF (ABS(DM).GT.0.001.AND.ITER.LE.10) THEN
          ITER = ITER + 1
          IF (ABS(DM).GT.0.05) DM = DM*0.05/ABS(DM)
          NM = DM*FLOAT(N-NPAIRS)
          IF (NM.NE.0) THEN
              NAM1 = NAM1 - NM
              IF (KZ(5).EQ.2) THEN
                  NAM2 = NAM2 - NM*(NZERO - N1)/(5*NAM2)
              END IF
              GO TO 1
          END IF
      END IF
*
*       Obtain sum of inverse separations to evaluate mass segregation.
      POT1 = 0.0
      DO 14 L = 1,NP-1
          I = JLIST(L)
          DO 12 LL = L+1,NP
              J = JLIST(LL)
              RIJ2 = 0.0
              DO 11 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   11         CONTINUE
              POT1 = POT1 + 1.0/SQRT(RIJ2)
   12     CONTINUE
   14 CONTINUE
*       Define mean geometric distance.
      RP1 = FLOAT(NP)*(NP - 1)/(2.0*POT1)
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)

*       Determine Lagrangian radii for specified mass fractions.
      DO 20 IL = 1,LX
          ZM = 0.0
          ZMH = FLAGR(IL)*ZM1
          I = 0
   15     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
          IF (ZM.LT.ZMH) GO TO 15
          RLAGR(IL) = SQRT(R2(I))
   20 CONTINUE
*
*       Obtain half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5*ZM1
      I = 0
   25 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 25

      RH1 = SQRT(R2(I))
      NP1 = NP
*
      WRITE (31,30)  TIME+TOFF, (LOG10(RLAGR(K)),K=1,LX)
   30 FORMAT (' ',F8.1,13F7.3)
      CALL FLUSH(31)
*
*       Sort masses from first group in increasing order.
      CALL SORT1(NP,ZMI,ILIST)
      ZMX1 = ZMI(1)
      ZMX0 = 0.5*(ZMX1 + BODY1)
      ITRY = 0
   32 NP0 = 0
      LX0 = 0
*       Count dominant masses above half minimum + maximum single stars.
      DO 35 L = 1,NP
          IF (ZMI(L).GT.ZMX0) THEN
              IF (LX0.EQ.0) LX0 = L
              NP0 = NP0 + 1
          END IF
   35 CONTINUE
*       Ensure at least 1% in the most massive bin or use NBLACK > 2.
      ITRY = ITRY + 1
      NPCRIT = SQRT(FLOAT(N - NPAIRS))
      NPCRIT = MIN(NPCRIT,10)
      IF (NBLACK.GT.2) NPCRIT = NBLACK
      IF (NP0.LT.NPCRIT.AND.ITRY.LT.20) THEN
          ZMX0 = 0.9*ZMX0
          GO TO 32
      END IF
      IF (NP0.GT.NPCRIT.AND.ITRY.LT.20) THEN
          ZMX0 = 1.1*ZMX0
          GO TO 32
      END IF
*
*       Form sum of inverse separations to evaluate mass segregation.
      POT0 = 0.0
      LX0 = MAX(LX0,1)
      DO 44 L = LX0,NP-1
          I = ILIST(L)
          DO 42 LL = L+1,NP
              J = ILIST(LL)
              RIJ2 = 0.0
              DO 41 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   41         CONTINUE
              POT0 = POT0 + 1.0/SQRT(RIJ2)
   42     CONTINUE
   44 CONTINUE
*       Define mean geometric distance.
      IF (NP0.LE.1) POT0 = 1.0
      RP0 = FLOAT(NP0)*(NP0 - 1)/(2.0*POT0)
*
*       Search for single black holes and include KS binary.
      N14 = 0
      N14X = 0
      SEMI = 0
      POTBH = 0.0
      ZMBH = 0.0
      DO 45 I = IFIRST,NTOT
          IF (KSTAR(I).EQ.14.AND.BODY(I).GT.0.0D0) THEN
              N14 = N14 + 1
              ILIST(N14) = I
              ZMBH = ZMBH + BODY(I)*SMU
          ELSE IF (I.GT.N) THEN
*       Note: the following lines must be used for KS binary.
              J1 = 2*(I - N) - 1
              IF (KSTAR(J1).EQ.14.AND.KSTAR(J1+1).EQ.14) THEN
*       Add only one KS component to avoid small RIJ2 but use extra count.
                  N14 = N14 + 1
                  N14X = N14X + 1
                  ILIST(N14) = I
                  ZMBH = ZMBH + BODY(I)*SMU
                  IP = I - N
                  SEMI = -0.5*BODY(I)/H(IP)
              END IF
          END IF
   45 CONTINUE
*
*       Form sum of inverse separations to evaluate mass segregation.
      POT12 = 0.0
      RBH = 0.0
      DO 48 L = 1,N14-1
          I = ILIST(L)
          DO 47 LL = L+1,N14
              J = ILIST(LL)
              RIJ2 = 0.0
              DO 46 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   46         CONTINUE
              POTBH = POTBH + 1.0/SQRT(RIJ2)
*       Note internal terms not included at present.
              IF (I.LT.IFIRST) POT12 = POT12 + 1.0/SQRT(RIJ2)
   47     CONTINUE
   48 CONTINUE
*       Define mean geometric distance (modified by binary 10/2/11).
      N14 = N14 + N14X
      IF (N14.GT.1.AND.POTBH.GT.0.0D0) THEN
          RBH = FLOAT(N14)*(N14 - 1)/(2.0*POTBH)
      END IF
*       Omit internal binary contribution for secondary definition.
      IF (POT12.GT.0.0) THEN
          RX = FLOAT(N14-1)*(N14 - 2)/(2.0*(POTBH-POT12))
      ELSE
          RX = 0.0
      END IF
*
*       Search for single neutron stars and include KS binary.
      N13 = 0
      SEMI2 = 0
      POT13 = 0.0
      R13 = 0.0
      ZM13 = 0.0
      DO 145 I = IFIRST,NTOT
          IF (KSTAR(I).EQ.13) THEN
              N13 = N13 + 1
              ILIST(N13) = I
              ZM13 = ZM13 + BODY(I)*SMU
          ELSE IF (I.GT.N) THEN
              J1 = 2*(I - N) - 1
              IF (KSTAR(J1).EQ.13.AND.KSTAR(J1+1).EQ.13) THEN
                  N13 = N13 + 1
                  ILIST(N13) = I
                  ZM13 = ZM13 + BODY(I)*SMU
                  IP = I - N
                  SEMI2 = -0.5*BODY(I)/H(IP)
              END IF
          END IF
  145 CONTINUE
*
*       Form sum of inverse separations to evaluate mass segregation.
      DO 148 L = 1,N13-1
          I = ILIST(L)
          DO 147 LL = L+1,N13
              J = ILIST(LL)
              RIJ2 = 0.0
              DO 146 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
  146         CONTINUE
              POT13 = POT13 + 1.0/SQRT(RIJ2)
  147     CONTINUE
  148 CONTINUE
*       Define mean geometric distance (modified by binary 10/2/11).
      IF (N13.GT.1) R13 = FLOAT(N13)*(N13 - 1)/(2.0*POT13)
*
*       Treat the low-mass particles in the same way (exclude c.m. bodies).
      NP = 0
      ZM2 = 0.0
      DO 50 I = IFIRST,N
          IF (KZ(5).EQ.1) THEN
              IF (NAME(I).LE.NAM1) GO TO 50
          ELSE IF (KZ(5).EQ.2) THEN
              IF (NAME(I).LE.NAM1.OR.
     &           (NAME(I).GT.N1.AND.NAME(I).LE.NAM2)) GO TO 50
          END IF
          ZM2 = ZM2 + BODY(I)
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   50 CONTINUE
*
*       Treat second mass group in the same way.
      POT2 = 0.0
      DO 54 L = 1,NP-1
          I = JLIST(L)
          DO 52 LL = L+1,NP
              J = JLIST(LL)
              RIJ2 = 0.0
              DO 51 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   51         CONTINUE
              POT2 = POT2 + 1.0/SQRT(RIJ2)
   52     CONTINUE
   54 CONTINUE
      RP2 = FLOAT(NP)*(NP - 1)/(2.0*POT2)
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)

      DO 60 IL = 1,LX
          ZM = 0.0
          ZMH = FLAGR(IL)*ZM2
          I = 0
   55     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
          IF (ZM.LT.ZMH) GO TO 55
          RLAGR(IL) = SQRT(R2(I))
   60 CONTINUE
*
*       Obtain half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5*ZM2
      I = 0
   65 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 65
*
      RH2 = SQRT(R2(I))
      NP2 = NP
*
      WRITE (32,30)  TIME+TOFF, (LOG10(RLAGR(K)),K=1,LX)
      CALL FLUSH(32)
*
      WRITE (6,70)  TIME+TOFF, NP0, NP1, NP2, RH1, RH2, RP0, RP1, RP2
   70 FORMAT(/,' MASS GROUPS:    T NP0 NP1 NP2 RM1 RM2 RP0 RP1 RP2 ',
     &                           F7.1,I5,I6,I7,1X,5F7.3)
*
      IF (N14.GT.1.OR.N13.GT.1) THEN
          N14 = MAX(N14,1)
          WRITE (6,80)  TIME+TOFF, TPHYS, N14, ZMBH/FLOAT(N14), RBH,
     &                  SEMI, N13, R13, RX
   80     FORMAT (/,' BH/NS SUBSYSTEM:    T TPHYS NBH <MBH> RBH A ',
     &                                    'NS R13 RX ',
     &                       2F7.1,I4,F6.1,F7.3,1P,E10.2,0P,I5,2F7.3)
      END IF
*
      RETURN
*
      END
