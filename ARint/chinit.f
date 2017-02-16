      SUBROUTINE CHINIT(ISUB)
*
*
*       Initialization of chain system.
*       -------------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC
        PARAMETER (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/ECHAIN/  ECH
      COMMON/POSTN/  CVEL,TAUGR,RZ1,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
      REAL*8  ANG(3),FIRR(3),FD(3)
      LOGICAL FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Initialize GR parameters and current time (in case #28 = 0).
      IF (FIRST) THEN
          CVEL = 3.0D+05/VSTAR
          CLIGHT = CVEL
          FIRST = .FALSE.
      END IF
      KZ24 = KZ(24)
      TAUGR = 1.3D+18*RAU**4/SMU**3
      TSP = TOFF + TIME
*
*       Define chain membership.
      CALL SETSYS
*
*     DO 15 L = 1,NCH
*         J = JLIST(L)
*         WRITE (6,1)  J, NAME(J), KSTAR(J), BODY(J)*SMU, RADIUS(J)*SU
*   1     FORMAT (' INITIAL MEMBER    J NM K* M R* ',2I6,I4,2F7.1)
*  15 CONTINUE
*     CALL FLUSH(6)
*
*       Initialize c.m. variables.
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      BX = 0.0
      DO 4 L = 1,NCH
          J = JLIST(L)
          SIZE(L) = RADIUS(J)
          ISTAR(L) = KSTAR(J)
*       Place the system in first single particle locations.
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X4(K,L) = X(K,J)
              XDOT4(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
          IF (M(L).GT.BX) THEN
              BX = M(L)
              LX = L
          END IF
    4 CONTINUE
*
*       Produce output for relevant parameters (RZ = 3*Schwarzschild).
      IF (CLIGHT.GT.0.0D0) THEN
          RZ1 = 6.0*CM(7)/CLIGHT**2
      ELSE
          RZ1 = 0.0
      END IF
      L1 = NCH + 1 - LX
      RCOLL = (CM(7)/M(L1))**0.3333*SIZE(L1)
*
*       Set c.m. coordinates & velocities of subsystem.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for chain regularization.
      LK = 0
      DO 8 L = 1,NCH
          DO 7 K = 1,3
              LK = LK + 1
              X4(K,L) = X4(K,L) - CM(K)
              XDOT4(K,L) = XDOT4(K,L) - CM(K+3)
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and and save in chain energy.
      CALL CONST(XCH,VCH,M,NCH,ENERGY,ANG,GAM)
      ECH = ENERGY
*
*       Find sum of mass products and individual separations (for CHLIST).
      SUM = 0.0D0
      RSUM = 0.0D0
      DO 10 L = 1,NCH-1
          DO 9 K = L+1,NCH
              SUM = SUM + M(L)*M(K)
              RLK2 = (X4(1,L) - X4(1,K))**2 + (X4(2,L) - X4(2,K))**2 +
     &                                        (X4(3,L) - X4(3,K))**2
              RSUM = RSUM + SQRT(RLK2)
    9     CONTINUE
   10 CONTINUE
*
*       Reduce RSUM by geometrical factor and check upper limit from IMPACT.
      IF (NCH.EQ.4) RSUM = 0.5*RSUM
      RSUM = MIN(FLOAT(NCH-1)*RSUM/FLOAT(NCH),RMIN)
*
*       Define gravitational radius for initial perturber list.
      RGRAV = SUM/ABS(ENERGY)
*
*       Avoid small value after collision (CHTERM improves perturbers).
      IF (NCH.GT.2) THEN
          RGRAV = MIN(RGRAV,0.5*RSUM)
      END IF
*
*       Set global index of c.m. body and save name (SUBSYS sets NAME = 0).
      IF (TIMEC.GT.0.0D0) ICH0 = ICH
      ICH = JLIST(LX)
      NAME0 = NAME(ICH)
*
*       Define subsystem indicator (ISYS = 1,2,3,4 for triple, quad, chain).
      ISYS(NSUB+1) = 4
*
*       Form ghosts and initialize c.m. motion in ICOMP = JLIST(NCH).
      CALL SUBSYS(NCH,CM)
*
*       Copy neighbour list for ghost removal.
      NNB = LIST(1,ICH)
      DO 20 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   20 CONTINUE
*
*       Remove any ghosts from neighbour list.
      CALL NBREM(ICH,NCH,NNB)
*
*       Initialize perturber list for integration of chain c.m.
      CALL CHLIST(ICH)
*
*       Perform differential F & FDOT corrections due to perturbers.
      DO 25 K = 1,3
          FIRR(K) = 0.0D0
          FD(K) = 0.0
   25 CONTINUE
      CALL CHFIRR(ICH,2,X(1,ICH),XDOT(1,ICH),FIRR,FD)
      DO 30 K = 1,3
          F(K,ICH) = F(K,ICH) + 0.5*FIRR(K)
          FDOT(K,ICH) = FDOT(K,ICH) + ONE6*FD(K)
   30 CONTINUE
*
*       Take maximum integration interval equal to c.m. step.
      TMAX = STEP(ICH)
*
*       Check next treatment time of perturbers & output time.
      CALL TCHAIN(NSUB,TSMIN)
      TMAX = MIN(TMAX,TSMIN)
      TMAX = MIN(TMAX,TADJ - TIME)
*
*       Copy total energy and output & capture option for routine CHAIN.
      CM(8) = BE(3)
      KZ27 = KZ(27)
      KZ30 = KZ(30)
*       Copy velocity scale factor to VSTAR1.
      VSTAR1 = VSTAR
*
*       Assign new subsystem index and begin chain regularization.
      ISUB = NSUB
      NCHAIN = NCHAIN + 1
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
