      SUBROUTINE HRPLOT(luminosities,teff, kstara)
*
*
*       HR diagnostics of evolving stars.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  M0,M1,M2,LUM,LUM2,MC,ME,K2
*
      real*4 luminosities(nmax), teff(nmax)
      integer kstara(nmax)
*
      TPLOT = TIME
    1 FORMAT (I8,F12.4,F12.4)
      NS = N - 2*NPAIRS
      IMERGE = 0
*
      DO 20 I = 1,N
          M0 = BODY0(I)*ZMBAR
          M1 = BODY(I)*ZMBAR
*       Replace ghost mass of single star with original value from merged KS.
          IF (M1.EQ.0.0.AND.I.GE.IFIRST) THEN
              IM = 0
*       Search merger table for ghost to identify corresponding index.
              DO 2 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(I)) THEN
                      IM = K
                  END IF
    2         CONTINUE
*       Skip any ghosts associated with chain regularization.
              IF (IM.EQ.0) THEN
                  WRITE (6,3)  I, NCH
    3             FORMAT (' WARNING!    HRPLOT   I NCH ',I6,I4)
                  GO TO 20
              END IF
              M1 = CM(2,IM)*ZMBAR
          END IF
*
*       Obtain stellar parameters at current epoch.
          KW = KSTAR(I)
          AGE = MAX(TPLOT,TEV0(I))*TSTAR - EPOCH(I)
          CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
          CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                RM,LUM,KW,MC,RCC,ME,RE,K2)
          IF (I.LT.IFIRST) THEN
              JPAIR = KVEC(I)
              J2 = 2*JPAIR
              IF (I.EQ.J2) GO TO 20
              J1 = 2*JPAIR - 1
              ICM = N + JPAIR
              M2 = BODY(J2)*ZMBAR
              RI = (X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2 +
     &                                        (X(3,ICM) - RDENS(3))**2
*       Check for ghost binary.
              IF (M1.EQ.0.0) THEN
                  IM = 0
*       Search merger table to identify corresponding index of c.m. name.
                  DO 4 K = 1,NMERGE
                      IF (NAMEM(K).EQ.NAME(ICM)) THEN
                          IM = K
                      END IF
    4             CONTINUE
                  IF (IM.EQ.0) GO TO 20
*       Copy masses and obtain evolution parameters for first component.
                  M1 = CM(3,IM)*ZMBAR
                  M2 = CM(4,IM)*ZMBAR
                  KW = KSTAR(J1)
                  AGE = MAX(TPLOT,TEV0(J1))*TSTAR - EPOCH(J1)
                  CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
                  CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                        RM,LUM,KW,MC,RCC,ME,RE,K2)
              END IF
              RJ = R(JPAIR)
              HJ = H(JPAIR)
*       Determine merger & ghost index for negative c.m. name.
              IF (NAME(N+JPAIR).LT.0) THEN
                  CALL FINDJ(J1,J,IMERGE)
*       Skip second binary of quadruple.
                  IF (NAME(J).GT.NZERO) GO TO 20
                  M1 = CM(1,IMERGE)*ZMBAR
                  M2 = CM(2,IMERGE)*ZMBAR
                  HJ = HM(IMERGE)
                  RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                          XREL(3,IMERGE)**2)
*       Re-define index of second component and obtain parameters of M1.
                  J2 = J
                  AGE = MAX(TPLOT,TEV0(J1))*TSTAR - EPOCH(J1)
                  KW = KSTAR(J1)
                  CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
                  CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                        RM,LUM,KW,MC,RCC,ME,RE,K2)
              END IF
              M0 = BODY0(J2)*ZMBAR
              KW2 = KSTAR(J2)
              AGE = MAX(TPLOT,TEV0(J2))*TSTAR - EPOCH(J2)
              CALL STAR(KW2,M0,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
              CALL HRDIAG(M0,AGE,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                    RM2,LUM2,KW2,MC,RCC,ME,RE,K2) 
*     Mark: remove division by RC
              RI = SQRT(RI)
*       Specify relevant binary mass.
              IF (BODY(J1).GT.0.0D0) THEN
                  BODYI = (M1 + M2)/ZMBAR
              ELSE
                  BODYI = CM(3,IMERGE) + CM(4,IMERGE)
              END IF
              SEMI = -0.5*BODYI/HJ
              ECC2 = (1.0 - RJ/SEMI)**2
              ECC = SQRT(ECC2)
              PB = DAYS*SEMI*SQRT(ABS(SEMI)/BODYI)
              PB = MIN(PB,99999.9D0)
              PB = LOG10(ABS(PB))
              SEMI = LOG10(ABS(SEMI*SU))
              R1 = LOG10(RM)
              R2 = LOG10(RM2)
              ZL1 = LOG10(LUM)
              ZL2 = LOG10(LUM2)
              TE1 = 0.25*(ZL1 - 2.0*R1) + 3.7
              TE2 = 0.25*(ZL2 - 2.0*R2) + 3.7
              luminosities(j1) = zl1
              luminosities(j2) = zl2
              luminosities(icm) = log10(10**zl1 + 10**zl2)

              teff(j1) = te1
              teff(j2) = te2
              teff(icm) = 0.0

              kstara(j1) = kw
              kstara(j2) = kw2
              kstara(icm) = kstar(icm)

    5         FORMAT (2I7,2I3,I4,6F10.4,F9.2,F10.5,10F8.3)
          ELSE
*       Create output file for single stars (skip chain subsystem or ghost).
              IF (NAME(I).EQ.0.OR.BODY(I).EQ.0.0D0) GO TO 20
              RI = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
*     Mark: remove division by RC and max(Ri,99)
              RI = SQRT(RI)
*              RI = MIN(RI,99.0D0)
              R1 = LOG10(RM)
              ZL1 = LOG10(LUM)
*       Form LOG(Te) using L = 4*pi*R**2*\sigma*T**4 and solar value 3.7.
              TE = 0.25*(ZL1 - 2.0*R1) + 3.7

              luminosities(i) = zl1 
              teff(i) = TE
              kstara(i) = kw
   10         FORMAT (I10,I4,6F12.4,5F10.3)
          END IF
   20 CONTINUE
*
*       Update next plotting time.
      TPLOT = TPLOT + DTPLOT
      CALL FLUSH(82)
      CALL FLUSH(83)
*
      RETURN
*
      END
