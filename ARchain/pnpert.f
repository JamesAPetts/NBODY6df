      SUBROUTINE PNPERT(X,V,R,RDOTV,V2,M1,M2,ACC,DVGR,SPINA,DSPIN)
*
*
*       Post-Newtonian terms.
*       ---------------------
*
*       L. Blanchet & B. Iyer, Class. Quantum Grav. 20 (2003), 755.
*       T. Mora & Clifford Will, gr-qc/0312082.
*       ---------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/POSTN/  CLIGHT,TAUGR,RZ,GAMMA,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
      REAL*8  ACC(3),X(3),V(3),DVGR(3),SPINA(3),DSPIN(3),DVQ(3)
      SAVE IT,IT2,PI2
      DATA IT,IT2 /0,0/
*
*
      IF (IT2.EQ.0) THEN
          PI = 4.0*ATAN(1.0D0)
          PI2 = PI**2
          IT2 = 1
      END IF
*       Set total mass and reduced mass parameter.
      M = M1 + M2
      ETA = M1*M2/M**2
      ETA2 = ETA**2
      ETA3 = ETA2*ETA
*
      R2 = R**2
      MR = M/R
      RD = RDOTV/R
*
*       Form the post-Newtonian scalar A & B-terms (skip A2/B2 if IPN < 2).
      IF (IPN.GE.1) THEN
          A1 = 2.0*(2.0 + ETA)*MR - (1.0 + 3.0*ETA)*V2 + 1.5*ETA*RD**2
          B1 = 2.0*(2.0 - ETA)*RD
      ELSE
          A1 = 0.0
          B1 = 0.0
      END IF
*
      IF (IPN.GE.2) THEN
          A2 = - 0.75*(12.0 + 29.0*ETA)*MR**2 - ETA*(3.0-4.0*ETA)*V2**2
     &         - 15.0/8.0*ETA*(1.0 - 3.0*ETA)*RD**4
     &         + 0.5*ETA*(13.0 - 4.0*ETA)*MR*V2
     &         + (2.0 + 25.0*ETA + 2.0*ETA**2)*MR*RD**2
     &         + 1.5*ETA*(3.0 - 4.0*ETA)*V2*RD**2
          B2 = - 0.5*RD*((4.0 + 41.0*ETA + 8.0*ETA**2)*MR
     &         - ETA*(15.0 +4.0*ETA)*V2 + 3.0*ETA*(3.0 + 2.0*ETA)*RD**2)
          A2 = A2/CLIGHT**2
          B2 = B2/CLIGHT**2
      ELSE
          A2 = 0.0
          B2 = 0.0
      END IF
*
      IF (IPN.GE.3) THEN
          A3 = (16. + (1399./12. - 41./16.*PI2)*ETA + 35.5*ETA2)*MR**3
     &         + ETA*(20827./840. + 123./64.*PI2 - ETA2)*MR**2*V2
     &         - (1. + (22717./168. + 615./64.*PI2)*ETA + 11./8.*ETA2
     &         - 7.*ETA3)*MR**2*RD**2
     &         - 0.25*ETA*(11. - 49.*ETA + 52.*ETA2)*V2**3
     &         + 35./16.*ETA*(1. -5.*ETA + 5.*ETA2)*RD**6
     &         - 0.25*ETA*(75. + 32.*ETA - 40.*ETA2)*MR*V2**2
     &         - 0.5*ETA*(158. - 69.*ETA - 60.*ETA2)*MR*RD**4
     &         + ETA*(121. - 16.*ETA - 20.*ETA2)*MR*V2*RD**2
     &         + 0.375*ETA*(20. - 79.*ETA + 60.*ETA2)*V2**2*RD**2
     &         - 15./8.*ETA*(4. - 18.*ETA + 17.*ETA2)*V2*RD**4
          B3 = (4. + (5849./840. + 123./32.*PI2)*ETA - 25.*ETA2
     &         - 8.*ETA3)*MR**2
     &         + 0.125*ETA*(65. - 152.*ETA - 48.*ETA2)*V2**2
     &         + 15./8.*ETA*(3. - 8.*ETA - 2.*ETA2)*RD**4
     &         + ETA*(15. + 27.*ETA + 10.*ETA2)*MR*V2
     &         - 1./6.*ETA*(329. + 177.*ETA + 108.*ETA2)*MR*RD**2
     &         - 0.75*ETA*(16. - 37.*ETA - 16.*ETA2)*V2*RD**2
          B3 = B3*RD
          A3 = A3/CLIGHT**4
          B3 = B3/CLIGHT**4
      ELSE
          A3 = 0.0
          B3 = 0.0
      END IF
*
      IF (IPN.GE.4) THEN
*       Note: copied from Seppo Mikkola code 3/1/2012 with rd = vr.
          A4=-8./5*eta*(m/r)*rd*(23./14*(43+14*eta)*(m/r)**2
     &       +3./28*(61+70*eta)*vv**2
     &       +70*rd**4+1./42*(519-1267*eta)*(m/r)*vv
     &       +.25d0*(147+188*eta)*(m/r)*rd**2-15/4.*(19+2*eta)*vv*rd**2)
          B4=8./5.*eta*(m/r)*(1./42.*(1325+546*eta)*(m/r)**2
     &       +1./28.*(313+42*eta)*vv**2+75*rd**4
     &       -1./42.*(205+777*eta)*(m/r)*vv
     &       +1./12.*(205+424*eta)*(m/r)*rd**2-.75*(113+2*eta)*vv*rd**2)
          A4 = A4/CLIGHT**5
          B4 = B4/CLIGHT**5
      ELSE
          A4 = 0.0
          B4 = 0.0
      END IF
*
*       Include PN2.5 and add all contributions.
      A25 = 8.0/5.0*ETA*MR*RD*(17.0/3.0*MR + 3.0*V2)
      A25 = A25/CLIGHT**3
      A = A1 + A2 + A25 + A3 + A4
*
      B25 = -8.0/5.0*ETA*MR*(3.0*MR + V2)
      B25 = B25/CLIGHT**3
      B = B1 + B2 + B25 + B3 + B4
*
      IF (ISPIN.GT.0) then
          MX = MAX(M1,M2)
          MY = MIN(M1,M2)
         call gopu_SpinTerms(X,V,r,MX,MY,CLIGHT,spina,dvq,dspin)
      else
         do k=1,3
         dvq(k)=0.0
         dspin(k)=0.0
         end do
      end if
*
*       Note scaling by M/R2 with individual masses used after RETURN.
      GMC = 1.0/(CLIGHT**2*R2)
*       Form PN contributions (negative sign for equations of motion).
      DO 20 K = 1,3
*       Note minus sign removed 6/2010 consistent with Seppo's new convention.
*         ACC(K) = GMC*(A*X(K)/R + B*V(K)) + DVQ(K)*GMC
          ACC(K) = GMC*(A*X(K)/R + B*V(K)) + DVQ(K)/M
*         DVGR(K) = GMC*(A25*X(K)/R + B25*V(K))
*       Include all GR terms instead of just radiation in EnerGR (4/11/09).
          DVGR(K) = ACC(K)
   20 CONTINUE
*
*     IT = IT + 1
*     WRITE (6,22)  IT, R, ACC
*  22 FORMAT (' PNPERT    IT R ACC  ',I9,1P,4E10.2)
*       Obtain the relative perturbation for diagnostics.
*     DF = SQRT(ACC(1)**2 + ACC(2)**2 + ACC(3)**2)
*     GAMMA = DF*R2/M
*
      RETURN
      END
