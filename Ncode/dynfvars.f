      SUBROUTINE DYNFVARS
*
*       Jpetts - Calculate DYNFCOEF et al
*       Recalculate variables used for dynamical friction formula every
*       time parameters are adjusted. Note if one changes the potential
*       one must also update the calculations of BDENS, VDIST and bmax.
*       One must note that the Maxwellian assumption for the velocity
*       distribution underestimates the frictional force for a power law
*       slope shallower than gamma = 0.5. The stalling is also not
*       not correctly captured in constant density cores. A dynamical
*       friction model for cores will be added in the future
*       (see Petts et al 2015 & Petts et al 2016a).
*       ----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8 XRA,XRA2,XRA3,XRA4,MASSENC
      REAL*8 RG_RT,DRG_RT,BDENS_RT1,BDENS_RT2,RT_D2PHIDR2
      REAL*8 RT_OMEGA2,RDIST_RT,RTIDE_TEMP, FVSTAR, E_APO, L2_APO, RG_APO
      REAL*8 VESC, FV_VS, FV_ESC
      INTEGER*4 RT_CHECK
*
*
*        Calculate r
      DO 10 K = 1,3
          RGDENS(K) = RG(K) + RDENS(K)
   10 CONTINUE   
      RGDENSMAG = SQRT(RGDENS(1)**2.0+RGDENS(2)**2.0+RGDENS(3)**2.0)
*
*       Calculate v
      DO 101 K = 1,3
          VGCORE(K) = VG(K) + VCORE(K)
  101 CONTINUE 
      VGCOREMAG = SQRT(VGCORE(1)**2.0+VGCORE(2)**2.0+VGCORE(3)**2.0)

*       Calculate the coulomb logarithm.
*     ----Calculate bmax = min((rho/[drho/dr]),r)
      COULOG = min(RGDENSMAG*(RGDENSMAG+AR)/(AR*GAM + 4.d0*RGDENSMAG)
     &     ,RGDENSMAG)
*     ----Calculate bmin = max(R_hm, G*M_s/V_s**2)
      COULOG = COULOG/max(RSCALE,MASSCL/VGCOREMAG**2)
      COULOG = log(COULOG)
*
*       Calculate the local background density.
      BDENS = ((3.0-GAM)*GMB*AR)/(4.0*ONEPI)
      BDENS = BDENS/(RGDENSMAG**GAM)
      BDENS = BDENS/((RGDENSMAG+AR)**(4.0-GAM))

***********************Calculate f(v*<vs)******************************
* If gamma = 0.0, integrate the self-consistent velocity distribution,
* as the Maxwellian approximation of the velocity distribution is inaccurate
* at ~< 0.5 a. for gam >= 0.5 the Maxwellian approximation is accurate
* within ~10%, and is employed for speed of calculation.
*
********Calculate the local velocity dispersion*************************
*
*  Included here is the general velocity disperion profile for a Dehnen
*  model with (or without) a central point mass black hole. The 
*  derivation is general in that it allows one to use any SMBH/
*  background mass ratio, any slope gamma and any scale radius. The 
*  functions f(x,gam) and h(x,gam) must be pre-calculated and included 
*  in the code. gam = 0.5,1.0, 1.5, 1.75 and 2.0 are included by 
*  default. Adding another value of gamma is straightforward and described
*  in the README.
*
*  The presence of a SMBH is optional, setting GMG to 0.0 in the
*  input file will ommit a SMBH.
*
*  One needs to derive the velocity disperion profile of the
*  background distribution themselves if they want to use a density
*  distribution other than a Dehnen model. The functions f(x,gam) and h(x,gam)
*  are only analytic when 4*gamma is an integer.
*  
*
*  V_r^2(r) = ((GMB/a)*(x+1)^(4-gam)) * (f(x,gam)  + mu*h(x,gam))
*
*  where:       x = r/a,                 f(x,gam) = BVDISPBULGE
*               h(x,gam) = BVDISPBH      mu = GMG/GMB
*
************************************************************************       
      XRA = RGDENSMAG/AR
      XRA2 = XRA**2.d0
      XRA3 = XRA**3.d0
      XRA4 = XRA**4.d0

*
*     ----select apropriate calculated pre-calculated integrals
      IF(GAM == 0.0) THEN
         IF (GMG .ne. 0.0) THEN
            WRITE(6,*) "No isotropic solution for gamma = 0.0 + BH."
            WRITE(6,*) "See Antonnini+ 2012"
            STOP
         END IF
         !Calculate the self consistent f(v*<vs)
         GO TO 4

      ELSE IF(GAM == 0.25) THEN
         IF (GMG .ne. 0.0) THEN
            WRITE(6,*) "No isotropic solution for gamma = 0.25 + BH."
            WRITE(6,*) "See Antonnini+ 2012"
            STOP
         END IF

         BVDISPBULGE = 0.0005772005772d0*XRA**(3.d0/4.d0) * (128.d0*XRA4
     &        + 704.d0*XRA3 + 1584.d0*XRA2 + 1848.d0*XRA + 1155.d0)
     &        *(XRA/(XRA+1.d0))**(3.d0/4.d0)
     &        /(XRA+1.d0)**(19.d0/4.d0)


         BVDISPBH = 0.0


      ELSE IF(GAM == 0.5) THEN
          BVDISPBULGE = 0.2d0/(1.d0+XRA)**5.d0
*     
*
          BVDISPBH = (0.1333333333d0*(128.d0*XRA**(3.d0/2.d0)*(1.d0+XRA)
     &    **(5.d0/2.d0)-128.d0*XRA4-320.d0*XRA3-240.d0*XRA2-40.d0*XRA
     &    +5.d0))/(XRA**(3.d0/2.d0)*(1.d0+XRA)**(5.d0/2.d0))
*
*     ----multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE IF(GAM == 1.0) THEN
          BVDISPBULGE = (1.d0/12.d0)*(-12.d0*XRA3-25.d0+12.d0*
     &    log(XRA+1.d0)-12.d0*log(XRA)-52.d0*XRA-42.d0*XRA2+12.d0*
     &    log(XRA+1.d0)*XRA4-12.d0*log(XRA)*XRA4+48.d0*log(XRA+1.d0)
     &    *XRA3-48.d0*log(XRA)*XRA3+72.d0*log(XRA+1.d0)*XRA2-72.d0*
     &    log(XRA)*XRA2+48.d0*log(XRA+1.d0)*XRA-48.d0*log(XRA)*XRA)/
     &    (XRA4+4.d0*XRA3+6*XRA2+4.d0*XRA+1.d0)
*     
*
          BVDISPBH = -(1.d0/2.d0)*(12.d0*log(XRA)*XRA4-12.d0*
     &    log(XRA+1.d0)*XRA4+24.d0*log(XRA)*XRA3-24.d0*log(XRA+1.d0)*
     &    XRA3+12.d0*log(XRA)*XRA2-12.d0*log(XRA+1.d0)*XRA2+12.d0*XRA3+
     &    18.d0*XRA2+4.d0*XRA-1.d0)/(XRA2*(XRA2+2.d0*XRA+1.d0))
*
*     ----multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE IF(GAM == 1.5) THEN
          BVDISPBULGE = -(0.3333333333d0*(-3.d0-12.d0*XRA3+12.d0*
     &    log(1.d0+XRA)*XRA4-12.d0*log(XRA)*XRA4+36.d0*log(1.d0+XRA)
     &    *XRA3-36.d0*log(XRA)*XRA3+36.d0*log(1.d0+XRA)*XRA2-36.d0*
     &    log(XRA)*XRA2+12.d0*log(1.d0+XRA)*XRA-12.d0*log(XRA)*XRA-
     &    30.d0*XRA2-22.d0*XRA))/((1.d0+XRA)**3*XRA)
*     
*
          BVDISPBH = -(0.1333333333d0*(128.d0*XRA**(5.d0/2.d0)*
     &    (1.d0+XRA)**(3.d0/2.d0)-128.d0*XRA4-192.d0*XRA3-48.d0*XRA2+
     &    8.d0*XRA-3.d0))/(XRA**(5.d0/2.d0)*(1.d0+XRA)**(3.d0/2.d0))
*
*     ----multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE IF(GAM == 1.75) THEN
          BVDISPBULGE = (0.1333333333d0*(128.d0*XRA**(3.d0/4.d0)*
     &    (1.d0+XRA)**(7.d0/4.d0)*(XRA*(1.d0+XRA))**(3.d0/4.d0)-128.d0
     &    *XRA4-320.d0*XRA3-240.d0*XRA2-40.d0*XRA+5.d0))/
     &    (XRA**(3.d0/4.d0)*(1.d0+XRA)**(7.d0/4.d0)*(XRA*(1.d0+XRA))
     &    **(3.d0/4.d0))
*     
*
          BVDISPBH = -(0.01038961039d0*(2048.d0*XRA**(11.d0/4.d0)*
     &    (1.d0+XRA)**(5.d0/4.d0)-2048.d0*XRA4-2560.d0*XRA3-320.d0*XRA2
     7    +80.d0*XRA-35.d0))/(XRA**(11.d0/4.d0)*(1.d0+XRA)**(5.d0/4.d0))
*
*     ----multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE IF(GAM == 2.0) THEN
          BVDISPBULGE = (0.5000000000d0*(1.d0-12.d0*XRA3-18.d0*XRA2+
     &    12.d0*log(1.d0+XRA)*XRA4-12.d0*log(XRA)*XRA4+24.d0*
     &    log(1.d0+XRA)*XRA3-24.d0*log(XRA)*XRA3+12.d0*log(1.d0+XRA)*
     &    XRA2-12.d0*log(XRA)*XRA2-4.d0*XRA))/
     &    (XRA2*(XRA2+2.d0*XRA+1.d0))
*     
*
          BVDISPBH = -(1.d0/3.d0)*(-1.d0-6.d0*XRA2+12.d0*log(XRA+1.d0)*
     &    XRA4-12.d0*log(XRA)*XRA4+12.d0*log(XRA+1.d0)*XRA3-12.d0*
     &    log(XRA)*XRA3-12.d0*XRA3+2.d0*XRA)/((XRA+1.d0)*XRA3)
*
*     ----multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE
          write(6,*) "f(x,gamma) and h(x,gamma) functions not defined"
          write(6,*) "for gamma =", GAM,". These needed to be added to"
          write(6,*) "dynfvars.f. Please consult the README file."
          stop
      END IF
*
*
      BVDISP = ((GMB/AR)*(XRA**GAM)*(XRA+1.d0)**(4.d0-GAM))*
     &        (BVDISPBULGE+BVDISPBH)
      BVDISP = SQRT(BVDISP)
*
*     ----Calcuate BRATIO for BDIST
      BRATIO = VGCOREMAG/(SQRT(2.0)*BVDISP)
*    
*       Calculate fraction of stars with V*<VCORE for Gamma > 0.0
      BDIST = erf(BRATIO) - (2.0*BRATIO/SQRT(ONEPI))*exp(-(BRATIO**2.0))

*     --for gamma = 0 integrate the true distribution function.
    4 IF (GAM == 0.0) THEN

         !Calculate vesc
         VESC = sqrt(2.d0*(GMB/AR)*(1.d0/(2.d0-GAM))*
     &         ( 1.d0 - (RGDENSMAG/(RGDENSMAG+AR))**(2.d0-GAM)))
         
         !integrate f(v) from 0 to vs (Trapezium rule)
         FV_VS = 0.d0
         do I = 1,999
            FVSTAR = 1d-12 + ((VGCOREMAG-2d-12)*real(I)/1000.0)
            FV_VS = FV_VS + 2.0*fv(FVSTAR)
         end do
         FVSTAR = 1d-12
         FV_VS = FV_VS + fv(FVSTAR)
         VSTAR = VGCOREMAG-1d-12
         FV_VS = FV_VS + fv(FVSTAR)

         FV_VS = FV_VS * ((VGCOREMAG-2d-12)/1000.0)/2.0

         !integrate f(v) from 0 to vesc (Trapezium rule)
         FV_ESC = 0.d0
         do I = 1,999
            FVSTAR = VGCOREMAG+1d-12 + ((VESC-VGCOREMAG-2d-12)*real(I)
     &               /1000.0)
            FV_ESC = FV_ESC + 2.0*fv(FVSTAR)
         end do
         FVSTAR = VGCOREMAG + 1d-12
         FV_ESC = FV_ESC + fv(FVSTAR)
         FVSTAR = VESC-1d-12
         FV_ESC = FV_ESC + fv(FVSTAR)

         FV_ESC = FV_ESC * ((VESC-VGCOREMAG-2.0d-12)/1000.0)/2.0

         !calc BDIST = f(v)_vs/ f(v)_vesc
         BDIST = FV_VS/(FV_VS+FV_ESC)
         write(205,*) RGDENSMAG, BDIST

      END IF
*
*

*
*        Calculate VDOT for dynfri calc
      
      DO 1 K = 1,3
          IF (TTOT >= DTADJ) THEN
              VDOT(K) = (VGCORE(K) - VPREV(K))  /  (TTOT - DYNFTPREV)
          END IF
          VPREV(K) = VGCORE(K)
    1 CONTINUE
      DYNFTPREV = TTOT



*       Jpetts- calculate tidal radius, itterate and update membership

*       Jpetts - Asume Dehnen profile
      RG_RT = RGDENSMAG
      RT_D2PHIDR2 = -AR*((RG_RT/(RG_RT+AR))**(-GAM))
     &            *(AR*GAM - AR + 2.0*RG_RT)
      RT_D2PHIDR2 = RT_D2PHIDR2/(AR+RG_RT)**4
      RT_D2PHIDR2 = RT_D2PHIDR2*GMB/AR
      RT_OMEGA2 = (GMB*(RG_RT/(RG_RT+AR))**(3.0-GAM))/RG_RT**3


      !take 90% (arbitrarily) of old mass as starting point
      RTIDE_TEMP = (0.9*MASSCL/(RT_OMEGA2-RT_D2PHIDR2))**(1.0/3.0)
      MASSCL = 0.0
      NBOUND = 0
      DO 12 I = IFIRST,NTOT
         BOUND(I) = 0
 12   CONTINUE


      DO 2 I = IFIRST,NTOT
         RDIST_RT = SQRT((X(1,I)-RDENS(1))**2 + 
     &      (X(2,I)-RDENS(2))**2 + (X(3,I)-RDENS(3))**2)
         IF (RDIST_RT < RTIDE_TEMP) THEN
            BOUND(I) = 1
            NBOUND = NBOUND + 1
            MASSCL = MASSCL + BODY(I)
         END IF
    2 CONTINUE

      DO 3
        RTIDE_TEMP = (MASSCL/(RT_OMEGA2-RT_D2PHIDR2))**(1.0/3.0)
        RT_CHECK = 0
        DO 31 I = IFIRST,NTOT
            RDIST_RT = SQRT((X(1,I)-RDENS(1))**2 + 
     &      (X(2,I)-RDENS(2))**2 + (X(3,I)-RDENS(3))**2)
            IF (RDIST_RT < RTIDE_TEMP .and. BOUND(I)==0) THEN
                BOUND(I) = 1
                NBOUND = NBOUND + 1
                MASSCL = MASSCL + BODY(I)
                RT_CHECK = 1
            END IF
   31 CONTINUE
      
      IF (RT_CHECK == 0) EXIT

    3 CONTINUE

*     ----Calculate the dynamical friction coeffecient - zero in stalling state


* ----- Calculate apocentre -----
*

      !Check not at apocentre
      if ((((RGDENS(1)*VGCORE(1)) + (RGDENS(2)*VGCORE(2)) + 
     &      (RGDENS(3)*VGCORE(3))) .eq. 0.d0) .AND. 
     &  (VGCOREMAG .le. sqrt(GMB*RGDENSMAG**(2.0-GAM)/
     &   (RGDENSMAG+AR)**(3.0-GAM)))) THEN
      ELSE
*        Jpetts - NOTE 10.0*AR is arbitrary, code will break if apo further than 5.0
*                 unless one increases the limit past the maximum apocentre distance below.
         E_APO = 0.5*VGCOREMAG**2 -(GMB/AR)*(1.d0/(2.d0-GAM))*
     &   (1.d0 - (RGDENSMAG/(RGDENSMAG+AR))**(2.d0-GAM))
         L2_APO = (RGDENS(2)*VGCORE(3) - RGDENS(3)*VGCORE(2))**2 +
     &            (RGDENS(3)*VGCORE(1) - RGDENS(1)*VGCORE(3))**2 +
     &            (RGDENS(1)*VGCORE(2) - RGDENS(2)*VGCORE(1))**2
         RG_APO = dfapo(1d-9,RGDENSMAG-1d-9,10.0*AR,GMB,AR,GAM,
     &                 E_APO,L2_APO)
      END IF

*       Jpetts - Calculate circular orbit tidal radius at apocentre of orbit for stalling

*       Jpetts - Asume Dehnen profile
      RT_D2PHIDR2 = -AR*((RG_APO/(RG_APO+AR))**(-GAM))
     &            *(AR*GAM - AR + 2.0*RG_APO)
      RT_D2PHIDR2 = RT_D2PHIDR2/(AR+RG_APO)**4
      RT_D2PHIDR2 = RT_D2PHIDR2*GMB/AR
      RT_OMEGA2 = (GMB*(RG_APO/(RG_APO+AR))**(3.0-GAM))/RG_APO**3
      RTIDE_APO = (MASSCL/(RT_OMEGA2-RT_D2PHIDR2))**(1.0/3.0)
*
      IF (RG_APO .lt. RTIDE_APO .AND. DFOFF == 0) THEN
         DYNFCOEF = 0.d0
         write(212,*) "RTIDE_APO > RGDENSMAG:", RTIDE_APO, ">"
     &        , RG_APO, "DF SWITCHED OFF"
         DFOFF = 1
      ELSEIF (DFOFF == 1) THEN
         DYNFCOEF = 0.d0
      ELSE
         DYNFCOEF = 4.d0*ONEPI*MASSCL*COULOG*BDENS*BDIST
      END IF
*


*
      RETURN
*
      END
