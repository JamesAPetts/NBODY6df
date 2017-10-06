      SUBROUTINE DYNFVARS
*
*       Jpetts - Calculate DYNFCOEF et al
*       Recalculate variables used for dynamical friction formula every
*       time parameters are adjusted
*       ----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8 XRA,XRA2,XRA3,XRA4,MASSENC
*
*
*        Calculate r
      DO 10 K = 1,3
          RGDENS(K) = RG(K) + RDENS(K)
   10 CONTINUE   
      RGDENSMAG = SQRT(RGDENS(1)**2.0+RGDENS(2)**2.0+RGDENS(3)**2.0)
*
*       Calculate the circular velocity at RG.

      MASSENC = GMG + GMB*(RGDENSMAG/(RGDENSMAG+AR))**(3.0-GAM)
      VCIRC2 = MASSENC/RGDENSMAG

*       Calculate the coulomb logarithm.
      !make a selection for bmax based on gamma
      COULOG = min(RGDENSMAG/GAM,RGDENSMAG) !b_max
      !COULOG = RGDENSMAG*(RGDENSMAG+AR)/(AR*GAM + 4.d0*RGDENSMAG)
      COULOG = COULOG/max(1.4d0*RSCALE/1.3d0,MASSCL/VCIRC2) !/b_min
      COULOG = log(max(COULOG,1.0)) !(if b_min >= b_max, set to 0.0)
*
*       Calculate the local background density.
      BDENS = (3.0*GMB*AR)/(8.0*ONEPI)
      BDENS = BDENS/(RGDENSMAG**GAM)
      BDENS = BDENS/((RGDENSMAG+AR)**(4.0-GAM))
*
********Calculate the local velocity dispersion*************************
*
*  Included here is the general velocity disperion profile for a Dehnen
*  model with (or without) a central point mass black hole. The 
*  derivation is general in that it allows one to use any SMBH/
*  background mass ratio, any slope gamma and any scale radius. The 
*  functions f(x,gam) and h(x,gam) must be pre-calculated and included 
*  in the code. gam = 0.0,0.5,1.0, 1.5 and 2.0 are included by default. 
*  Adding another value of gamma is straightforward and described in the
*  README.
*
*  If one wants a simulation with no SMBH, simply set GMG to 0.0 in the
*  input file.
*
*  One needs to derive the velocity disperion profile of the
*  background distribution themselves if they want to use a density
*  distribution other than a Dehnen model. The functions f(x,gam) and h(x,gam)
*  are only analytic
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
*     select apropriate calculated pre-calculated integrals
      IF(GAM == 0.0) THEN
          BVDISPBULGE = (1.d0/30.d0)*(1.d0+6.d0*XRA)/(XRA+1.d0)**6.d0
*     
*
          BVDISPBH = -(1.d0/3.d0)*(-3.d0-22.d0*XRA+12.d0*log(XRA+1.d0)*
     &    XRA4-12.d0*log(XRA)*XRA4+36.d0*log(XRA+1.d0)*XRA3-36.d0*
     &    log(XRA)*XRA3+36.d0*log(XRA+1.d0)*XRA2-36.d0*log(XRA)*XRA2+
     &    12.d0*log(XRA+1.d0)*XRA-12.d0*log(XRA)*XRA-30.d0*XRA2-12.d0
     &    *XRA3)/(XRA*(XRA3+3.d0*XRA2+3.d0*XRA+1.d0))
*
*     multiply by mu
          BVDISPBH = BVDISPBH*(GMG/GMB)
      ELSE IF(GAM == 0.5) THEN
          BVDISPBULGE = 0.2d0/(1.d0+XRA)**5.d0
*     
*
          BVDISPBH = (0.1333333333d0*(128.d0*XRA**(3.d0/2.d0)*(1.d0+XRA)
     &    **(5.d0/2.d0)-128.d0*XRA4-320.d0*XRA3-240.d0*XRA2-40.d0*XRA
     &    +5.d0))/(XRA**(3.d0/2.d0)*(1.d0+XRA)**(5.d0/2.d0))
*
*     multiply by mu
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
*     multiply by mu
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
*     multiply by mu
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
*     multiply by mu
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
*     multiply by mu
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
*       Calculate VGCORE & VGCOREMAG
      DO 20 K = 1,3
          VGCORE(K) = VG(K) + VCORE(K)
   20 CONTINUE 
      VGCOREMAG = SQRT(VGCORE(1)**2.0+VGCORE(2)**2.0+VGCORE(3)**2.0)

*       Calcuate BRATIO for BDIST
      BRATIO = VGCOREMAG/(SQRT(2.0)*BVDISP)
*    
*       Calculate fraction of stars with V*<VCORE
      BDIST = erf(BRATIO) - (2.0*BRATIO/SQRT(ONEPI))*exp(-(BRATIO**2.0))
*
*
*       Calculate the dynamical friction coeffecient
*
      DYNFCOEF = 4.0*ONEPI*MASSCL*COULOG*BDENS*BDIST
*
*        Calculate VDOT for dynfri calc
      
      DO 42 K = 1,3
          IF (TTOT >= DTADJ) THEN
              VDOT(K) = (VGCORE(K) - VPREV(K))  /  (TTOT - DYNFTPREV)
          END IF
          VPREV(K) = VGCORE(K)
   42 CONTINUE
      DYNFTPREV = TTOT

*
      RETURN
*
      END
