      function fv(VS)
        include "common6.h"
        COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
        real*8 :: VS, E_DIM, fv
*
*       Jpetts - calculate fv(vs). (note this is only for a Dehnen model
*                with gamma = 0.
        !Calculate dimensionless energy for star at v and r.
        E_DIM = 0.5*VS**2 -(GMB/AR)*(1.d0/(2.d0-GAM))*
     & ( 1.d0 - (RGDENSMAG/(RGDENSMAG+AR))**(2.d0-GAM))
        E_DIM = -E_DIM*(AR/GMB)
*
        !Calculate f(E) (gamma = 0) (Dehnen 1993)
        fv = ((3.d0*GMB/(2.d0*(ONEPI**3)*(GMB*AR)**1.5))*
     &  ((sqrt(2.0*E_DIM) * (3.0-4.0*E_DIM)/(1.0-2.0*E_DIM))
     &   - 3.0*asinh(sqrt(2.0*E_DIM/(1.0-2.0*E_DIM)))))
*
        fv = fv*4.0*ONEPI*VS**2
*
       end function
