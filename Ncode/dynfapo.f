      function dfapo(tol,rlo,rhi,GMB,AR,GAM,E,L2)
         real*8 tol,rlo,rhi,GMB,AR,GAM,E,L2
         real*8 y0,rm,yy,dfapo
         integer*8 steps

         !TEMP VARS
         !real*8 ri
         !integer*8 i
         !TEMP VARS
         steps=0
         !write(305,*) "Start:", tol, rlo, rhi, GMB, AR, GAM, E, L2

 100     y0=(1/rlo)**2 + 2.0*((-(GMB/AR)*(1.0/(2.0-gam)) * 
     &          (1.0- (rlo/(rlo+AR))**(2.0-gam)))    - E)/L2
         rm=(rlo+rhi)/2.d0
         yy=(1/rm)**2 + 2.0*((-(GMB/AR)*(1.0/(2.0-gam)) * 
     &          (1.0- (rm/(rm+AR))**(2.0-gam)))    - E)/L2
         steps=steps+1
         !write(305,*) steps
         !write(305,*) rlo, rhi
         !write(305,*) y0, yy, y0*yy
         if (abs(yy*y0).eq.0) goto 101
         if ((yy*y0)<0.d0) then
            rhi=rm
         else if ((yy*y0)>0.d0) then
            rlo=rm
         end if
         !write(305,*) rhi, rlo, rhi-rlo
         if (abs(rhi-rlo)>tol) goto 100
         !write(305,*)
 101     dfapo = rm
      end function
