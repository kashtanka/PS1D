      subroutine baroclinity
      use alloc_1d
      implicit none
      integer iz,it
      real*8 beta,dthdy,t0
      real*8 C1,C2,F,Fwt,pm
      
 !     p(nz,2)=ptop      !known value
 !     do iz=nz-1,1,-1
      
 !     C1=g*(z(iz+1)-z(iz))*0.5*p00**akapa/r/(0.5*(th(iz+1,2)+th(iz,2))) 
 !     C2=p(iz+1,2)
 !     pm=p(iz+1,2)-1000.                                 !first guess
 !     F=pm-C1*pm**(1.-akapa)-C2                                 ! function value with first guess
 !     Fwt=1.-C1*(1.-akapa)*pm**(-akapa)                         ! derivative value with first guess
      !-----iteration------------!
 !     do it=1,4
 !     pm=pm-F/Fwt
 !     F=pm-C1*pm**(1.-akapa)-C2
 !     Fwt=1.-C1*(1.-akapa)*pm**(-akapa)
 !     enddo
 !     p(iz,2)=2.*pm-p(iz+1,2)
      !write(0,*)z(iz), p(iz,2)
 !     enddo
      
 !     p0=p(1,2) 
        
 !     t=th(1,2)*(p0/p00)**akapa	
 !     ro=p0/(r*t)     
 !     ug_bar=ug
 !     t0=0.5*(th(1,2)+th(1,1))*(p0/p00)**akapa
 !     if (dt*nstep.gt.3600.*icetime) then
 !     do iz=nz,1,-1
 !     if (z(iz).lt.hbl) then
 !     beta=g/t0
 !     dthdy=0.5*(th(iz,2)+th(iz+1,2)-th(iz,1)-th(iz+1,1))/dy
 !     ug_bar(iz)=ug_bar(iz+1)+beta/fcor*dthdy*(z(iz+1)-z(iz))
 !     endif
 !     enddo
 !     dpdy_d=ro*fcor*(ug-ug_bar(1))
      !p0=p0+dpdy_d*dy
 !     endif
      
      end