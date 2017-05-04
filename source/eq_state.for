      subroutine eq_state
      use alloc_1d
      implicit none
      integer iz,it
      real*8 beta,dthdy,t0
      real*8 C1,C2,F,Fwt,pm
      real,external:: qsati
      real*8 z_hl(nz+1)
      
      ! calculation of pressure from hydrostatic equation (Newton iteration method 
      ! - must be avoided in the future 
!      do iz=nz-1,1,-1
      
!      C1=g*(z(iz+1)-z(iz))*0.5*p00**akapa/r/(0.5*(th(iz+1,2)+th(iz,2))) 
!      C2=p(iz+1,2)
!      pm=p(iz+1,2)-1000.                        !first guess
!      F=pm-C1*pm**(1.-akapa)-C2                 ! function value with first guess
!      Fwt=1.-C1*(1.-akapa)*pm**(-akapa)         ! derivative value with first guess
      !-----iteration------------!
!      do it=1,4
!      pm=pm-F/Fwt
!      F=pm-C1*pm**(1.-akapa)-C2
!      Fwt=1.-C1*(1.-akapa)*pm**(-akapa)
!      enddo
!      p(iz,2)=2.*pm-p(iz+1,2)
!      enddo
      
      !---------------------------------------------------!
cccc-pressure at half levels used in the Goddard radiation scheme-cccc
      z_hl(1) = 0.
      z_hl(nz+1) = z(nz)+0.5*(z(nz) -z(nz-1))
      do iz = 2,nz
         z_hl(iz) = 0.5*(z(iz) + z(iz-1)) ! height at half levels
      enddo
      phl(1) = pa               ! pressure right at the surface
      do iz = 2,nz+1
         phl(iz) = (phl(iz-1)**akapa - p00**akapa*g/cp/th(iz-1,2)
     :             *(z_hl(iz)-z_hl(iz-1)))**(1./akapa)
      enddo
cccc- pressure at full levels- cccc
      do iz = 1,nz
         p(iz,2) = 0.5*(phl(iz) + phl(iz+1))
         t(iz)=th(iz,2)*(p(iz,2)/p00)**akapa
!          qv(iz,1)=0.8*qsati(t(iz),p(iz,1))
!         qv(iz,2)=0.8*qsati(t(iz),p(iz,1))
!          qv(iz,3)=0.8*qsati(t(iz),p(iz,1))
         !write(0,*) z(iz),qv(iz,2),0.8*qsati(t(iz),p(iz,1))
         
         if (qif.ne.0) then
            ro(iz)=p(iz,2)/(r*t(iz)) ! ro(iz)=p(iz,2)/(r*t(iz)*(1.+0.61*qv(iz,2)))          
         else
            ro(iz)=p(iz,2)/(r*t(iz))
         endif	
      enddo
     
      end
      
