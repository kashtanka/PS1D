      subroutine eq_state
      use alloc_1d
      implicit none
      integer iz,it
      real*8 beta,dthdy,t0
      real*8 C1,C2,F,Fwt,pm
      
      ! calculation of pressure from hydrostatic equation (Newton iteration method 
      ! - can be avoided ) 
      do iz=nz-1,1,-1
      
      C1=g*(z(iz+1)-z(iz))*0.5*p00**akapa/r/(0.5*(th(iz+1,2)+th(iz,2))) 
      C2=p(iz+1,2)
      pm=p(iz+1,2)-1000.                        !first guess
      F=pm-C1*pm**(1.-akapa)-C2                 ! function value with first guess
      Fwt=1.-C1*(1.-akapa)*pm**(-akapa)         ! derivative value with first guess
      !-----iteration------------!
      do it=1,4
      pm=pm-F/Fwt
      F=pm-C1*pm**(1.-akapa)-C2
      Fwt=1.-C1*(1.-akapa)*pm**(-akapa)
      enddo
      p(iz,2)=2.*pm-p(iz+1,2)
      !write(0,*)z(iz), p(iz,2)
      enddo
      
      !---------------------------------------------------!
      
      do iz=1,nz
        t(iz)=th(iz,1)*(p(iz,1)/p00)**akapa
        if (qif.ne.0) then
           ro(iz)=p(iz,2)/(r*t(iz)) ! ro(iz)=p(iz,2)/(r*t(iz)*(1.+0.61*qv(iz,2)))          
        else
          ro(iz)=p(iz,2)/(r*t(iz))
        endif	

    
      enddo
  
      
      end
      
