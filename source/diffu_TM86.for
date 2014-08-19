      subroutine diffu_TM86
      use alloc_1d,only : th,qv,u,v,dift,difk,Fv,dzits,
     :                    ust_s,tst_s,nz,z,dz,def13,def23,
     :                    g,hbl,z_sl,qif,h3,h3c,difunt,difunv,difunu,
     :                    cdm
      implicit none
      integer iz
      real*8 def, ws0,wstar,eps,thv,fi_m,fi_h,Pr0,gammah
      real*8 lmax,mixl,xn,rich,dthdz,dudz
      real, parameter:: karm=0.4
      real, parameter:: b_0=6.5
      real, parameter:: b = 0.
      real, parameter:: zero = 1.e-8
            
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
      
      eps=z_sl/hbl
      wstar=(g/thv*Fv*hbl)**(1./3.)
      ws0 = (ust_s**3. + 7.*eps*karm*wstar**3.)**(1./3.)
      gammah=b_0*Fv/ws0/hbl
      if (Fv.gt.0) then
        fi_m=(1.-7.*dzits)**(-1./3.)
        fi_h=(1.-16.*dzits)**(-1./2.)
      else
        fi_m=1+4.7*dzits
        fi_h=fi_m
      endif  
      Pr0=fi_h/fi_m+b_0*eps*karm
!-------------maximum mixing lenth (for use in Blackadar f-la)------!
       if(ust_s*tst_s.ge.0)then      ! stable stratification
         lmax=40.         
       else                          !  unstable stratification
         lmax=0.15*hbl
       endif
!-------------------------------------------------------------------!
!---------for use in local closure / or Prandtl number--------------!
      do iz=2,nz
        def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
        def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
      enddo
!------------start of the loop for coefficients---------------------!       
      do iz=2,nz-1
        def=dsqrt(((u(iz+1,2)-u(iz-1,2))/(z(iz+1)-z(iz-1)))**2.+
     :         ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.)        
        mixl=karm*z(iz)/(1.+karm*z(iz)/lmax)
        xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
        rich=xn/(def*def+zero)
      
          if (Fv.le.0) then 
	      if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk(iz)=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift(iz)=difk(iz)
	        else
	          difk(iz)=mixl**2*def*b
	          dift(iz)=difk(iz)
	        endif  
	      else
	        difk(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
	      endif
	    else
	    
	    if(z(iz).gt.hbl) then
	
	      if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk(iz)=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift(iz)=difk(iz)
	        else
	          difk(iz)=mixl**2*def*b
	          dift(iz)=difk(iz)
	        endif  
	      else
	        difk(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
	      endif
	    
	    else      
	   
            difk(iz)=karm*ws0*z(iz)*(1.-z(iz)/hbl)**2.
            dift(iz)=difk(iz)/Pr0
            
          endif
        endif
      enddo
!-------------at the lowest level ------------------------!         
      if(-tst_s.le.0) then
      dudz=(1.+5.*dzits)*ust_s/(karm*z_sl)
      dthdz=(1.+5.*dzits)*tst_s/(karm*z_sl)
      endif
	if(-tst_s.gt.0) then
	dudz=(1.-16.*dzits)**(0.25)*ust_s/(karm*z_sl)
      dthdz=(1.-16.*dzits)**(0.5)*tst_s/(karm*z_sl)
      endif
      mixl=karm*z_sl/(1.+karm*z_sl/lmax)
      def=dudz
      xn=g*dthdz/th(1,2)
      rich=xn/(def*def+zero)
      if(-tst_s.le.0) then
          difk(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
	    dift(1)=difk(1)
	else
	  difk(1)=karm*ws0*z_sl*(1.-z_sl/hbl)**2.
        dift(1)=difk(1)/Pr0
	endif
	
	do iz=2,nz
        def13(iz)=def13(iz)*0.5*(difk(iz)+difk(iz-1))
        def23(iz)=def23(iz)*0.5*(difk(iz)+difk(iz-1))
      enddo
      
      def13(1)=cdm*u(1,2)
      def23(1)=cdm*v(1,2)
      
      
      do iz=1,nz-1
        difunu(iz)=(def13(iz+1)-def13(iz))/(0.5*(dz(iz+1)+dz(iz)))
        difunv(iz)=(def23(iz+1)-def23(iz))/(0.5*(dz(iz+1)+dz(iz)))
      enddo
      
       do iz=2,nz
         h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
         h3(iz)=h3(iz)*0.5*(dift(iz)+dift(iz-1))
         if(z(iz).le.hbl.and.-tst_s.gt.0) then
	     h3c(iz)=-gammah*0.5*(dift(iz)+dift(iz-1))
	   else
           h3c(iz)=0.
	   endif     
       enddo
      
      h3(1)=ust_s*tst_s
      h3c(1)=0.
      
      do iz=1,nz-1
        difunt(iz)=(h3(iz+1)+h3c(iz+1)-h3(iz)-h3c(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
      enddo
      
      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
	
	
      end