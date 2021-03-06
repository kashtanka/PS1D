      subroutine diffu_Lock
      use alloc_1d
      implicit none
      integer iz
      real*8 karm,def,lmax,mixl,rich,zero,xn,dudz,dthdz,b,thv,zb,
     :wstar,wm,w_2s,w_2,tstar_f,gammah,gammah_s,zf,Vsc,Vrad,Vbr
      real,parameter:: b_hm=3.
      karm=0.4
      zero=1.e-8
      !b=0.005
       b=0.
       if(qif.gt.0) then
       thv=th(1,2)+0.61*th(1,2)*qv(1,2)
       else
       thv=th(1,2)
       endif
       ! scaling stuff: w*, <w2>_s, th*, and countergradient term at the surface
       wstar=(g/thv*Fv*hbl)**(1./3.)    !what are the names of all these??
	 wm=(ust_s**3.+7.*z_sl/hbl*0.4
     :	*wstar**3.)**(1./3.)
       w_2s=(1.6*ust_s**2.*(1-z_sl/hbl))**(3./2.)+
     :       1.2*wstar**3.*z_sl/hbl*
     :      (1-0.9*z_sl/hbl)**(3./2.)
	 w_2s=w_2s**(2./3.)
	 tstar_f=-ust_s*tst_s/wstar	
	 gammah_s=b_hm*wstar**2*tstar_f/w_2s/hbl
         Vsc=(Vrad**3.+Vbr**3.)**(1./3.)
     
      do iz=2,nz
        def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
        def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
      enddo
       if(ust_s*tst_s.ge.0)then      ! stable stratification
         lmax=40.
         !lmax=0.15*hbl
         
       else                          !  unstable stratification
         !lmax=45.
         lmax=0.15*hbl
       endif
       do iz=2,nz-1
        def=dsqrt(((u(iz+1,2)-u(iz-1,2))/(z(iz+1)-z(iz-1)))**2.+
     :         ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.)
        
        mixl=karm*z(iz)/(1.+karm*z(iz)/lmax)
        xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
        rich=xn/(def*def+zero)
        
        if (Fv.le.0) then !(Fv(ix,iy).le.0) then
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
	
           if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk2(iz)=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift2(iz)=difk2(iz)
	        else
	          difk2(iz)=mixl**2*def*b
	          dift2(iz)=difk2(iz)
	        endif  
	      else
	        difk2(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift2(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
	      endif
     
!      w_2s=1.44*ust_s(ix,iy)**2.*(1.-1.5*dzits(ix,iy))**(2./3.)

!	w_2=(1.6*ust_s**2.*(1-z(iz)/hbl))**(3./2.)+
!     :1.2*wstar**3.*z(iz)/hbl*(1-0.9*z(iz)/hbl)**(3./2.)
!	w_2=w_2**(2./3.)	

	dift(iz)=karm*ust_s*z_sl/((1-16.*dzits)**(-0.5)
     :        -karm*z_sl/tst_s*gammah_s)*
     :        ((hbl-z(iz))/(hbl-z_sl))**2.*
     :        (ust_s*karm*z(iz)+wstar*hbl*
     :        (z(iz)/hbl)**(4./3.))/(ust_s*karm*z_sl+
     :        wstar*hbl*(z_sl/hbl)**(4./3.))
     
	difk(iz)=dift(iz)*((1-16.*dzits)**(-0.25)+b_hm*wstar*
     :        ust_s*karm*z_sl/hbl
     :        /((1-16.*dzits)**(-0.25))/w_2s)    
!     top-down diffusion coefs
        if (ifwr.ne.0) then
           if(z(iz).ge.zb.and.z(iz).le.hbl) then
           dift3(iz)=0.85*karm*Vsc*(z(iz)-zb)**2./(hbl-zb)
     :               *(1-(z(iz)-zb)/(hbl-zb))**0.5
           difk3(iz)=0.75*dift3(iz)
           endif
        endif

	dift(iz)=max(dift(iz),dift2(iz),dift3(iz))
	difk(iz)=max(difk(iz),difk2(iz),dift3(iz))

	    endif
	    endif	
	enddo    
	
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
	endif

	if(-tst_s.gt.0) then
	dift(1)=karm*ust_s*z_sl/((1-16.*dzits)**(-0.5)
     :        -karm*z_sl/tst_s*gammah_s)*
     :        ((hbl-z_sl)/(hbl-z_sl))**2.*
     :        (ust_s*karm*z_sl+wstar*hbl*
     :        (z_sl/hbl)**(4./3.))/(ust_s*karm*z_sl+
     :        wstar*hbl*(z_sl/hbl)**(4./3.))
     
	difk(1)=dift(1)*((1-16.*dzits)**(-0.25)+b_hm*wstar*
     :        ust_s*karm*z_sl/hbl
     :        /((1-16.*dzits)**(-0.25))/w_2s)
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
      
      
       do iz=2,nz-1
         h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
         h3(iz)=h3(iz)*0.5*(dift(iz)+dift(iz-1))
         if(z(iz-1).le.hbl.and.-tst_s.gt.0) then
           zf=z(iz-1)
	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
	     w_2=w_2**(2./3.)
	     gammah=b_hm*wstar**2*tstar_f/w_2/hbl
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
      
