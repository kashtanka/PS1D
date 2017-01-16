      subroutine diffu_LS96
      use alloc_1d
      implicit none
      integer iz
      real*8 karm,def,lmax,mixl,rich,zero,xn,dudz,dthdz,b,thv,
     :wstar,wm,w_2s,w_2,tstar_f,gammah_s,zf,surf_flux
      real,parameter:: b_hm=3.
      real, external:: qsat
      real*8 gammaq(1:nz),gammah(1:nz),gammaq_hl(1:nz),gammah_hl(1:nz)
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     : difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz)
      real*8 difk_st,f_ri
      dift_hl=0.
      difk_hl=0.
      gammah_hl=0.
      gammaq_hl=0.
      gammah=0.
      gammaq=0.
      karm=0.4
      zero=1.e-9
      b=0.
      if(qif.gt.0) then
         thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
         thv=th(1,2)
      endif
!----------scaling parameters--------------------!
c-------zatychka na vremya
      hbl = max(hbl,10.)
c-----------------------      
      wstar=(g/thv*abs(Fv)*hbl)**(1./3.)
      wm=(ust_s**3.+7.*z_sl/hbl*0.4
     :	*wstar**3.)**(1./3.)
      w_2s=(1.6*ust_s**2.*(1-z_sl/hbl))**(3./2.)+
     :       1.2*wstar**3.*z_sl/hbl*
     :      (1-0.9*z_sl/hbl)**(3./2.)
      w_2s=w_2s**(2./3.)
      tstar_f=-(ust_s*tst_s+0.61*thv*ust_s*qst_s)/wstar	
      gammah_s=b_hm*wstar**2*tstar_f/w_2s/hbl
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
         if(qif.ne.0) then 
           if(ifwr.eq.0.) then              ! virtual theta
           xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2))
     :      -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :      /(z(iz+1)-z(iz-1))/
     :     th(iz,2)/(1.+0.61*qv(iz,2))
           else   
            if(qv(iz,2)
     :             .ge.0.95*qsat(t(iz),p(iz,2))) then        
              xn=g*((th(iz+1,2)-hlat/cp*qc(iz+1,2))
     :        -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
     :        /(z(iz+1)-z(iz-1))/th(iz,2)
             else
            xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2))
     :      -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :      /(z(iz+1)-z(iz-1))/
     :     th(iz,2)/(1.+0.61*qv(iz,2))
             endif
           endif
           else
           xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
         endif

!        xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
        rich=xn/(def*def+zero)
        f_ri = (1.+5.*rich+44.*rich**2.)**(-2.) !(max(0. , 1.-5.*rich))**2.
        difk_st =  mixl**2.*def*f_ri !max(b,f_ri)
!        f_ri = (max(0. , 1.-5.*rich))**2.
!        difk_st = mixl**2.*def*max(b,f_ri)
        if (Fv.le.0) then 
           if(rich.ge.0) then
              difk(iz) = difk_st
              dift(iz)=difk(iz)
           else
              difk(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
              dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
           endif
        else
           if(z(iz).gt.hbl) then
              if(rich.ge.0) then
                 difk(iz)=difk_st
                 dift(iz)=difk(iz)
              else
                 difk(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
                 dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
              endif
           else      
              if(rich.ge.0) then
                 difk2(iz) = difk_st
                 dift2(iz) = difk2(iz)
              else
                 difk2(iz)=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
                 dift2(iz)=difk2(iz)*(min(3.,(1.-16.*rich)**0.25))
              endif
               
              dift(iz)=karm*ust_s*z_sl/((1-16.*dzits)**(-0.5)
     :             -karm*z_sl/tst_s*gammah_s)*
     :             ((hbl-z(iz))/(hbl-z_sl))**2.*
     :             (ust_s*karm*z(iz)+wstar*hbl*
     :             (z(iz)/hbl)**(4./3.))/(ust_s*karm*z_sl+
     :             wstar*hbl*(z_sl/hbl)**(4./3.))
               
              difk(iz)=dift(iz)*((1-16.*dzits)**(-0.25)+b_hm*wstar*
     :             ust_s*karm*z_sl/hbl
     :             /((1-16.*dzits)**(-0.25))/w_2s)  
         
              dift3(iz)=0.85*karm*(min(1.,dR*hbl))**(1./3.)*(z(iz)-0)**2.
     :             /(hbl-0)*(1.-(z(iz)-0)/(hbl-0))**0.5
              difk3(iz)=0.75*dift3(iz)
          
              dift(iz) = max(dift(iz),dift2(iz))
              difk(iz) = max(difk(iz),difk2(iz))
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
         if(z(iz).lt.hbl.and.-(tst_s+0.61*thv*qst_s).gt.0) then
           zf=z(iz)
	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
	     w_2=w_2**(2./3.)
	     gammah(iz)=b_hm*wstar**2*tstar_f/w_2/hbl
             gammaq(iz)=b_hm*wstar*(-ust_s*qst_s)/w_2/hbl
          else
             gammah(iz)=0.
             gammaq(iz)=0.
          endif
       enddo

      do iz = 2,nz
         dift_hl(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl(iz)=0.5*(difk(iz-1)+difk(iz))
         gammah_hl(iz)=0.5*(gammah(iz-1)+gammah(iz))
         gammaq_hl(iz)=0.5*(gammaq(iz-1)+gammaq(iz))
      enddo

!      if (implicit) then
      do iz=1,nz-1
      F(iz)=th(iz,1)+dtl/(0.5*dz(iz)+dz(iz+1))*
     :        (gammah_hl(iz)*dift_hl(iz)
     :         -gammah_hl(iz+1)*dift_hl(iz+1))         ! right-hand side
      enddo
      F(nz)=th(nz,1)+dtl/dz(nz)*gammah_hl(nz)*dift_hl(nz)
      if (seaice.eq.1.and.frac.lt.1) then
         surf_flux = frac*ust_s*tst_s + (1.-frac)*ust_s2*tst_s2
      else
         surf_flux = ust_s*tst_s
      endif
      call implicit_dif(F,th(1:nz,1),dift_hl,surf_flux,difunt2)
      F(1:nz)=u(1:nz,1)
      if (seaice.eq.1.and.frac.lt.1) then
         surf_flux = u(1,2)*(frac*cdm + (1.-frac)*cdm2)
      else
         surf_flux = cdm*u(1,2)
      endif
      call implicit_dif(F,u(1:nz,1),difk_hl,surf_flux,difunu2)
      F(1:nz)=v(1:nz,1)
      if (seaice.eq.1.and.frac.lt.1) then
         surf_flux = v(1,2)*(frac*cdm + (1.-frac)*cdm2)
      else
         surf_flux = cdm*v(1,2)
      endif
      call implicit_dif(F,v(1:nz,1),difk_hl,surf_flux,difunv2)
      do iz=1,nz-1
      F(iz)=qv(iz,1)+dtl/(0.5*dz(iz)+dz(iz+1))*
     :        (gammaq_hl(iz)*dift_hl(iz)
     :         -gammaq_hl(iz+1)*dift_hl(iz+1))         ! right-hand side
      enddo
      F(nz)=qv(nz,1)+dtl/dz(nz)*gammaq_hl(nz)*dift_hl(nz)
      if (seaice.eq.1.and.frac.lt.1) then
         surf_flux = frac*ust_s*qst_s + (1.-frac)*ust_s2*qst_s2
      else
         surf_flux = ust_s*qst_s
      endif
      call implicit_dif(F,qv(1:nz,1),dift_hl,surf_flux,difunqv2)
      F(1:nz)=qc(1:nz,1)
      call implicit_dif(F,qc(1:nz,1),dift_hl,0.,difunqc2)
      F(1:nz)=qr(1:nz,1)
      call implicit_dif(F,qr(1:nz,1),dift_hl,0.,difunqr2)
!         else
!         endif
	
      do iz=2,nz
        def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
        def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
        def13(iz)=def13(iz)*difk_hl(iz)
        def23(iz)=def23(iz)*difk_hl(iz)
      enddo
!--------------surface momentum flux-------------!
      if (frac.lt.1) then
         def13(1) = u(1,2)*(frac*cdm + (1.-frac)*cdm2)
         def23(1) = v(1,2)*(frac*cdm + (1.-frac)*cdm2)
      else
         def13(1) = cdm*u(1,2)
         def23(1) = cdm*v(1,2)
      endif
      
      do iz=1,nz-1
        difunu(iz)=(def13(iz+1)-def13(iz))/(0.5*(dz(iz+1)+dz(iz)))
        difunv(iz)=(def23(iz+1)-def23(iz))/(0.5*(dz(iz+1)+dz(iz)))
      enddo
      
      do iz=2,nz-1
         h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
         h3(iz)=h3(iz)*dift_hl(iz)
	 h3c(iz)=-gammah_hl(iz)*dift_hl(iz)
       enddo
      
!-----------surface heat flux ----------------!
      if (frac.lt.1) then
         h3(1) = frac*ust_s*tst_s + (1.-frac)*ust_s2*tst_s2
      else
         h3(1) = ust_s*tst_s
      endif
      h3c(1)=0.

       do iz=2,nz-1
         wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
         wq3(iz)=wq3(iz)*dift_hl(iz)
	 wq3_c(iz)=-gammaq_hl(iz)*dift_hl(iz)    
       enddo
 !-----------surface moisture flux-------------------------!     
       if (frac.lt.1) then
         wq3(1) = frac*ust_s*qst_s + (1.-frac)*ust_s2*qst_s2
      else
         wq3(1) = ust_s*qst_s
      endif
      wq3_c(1)=0.

       if (ifwr.ne.0) then
        do iz=2,nz
          wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
          wq3c(iz)=wq3c(iz)*dift_hl(iz)
          wq3r(iz)=(qr(iz,2)-qr(iz-1,2))/dz(iz)
          wq3r(iz)=wq3r(iz)*dift_hl(iz)
        enddo
      endif
 
      do iz=1,nz-1
        difunt(iz)=(h3(iz+1)+h3c(iz+1)-h3(iz)-h3c(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
        difunqv(iz)=(wq3(iz+1)+wq3_c(iz+1)-wq3(iz)-wq3_c(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
      enddo

       if (ifwr.ne.0) then
        do iz=1,nz-1
          difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
          difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
        enddo
      endif
      
!      do iz=1,nz
!      write(0,*) iz,difunt2(iz),difunt(iz)iu
!      enddo

      difunu(1:nz)=difunu2
      difunt(1:nz)=difunt2
      difunv(1:nz)=difunv2
      difunqv(1:nz)=difunqv2
      difunqc(1:nz)=difunqc2
      difunqr(1:nz)=difunqr2
      
      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
      difunqv(nz)=0.
      difunqc(nz)=0.
      difunqr(nz)=0.
      end
