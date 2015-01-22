      subroutine diffu_Noh03
      use alloc_1d,only : th,qv,u,v,dift,difk,Fv,dzits,
     :                    ust_s,tst_s,nz,z,dz,def13,def23,
     :                    g,hbl,z_sl,qif,h3,h3c,difunt,difunv,difunu,
     :                    cdm,h3e,def13c,def23c,difunqv,difunqc,difunqr,
     :                    qst_s,wq3,wq3c,wq3r,wq3_c,qc,qr,ifwr,cp,hlat,
     :                    t,p,dift2,difk2,dR,wq3e,dift3,difk3,ro,dtl,rad
     :                    ,vat,wth_h,wth_h2,we,condensat,hlat,akapa,p00,
     :                    the,tende
      implicit none
      integer iz
      real*8 def, ws0,wstar,eps,thv,fi_m,fi_h,Pr0,gammah
      real*8 lmax,mixl,xn,rich,dthdz,dudz,w_m3,ws,Pr
      real*8 delta,w_2,tstar_f,zf,gammah2,wsm,gamma_u,gamma_v
      real*8 gammaq,wq_h,eps2,dthv
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 dift_hl2(1:nz),difk_hl2(1:nz),dift_hl3(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     : difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz)
      real, external :: qsat
      real, parameter:: karm=0.4
      real, parameter:: b_0=6.5
      real, parameter:: b_hm=3.
      real*8, parameter:: b = 0.
      real, parameter:: zero = 1.e-4
      real, parameter:: B2 = 5.
      real, parameter:: A = 4.5
      real, parameter:: alpha = 3.
      real, parameter:: Sm=15.9
      h3e=0.
      dift_hl=0.
      difk_hl=0.
      dift_hl2=0.
      difk_hl2=0.
      dift_hl3=0.
      difunt2=0.
      difk=0.
      dift=0.
      difk2=0.
      dift2=0.
      dift3=0.
      difk3=0.
!--------virtual potential temperature---------------------            
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
!----surface layer fraction      
      eps=z_sl/hbl
!--------scaling parameters----------------------------------
      wstar=(g/thv*Fv*hbl)**(1./3.) 
      ws0 = (ust_s**3. + 7.*eps*karm*wstar**3.)**(1./3.)
      wsm=(ust_s**3. + 7.*karm*wstar**3.*0.5)**(1./3.)
      tstar_f=-ust_s*tst_s/wstar
      gammah=b_0*Fv/wsm/hbl
      gammaq=-b_0*qst_s*ust_s/wsm/hbl
!--------for the nonlocal transport of momentum--------------
!      gamma_u=-Sm*ust_s**2./wsm/hbl*(wstar/wsm)**3.
!     :        *u(1,2)/sqrt(u(1,2)**2.+v(1,2)**2.)
!      gamma_v=-Sm*ust_s**2./wsm/hbl*(wstar/wsm)**3.
!     :        *v(1,2)/sqrt(u(1,2)**2.+v(1,2)**2.)
!-------entrainment------------------------------------------
      w_m3=(wstar**3.+B2*ust_s**3.) +hbl*dR*g/thv
      wth_h= -6.*w_m3/hbl    
!---------surface Prandtl number-----------------!
      if (Fv.gt.0) then
        fi_m=(1.-7.*dzits)**(-1./3.)
        fi_h=(1.-16.*dzits)**(-1./2.)
      else
        fi_m=1+4.7*dzits
        fi_h=fi_m
      endif  
      Pr0=fi_h/fi_m+b_0*eps*karm
!------------------------------------------------!
!---------------matching K-profile with the entrainment flux------!
      do iz = 2,nz-1
      if (z(iz).gt.hbl.and.z(iz-1).le.hbl) then
         dthv=th(iz,2)+th(iz,2)*(0.61*qv(iz,2)-qc(iz,2)
     :        )-th(iz-1,2)-th(iz-1,2)*(0.61*qv(iz-1,2)-qc(iz-1,2))
         we= wth_h/dthv*100.
          Pr=1.+(Pr0-1.)*exp(-alpha*(z(iz-1)-eps*hbl)**2./hbl**2.)
         ws=(ust_s**3.+7.*karm*wstar**3.*z(iz-1)/hbl)**(1./3.) 
         eps2=hbl/z(iz-1)*
     :   (1.-sqrt(max(0.,-2.*Pr*(wth_h)*dz(iz)/
     :   (karm*ws*hbl*((dthv-gammah*dz(iz)))))))       
         eps2=min(1.,max(0.7,eps2))         
         endif
      enddo
!-------------maximum mixing length (for use in Blackadar f-la)------!
       if(Fv.lt.0)then      ! stable stratification
         lmax=40.         
       else                          !  unstable stratification
         lmax=300. !0.15*hbl
       endif
!-------------------------------------------------------------------!
!---------for use in local closure / or Prandtl number--------------!
      do iz=2,nz
        def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
        def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
      enddo
!------------start of the loop for coefficients---------------------!       
      do iz=2,nz-1
         def=max(zero,dsqrt(((u(iz+1,2)-u(iz-1,2))
     :           /(z(iz+1)-z(iz-1)))**2.+
     :         ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.))       
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
            xn=g*(th(iz+1,2)*(1.+0.61*qv(iz+1,2))
     :      -th(iz-1,2)*(1.+0.61*qv(iz-1,2)))
     :      /(z(iz+1)-z(iz-1))/
     :     th(iz,2)/(1.+0.61*qv(iz,2))
             endif
           endif
           else
           xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
         endif

        rich=xn/(def*def+zero)
	 if(rich.ge.0) then
            difk2(iz)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
	    dift2(iz)=difk2(iz)
	 else
	    difk2(iz)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift2(iz)=difk2(iz)*(min(3.,(1.-16.*rich)**0.25))
         endif
         
         if(Fv.gt.0.and.z(iz).le.hbl) then
	    ws=(ust_s**3.+7.*karm*wstar**3.*z(iz)/hbl)**(1./3.) 
	    Pr=1+(Pr0-1)*exp(-alpha*(z(iz)-eps*hbl)**2./hbl**2.)
            difk(iz)= karm*ws*z(iz)*(1.-eps2*z(iz)/hbl)**2.
            dift(iz)= difk(iz)/Pr
         else
            dift(iz)=0.
            difk(iz)=0.
         endif
         if(z(iz+1).le.hbl.and.z(iz).ge.hbl*0.75) then
           dift3(iz) = 0.85*karm*(dR*hbl*g/thv)**(1./3.)*
     :               (z(iz)-hbl*0.75)**2.
     :              /(hbl-hbl*0.75)*(1.-(z(iz)-hbl*0.75)
     :              /(hbl-hbl*0.75))**0.5
           difk3(iz) = 0.75*dift3(iz)
         else 
           dift3(iz)=0.
           difk3(iz)=0.
         endif
         dift2(iz)=max(dift2(iz),dift(iz)+dift3(iz))
         difk2(iz)=max(difk2(iz),difk(iz)+difk3(iz))
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
          difk2(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
	    dift2(1)=difk(1)
	else
	  difk2(1)=karm*ws0*z_sl*(1.-z_sl/hbl)**2.
        dift2(1)=difk(1)/Pr0
      endif
	
      do iz = 2,nz
         dift_hl(iz)=0.5*(dift2(iz-1)+dift2(iz))
         difk_hl(iz)=0.5*(difk2(iz-1)+difk2(iz))
         dift_hl2(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl2(iz)=0.5*(difk(iz-1)+difk(iz))
         dift_hl3(iz)=dift_hl(iz)
         if (z(iz-1).lt.hbl.and.z(iz).gt.hbl) then 
            dift_hl3(iz)=0.
!          dift_hl(iz)=dift(iz)
!          dift_hl(iz-1)=dift(iz-2)
!         h3e(iz)=-wth_h
!         h3e(iz-1)=-wth_h/5.
         endif
      enddo

       do iz=2,nz-1
         h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
         h3(iz)=h3(iz)*0.5*(dift(iz)+dift(iz-1))
         if(z(iz).le.hbl) then
 !        zf=z(iz)
!	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
!     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
!	     w_2=w_2**(2./3.)
!	     gammah2=b_hm*wstar**2*tstar_f/w_2/hbl
	     h3c(iz)=-gammah*dift_hl(iz)
!	     h3e(iz+1)=-wth_h*min(1.,(z(iz)/hbl)**7.)
!             write(0,*)z(iz), h3(iz),h3e(iz)
	   else
           h3c(iz)=0.
           h3e(iz+1)=0.
	   endif     
       enddo

      !      if (implicit) then
      do iz=1,nz-1
      F(iz)=th(iz,1)  +dtl/(0.5*dz(iz)+dz(iz+1))*
     :        (gammah*(dift_hl2(iz)-dift_hl2(iz+1)))         ! right-hand side
      enddo
      F(nz)=th(nz,1)+dtl/dz(nz)*gammah*dift_hl2(nz)
      call implicit_dif(F,th(1:nz,1),dift_hl,ust_s*tst_s,difunt2)
      F(1:nz)=u(1:nz,1)
      call implicit_dif(F,u(1:nz,1),difk_hl,cdm*u(1,2),difunu2)
      F(1:nz)=v(1:nz,1)
      call implicit_dif(F,v(1:nz,1),difk_hl,cdm*v(1,2),difunv2)
      do iz=1,nz-1
      F(iz)=qv(iz,1)+dtl/(0.5*dz(iz)+dz(iz+1))*
     :        (gammaq*(dift_hl2(iz)-dift_hl2(iz+1)))         ! right-hand side
      enddo
      F(nz)=qv(nz,1)+dtl/dz(nz)*gammaq*dift_hl2(nz)
      call implicit_dif(F,qv(1:nz,1),dift_hl,ust_s*qst_s,difunqv2)
      F(1:nz)=qc(1:nz,1)
      call implicit_dif(F,qc(1:nz,1),dift_hl3,0.,difunqc2)
      F(1:nz)=qr(1:nz,1)
      call implicit_dif(F,qr(1:nz,1),dift_hl,0.,difunqr2)
!         else
!         endif

	do iz=2,nz
        def13(iz)=def13(iz)*0.5*(difk(iz)+difk(iz-1))
        def23(iz)=def23(iz)*0.5*(difk(iz)+difk(iz-1))
        if(z(iz-1).le.hbl.and.-tst_s.gt.0) then
	     def13c(iz)=-gamma_u*0.5*(difk(iz)+difk(iz-1))
	     def23c(iz)=-gamma_v*0.5*(difk(iz)+difk(iz-1))
	     !write(0,*) def13(iz),def13c(iz)
	  else
           def13c(iz)=0.
           def23c(iz)=0.
	  endif     
      enddo
      
      def13(1)=cdm*u(1,2)
      def23(1)=cdm*v(1,2)
      
      
      do iz=1,nz-1
        difunu(iz)=(def13(iz+1)+def13c(iz+1)-def13(iz)-def13c(iz))
     :   /(0.5*(dz(iz+1)+dz(iz)))
        difunv(iz)=(def23(iz+1)+def23c(iz+1)-def23(iz)-def23c(iz))
     :   /(0.5*(dz(iz+1)+dz(iz)))
      enddo
      
       do iz=2,nz-1
         h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
         h3(iz)=h3(iz)*0.5*(dift(iz)+dift(iz-1))
         if(z(iz-1).le.hbl.and.-tst_s.gt.0) then
         zf=z(iz)
	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
	     w_2=w_2**(2./3.)
	     gammah2=b_hm*wstar**2*tstar_f/w_2/hbl
	     h3c(iz)=-gammah*0.5*(dift(iz)+dift(iz-1))
!	     h3e(iz+1)=-wth_h*min(1.,0.5*(z(iz)+z(iz))/hbl)**3.
!             write(0,*)z(iz), h3(iz),h3e(iz)
             
	   else
           h3c(iz)=0.
	   endif
                if (z(iz).gt.hbl.and.z(iz-1).le.hbl) then
                   wth_h2=difunt2(iz)*dz(iz)!h3(iz)+h3c(iz)
                   endif
       enddo

        do iz=2,nz-1
         wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
         wq3(iz)=wq3(iz)*0.5*(dift(iz)+dift(iz-1))
!         if(z(iz).le.hbl.and.-(tst_s+0.61*thv*qst_s).gt.0) then
!           zf=z(iz-1)
!	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
!     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
!	     w_2=w_2**(2./3.)
!	     gammaq=b_hm*wstar*(-ust_s*qst_s)/w_2/hbl
	     wq3_c(iz)=-gammaq*0.5*(dift(iz)+dift(iz-1))
!             wq3e(iz)=-wq_h*min(1.,0.5*(z(iz)+z(iz-1))/hbl)**3.
!	   else
!           wq3_c(iz)=0.
!	   endif     
       enddo
      
      h3(1)=ust_s*tst_s
      h3c(1)=0.
      wq3(1)=ust_s*qst_s
      wq3_c(1)=0.
      
      if (ifwr.ne.0) then
        do iz=2,nz
          wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
          wq3c(iz)=wq3c(iz)*0.5*(dift(iz)+dift(iz-1))
          wq3r(iz)=(qr(iz,2)-qr(iz-1,2))/dz(iz)
          wq3r(iz)=wq3r(iz)*0.5*(dift(iz)+dift(iz-1))
        enddo
      endif
      
      do iz=1,nz-1
        difunt(iz)=(h3(iz+1)+h3c(iz+1)+h3e(iz+1)
     :             -h3(iz)-h3c(iz)-h3e(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
!        write(0,*)z(iz),difunt(iz),difunt2(iz)
         difunqv(iz)=(wq3(iz+1)+wq3_c(iz+1)+wq3e(iz+1)
     :             -wq3(iz)-wq3_c(iz)-wq3e(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
 !        write(0,*) difunqv(iz)
          if(z(iz).gt.hbl.and.z(iz-1).lt.hbl) then
            write(0,*) 'rad = ', rad(iz)
            write(0,*) 'advection =',-vat(iz)
            write(0,*) 'diff =', difunt2(iz)
            write(0,*) 'cond =',
     :        hlat/cp*condensat(iz)
            tende= difunt2(iz)+rad(iz)-vat(iz)+
     :         hlat/cp*condensat(iz)
            the = th(iz,2)
            endif
      enddo

      if (ifwr.ne.0) then
        do iz=1,nz-1
          difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
          difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
        enddo
      endif
      
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
