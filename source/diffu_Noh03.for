      subroutine diffu_Noh03
      use alloc_1d,only : th,qv,u,v,dift,difk,Fv,dzits,
     :                    ust_s,tst_s,nz,z,dz,def13,def23,
     :                    g,hbl,z_sl,qif,h3,h3c,difunt,difunv,difunu,
     :                    cdm,h3e,def13c,def23c,difunqv,difunqc,difunqr,
     :                    qst_s,wq3,wq3c,wq3r,wq3_c,qc,qr,ifwr,cp,hlat,
     :                    t,p,dift2,difk2,dR,wq3e,dift3,difk3,ro
      implicit none
      integer iz
      real*8 def, ws0,wstar,eps,thv,fi_m,fi_h,Pr0,gammah
      real*8 lmax,mixl,xn,rich,dthdz,dudz,wth_h,w_m3,ws,Pr
      real*8 delta,w_2,tstar_f,zf,gammah2,wsm,gamma_u,gamma_v
      real*8 gammaq,wq_h
      real, external :: qsat
      real, parameter:: karm=0.4
      real, parameter:: b_0=6.5
      real, parameter:: b_hm=3.
      real, parameter:: b = 0.
      real, parameter:: zero = 1.e-8
      real, parameter:: B2 = 5.
      real, parameter:: A = 4.5
      real, parameter:: alpha = 3.
      real, parameter:: Sm=15.9
            
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
      
      
      
      eps=z_sl/hbl
!      write(0,*) eps
      wstar=(g/thv*Fv*hbl)**(1./3.)
      ws0 = (ust_s**3. + 7.*eps*karm*wstar**3.)**(1./3.)
      wsm=(ust_s**3. + 7.*karm*wstar**3.*0.5)**(1./3.)
      gammah=b_0*Fv/wsm/hbl
      gammaq=-b_0*qst_s*ust_s/wsm/hbl
!      gamma_u=-Sm*ust_s**2./wsm/hbl*(wstar/wsm)**3.
!     :        *u(1,2)/sqrt(u(1,2)**2.+v(1,2)**2.)
!      gamma_v=-Sm*ust_s**2./wsm/hbl*(wstar/wsm)**3.
!     :        *v(1,2)/sqrt(u(1,2)**2.+v(1,2)**2.)
      w_m3=wstar**3.+B2*ust_s**3.+hbl*dR
!      write(0,*) w_m3,150.*dR
      wth_h= 0. !-0.2*w_m3/hbl -0.2*dR
      wq_h=wth_h/9.*(-7.5e-3)
      write(0,*) wth_h/10*100
!      write(0,*)-A*w_m3/hbl,-0.2*dR
      
      delta=0.01*hbl
      tstar_f=-ust_s*tst_s/wstar
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
	      ws=(ust_s**3.+7.*karm*wstar**3.*z(iz)/hbl)**(1./3.) 
	      Pr=1+(Pr0-1)*exp(-alpha*(z(iz)-eps*hbl)**2./hbl**2.)
            difk(iz)= karm*ws*z(iz)*(1.-z(iz)/hbl)**2.
            dift(iz)= difk(iz)/Pr
            
!            dift(iz)=2.*dift(iz) !max(dift(iz),dift2(iz))
!	    difk(iz)=2.*difk(iz) !max(difk(iz),difk2(iz))

            if(z(iz).le.hbl) then
           dift3(iz)=0.85*karm*(min(1.,dR*hbl))**(1./3.)*(z(iz)-0)**2.
     :              /(hbl-0)*(1.-(z(iz)-0)/(hbl-0))**0.5
           difk3(iz)=0.75*dift3(iz)
         !  write(0,*) dift3(iz),dift(iz),(dR*hbl)**(1./3.)
           endif

           dift(iz)=max(dift(iz),dift3(iz))
	    difk(iz)=max(difk(iz),difk3(iz))

            
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
         if(z(iz).le.hbl.and.-tst_s.gt.0) then
         zf=z(iz)
	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
	     w_2=w_2**(2./3.)
	     gammah2=b_hm*wstar**2*tstar_f/w_2/hbl
	     h3c(iz)=-gammah*0.5*(dift(iz)+dift(iz-1))
	     h3e(iz)=-wth_h*min(1.,0.5*(z(iz)+z(iz))/hbl)**3.
!             write(0,*)z(iz), h3(iz),h3e(iz)
	   else
           h3c(iz)=0.
	   endif     
       enddo

        do iz=2,nz-1
         wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
         wq3(iz)=wq3(iz)*0.5*(dift(iz)+dift(iz-1))
         if(z(iz).le.hbl.and.-(tst_s+0.61*thv*qst_s).gt.0) then
!           zf=z(iz-1)
!	     w_2=(1.6*ust_s**2.*(1.-zf/hbl))**(3./2.)+
!     :         1.2*wstar**3.*zf/hbl*(1.-0.9*zf/hbl)**(3./2.)
!	     w_2=w_2**(2./3.)
!	     gammaq=b_hm*wstar*(-ust_s*qst_s)/w_2/hbl
	     wq3_c(iz)=-gammaq*0.5*(dift(iz)+dift(iz-1))
!             wq3e(iz)=-wq_h*min(1.,0.5*(z(iz)+z(iz-1))/hbl)**3.
	   else
           wq3_c(iz)=0.
	   endif     
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
         difunqv(iz)=(wq3(iz+1)+wq3_c(iz+1)+wq3e(iz+1)
     :             -wq3(iz)-wq3_c(iz)-wq3e(iz))
     :             /(0.5*(dz(iz+1)+dz(iz)))
 !        write(0,*) difunqv(iz)
      enddo

      if (ifwr.ne.0) then
        do iz=1,nz-1
          difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
          difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
        enddo
      endif
      
      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
      difunqv(nz)=0.
      difunqc(nz)=0.
      difunqr(nz)=0.
	
	
      end
