      subroutine diffu_local
      use alloc_1d
      implicit none
      integer iz
      real,external:: qsat
      real*8 karm,def,lmax,mixl,rich,zero,xn,dudz,dthdz,b
      real*8 A(1:nz),BB(1:nz),C(1:nz),F(1:nz),dift_hl(1:nz)
      real*8 difunt2(1:nz)
      real*8 D1,Dn,kappa1,nu1,kappa2,nu2
      dift_hl=0.
      difunt2=0.
      karm=0.4
      zero=1.e-4
!      b=0.005
       b=0.
       
      do iz=2,nz
        def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
        def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
      enddo
       if(ust_s*tst_s.ge.0)then      ! stable stratification
         lmax=40.       
       else                          !  unstable stratification
         lmax=300.
       endif
      do iz=2,nz-1
        def=max(zero,dsqrt(((u(iz+1,2)-u(iz-1,2))
     :          /(z(iz+1)-z(iz-1)))**2.+
     :         ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.))
 !       def=dsqrt(((u(iz,2)-u(iz-1,2))/(z(iz)-z(iz-1)))**2.+
 !    :         ((v(iz,2)-v(iz-1,2))/(z(iz)-z(iz-1)))**2.)
        
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
!            xn=g*(1+hlat*qsat(t(iz),p(iz,2))/r/t(iz))
!     :      *(1+0.622*hlat**2.*qsat(t(iz),p(iz,2))/cp/r/t(iz)**2)**(-1)
!     :      *((log(th(iz+1,2))
!     :      -log(th(iz-1,2)))/(z(iz+1)-z(iz-1))+
!     :       hlat/cp/t(iz)*(qsat(t(iz+1),p(iz+1,2))-
!     :       qsat(t(iz-1),p(iz-1,2))))/(z(iz+1)-z(iz-1))-
!     :       g*(qv(iz+1,2)+qc(iz+1,2)-qv(iz-1,2)-qc(iz-1,2))/
!     :       (z(iz+1)-z(iz-1))
              xn=g*((th(iz+1,2)-hlat/cp*qc(iz+1,2))
     :        -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
     :        /(z(iz+1)-z(iz-1))/th(iz,2)
 !              xn=g*(th(iz+1,2)
 !    :      -th(iz-1,2))
 !    :      /(z(iz+1)-z(iz-1))/
 !    :     th(iz,2)
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
  !      xn=g*(th(iz,2)-th(iz-1,2))/(z(iz)-z(iz-1))/
  !   :   (0.5*(th(iz,2)+th(iz-1,2)))
 !       xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
        rich=xn/(def*def+zero)
        ri(iz)=rich
        
	  if(rich.ge.0) then
          difk(iz)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
	    dift(iz)=difk(iz)
	  endif
	  if(rich.lt.0) then
	    difk(iz)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))	
	    
        endif
        
      enddo

      if(-tst_s.le.0) then
      dudz=(1.+5.*dzits)*ust_s/(karm*z_sl)
      dthdz=(1.+5.*dzits)*tst_s/(karm*z_sl)
      endif
	if(-tst_s.gt.0) then
	dudz=(1.-16.*dzits)**(0.25)*ust_s/(karm*z_sl)
      dthdz=(1.-16.*dzits)**(0.25)*tst_s/(karm*z_sl)
      endif
      !write(0,*) dudz,dthdz,dzits
      mixl=karm*z_sl/(1.+karm*z_sl/lmax)
      def=dudz
      xn=g*dthdz/th(1,2)
      rich=xn/(def*def+zero)
      
      if(-tst_s.le.0) then
         difk(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
         dift(1)=difk(1)
      endif
      if(-tst_s.gt.0) then
         difk(1)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
         dift(1)=difk(1)*(min(3.,(1.-16.*rich)**0.25))	
      endif

!      if (implicit) then
! 1. heat
! boundary conditions for forward phase
      do iz = 2,nz
         dift_hl(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl(iz)=0.5*(difk(iz-1)+difk(iz))
!         write(0,*) dift_hl(iz)
      enddo
         D1=dtl*dift_hl(2)/(0.5*(dz(1)+dz(2))*dz(2))
         kappa1=D1/(D1+1.)
         nu1=-D1*dz(2)/dift_hl(2)/(D1+1.)*ust_s*tst_s+1./(D1+1.)*th(1,1)
! boundary conditions for backward phase
         Dn=dtl*dift_hl(nz)/(dz(nz)**2.)
         kappa2=Dn/(Dn+1.)
         nu2=th(nz,1)/(Dn+1.)
         F(1:nz)=th(1:nz,1)
! coefficients
         do iz = 2,nz-1
            A(iz)=dtl*dift_hl(iz)/(0.5*(dz(iz)+dz(iz+1))*dz(iz))
            BB(iz)=dtl*dift_hl(iz+1)/(0.5*(dz(iz)+dz(iz+1))*dz(iz+1))
            C(iz)=A(iz)+BB(iz)+1.
         enddo
         call progonka(A,BB,C,F,kappa1,nu1,kappa2,nu2,difunt2)
! temperature

! etc
!         else
!         endif

        
	!difk(1)=ust_s**2./dudz
	!dift(1)=tst_s*ust_s/dthdz
	!write(0,*)rich
  !    difk(1)=difk(2)
  !    dift(1)=dift(2)
      !write(0,*)cdm
      do iz=2,nz
        def13(iz)=def13(iz)*0.5*(difk(iz)+difk(iz-1))
        def23(iz)=def23(iz)*0.5*(difk(iz)+difk(iz-1))
 !        def13(iz)=def13(iz)*difk(iz)
 !        def23(iz)=def23(iz)*difk(iz)
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
        
        wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
        wq3(iz)=wq3(iz)*0.5*(dift(iz)+dift(iz-1))
              
      enddo
      
      if (ifwr.ne.0) then
        do iz=2,nz
          wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
          wq3c(iz)=wq3c(iz)*0.5*(dift(iz)+dift(iz-1))
        
          wq3r(iz)=(qr(iz,2)-qr(iz-1,2))/dz(iz)
          wq3r(iz)=wq3r(iz)*0.5*(dift(iz)+dift(iz-1))
        enddo
      endif
      
      if (ifmf.ne.0) then
      do iz=2,nz
          wq3ci(iz)=(qci(iz,2)-qci(iz-1,2))/dz(iz)
          wq3ci(iz)=wq3ci(iz)*0.5*(dift(iz)+dift(iz-1))
        
          wq3sn(iz)=(qsn(iz,2)-qsn(iz-1,2))/dz(iz)
          wq3sn(iz)=wq3sn(iz)*0.5*(dift(iz)+dift(iz-1))
        enddo
      endif
       
      
 !     if (qif.ne.0) then
 !        h3(1)=ust_s*tst_s 
 !     else
         h3(1)=ust_s*tst_s
 !     endif
      wq3(1)=ust_s*qst_s
      !write(0,*)ust_s*tst_s, hlat/cp*ust_s*qst_s
      
      do iz=1,nz-1
        difunt(iz)=(h3(iz+1)-h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
        difunqv(iz)=(wq3(iz+1)-wq3(iz))/(0.5*(dz(iz+1)+dz(iz)))
      enddo
       
      
      if (ifwr.ne.0) then
        do iz=1,nz-1
          difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
          difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
        enddo
      endif
      
      if (ifmf.ne.0) then
        do iz=1,nz-1
          difunqci(iz)=(wq3ci(iz+1)-wq3ci(iz))/(0.5*(dz(iz+1)+dz(iz)))
          difunqsn(iz)=(wq3sn(iz+1)-wq3sn(iz))/(0.5*(dz(iz+1)+dz(iz)))
        enddo
      endif
      
      do iz=1,nz
      write(0,*) iz,difunt(iz),difunt2(iz)
      enddo
      
      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
      difunqv(nz)=0.
      difunqc(nz)=0.
      difunqr(nz)=0.
      difunqci(nz)=0.
      difunqsn(nz)=0.
      
      end
