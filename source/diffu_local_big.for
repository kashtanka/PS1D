      subroutine diffu_local
      use alloc_1d
      implicit none
      integer iz
      real,external:: qsat
      real*8 karm,def,lmax,mixl,rich,zero,xn,dudz,dthdz,b
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     : difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz),
     :     difunqsn2(1:nz),difunqci2(1:nz),xthl(1:nz),
     :     xqt(1:nz)
      real*8 t_hl, pot_hl, alpha, thv, wstar, w_m3
      real, parameter:: B2 = 5.
      real, parameter:: b_0=6.5
      real, parameter:: al = 3.
      real*8 gammah,tstar_f,gammaq,wsm,ws0,eps,eps2,ws,
     :       Pr,Pr0,fi_m,fi_h, fri
      integer lambda_b, lambda_s, mix_length
      dift_hl=0.
      difk_hl=0.
      h3e = 0.
      wq3e=0.
      eps2=1.
!%%%%%%%%% LOCAL CLOSURE %%%%%%%%%%%%%%%%%%%
!  def    - vertical wind shear
!  mixl   - mixing length (calc. using Blackadar f-la)
!  xn     - square of Brunt-Vaisala frequency
!  rich   - Richardson number
!  difk   - eddy diffusivity for momentum
!  dift   - eddy diffusivity for scalars
!  FLUXES:
!  def13  - u-component of the vertical turb flux of momentum (u'w')
!  def23  - v-component of the vertical turb flux of momentum (v'w')
!  h3     - vertical turb flux of heat  (th'w')
!  wq3    - vertical turb flux of spec humidity (q'w')
!  wq3c   - vertical turb flux of cloud water   (qc'w')
!  wq3r   - vertical turb flux of rain water    (qr'w')
!  wq3ci  - vertical turb flux of cloud ice     (qci'w')
!  wq3sn  - vertical turb flux of snow          (qsn'w')
!  TENDENCIES (vertical divergence of fluxes):
!  difunu - d(u'w')/dz
!  difunv - d(v'w')/dz
!  difunt - d(th'w')/dz
!  difunq - d(q'w')/dz
!  etc.
      karm=0.4      ! von Karman constant
      zero=1.e-8    ! minumum vertical wind shear
      b=0.0 ! 0.005  ! fixed value for F(Ri) for strongly-stable stratification
      mix_length = 2 ! 1 - Blackadar original; 2 - Sorbjan's modification
!--------virtual potential temperature---------------------            
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
!----surface layer fraction      
      eps=z_sl/hbl
!---------surface Prandtl number-----------------!
       if (Fv.gt.0) then
        fi_m=(1.-7.*dzits)**(-1./3.)
        fi_h=(1.-16.*dzits)**(-1./2.)
      else
        fi_m=1+4.7*dzits
        fi_h=fi_m
      endif                                                    
      Pr0=fi_h/fi_m+b_0*eps*karm
!--------scaling parameters----------------------------------
      wstar=(g/thv*Fv*hbl)**(1./3.) 
      ws0 = (ust_s**3. + 7.*eps*karm*wstar**3.)**(1./3.)
      wsm=(ust_s**3. + 7.*karm*wstar**3.*0.5)**(1./3.)
      tstar_f=-ust_s*tst_s/wstar
      gammah=b_0*Fv/wsm/hbl
      gammaq=-b_0*qst_s*ust_s/wsm/hbl

!--------ENTRAINMENT----------------------------------
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
      w_m3=(wstar**3.+B2*ust_s**3.) +hbl*dR*g/thv
      wth_h= -6.*w_m3/hbl*2.     ! flux of virtual potential temperature
!      wth_h = (-1.6e-2)*2
!      if(nstep.gt.15000) then
      do iz = 2, nz-1
         if (0.5*(z(iz)+z(iz+1))
     :      .ge.hbl.and.0.5*(z(iz-1)+z(iz)).lt.hbl) then
            w_e = wth_h/(delta_th+0.61*th(iz,2)*delta_qv)
            wthl_h = -w_e*(delta_th - hlatcp*delta_qc)
            wqt_h = -w_e*(delta_qv) ! + delta_qc)
        endif
      enddo
 !     endif
      write(0,*) '_____________',u(2,2)
     
                                   ! entrainment velocity
!-----------CONSERVATIVE VARIABLES--------------------------------!
      do iz = 1,nz
         xthl(iz) = th(iz,1)-th(iz,1)/t(iz)*hlatcp*qc(iz,1)
         xqt(iz) = qv(iz,1) + qc(iz,1)
      enddo

!-------asymptotic mixing length for use in Blackadar formula-----!
      if(Fv.lt.0) then      ! stable stratification
         lmax=40.       
      else                  ! unstable stratification
         lmax=0.15*hbl
      endif

      do iz=2,nz-2
         def=dsqrt(((u(iz+1,1)-u(iz-1,1))
     :       /(z(iz+1)-z(iz-1)))**2.+
     :       ((v(iz+1,1)-v(iz-1,1))/(z(iz+1)-z(iz-1)))**2.)
 !       def=dsqrt(((u(iz,2)-u(iz-1,2))/(z(iz)-z(iz-1)))**2.+
 !    :         ((v(iz,2)-v(iz-1,2))/(z(iz)-z(iz-1)))**2.)
         if(qif.ne.0) then 
             if(ifwr.eq.0.) then              ! virtual theta
                xn=g*(th(iz+1,1)*(1.+0.61*qv(iz,1))
     :             -th(iz-1,1)*(1.+0.61*qv(iz,1)))
     :             /(z(iz+1)-z(iz-1))/
     :             th(iz,1)/(1.+0.61*qv(iz,1))
             else   
!                if(qv(iz,2)
!     :             .ge.0.9*qsat(t(iz),p(iz,2))) then
               if(qc(iz,1).gt.0.or.qc(iz+1,1).gt.0.or.qc(iz-1,1).gt.0)
     :          then
            xn=g*(1.+hlat*qsat(t(iz),p(iz,1))/r/t(iz))
     :    *(1.+0.622*hlat**2.*qsat(t(iz),p(iz,1))/cp/r/t(iz)**2.)**(-1.)
     :      *((log(th(iz+1,1))
     :      -log(th(iz-1,1)))/(z(iz+1)-z(iz-1))+
     :       hlat/cp/t(iz)*(qsat(t(iz+1),p(iz+1,1))-
     :       qsat(t(iz-1),p(iz-1,1))))/(z(iz+1)-z(iz-1))-
     :       g*(qv(iz+1,1)+qc(iz+1,1)-qv(iz-1,1)-qc(iz-1,1))/
     :       (z(iz+1)-z(iz-1))
!                    xn=g*((th(iz+1,2)-hlat/cp*qc(iz+1,2))
!     :                 -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
!     :                 /(z(iz+1)-z(iz-1))/th(iz,2)
 !              xn=g*(th(iz+1,2)
 !    :      -th(iz-1,2))
 !    :      /(z(iz+1)-z(iz-1))/
 !    :     th(iz,2)
                else
                   xn=g*(th(iz+1,1)*(1.+0.61*qv(iz+1,1))
     :                -th(iz-1,1)*(1.+0.61*qv(iz-1,1)))
     :                /(z(iz+1)-z(iz-1))/
     :                th(iz,1)/(1.+0.61*qv(iz,1))
                endif
            endif
         else
            xn=g*(th(iz+1,1)-th(iz-1,1))/(z(iz+1)-z(iz-1))/th(iz,1)
         endif
         rich=xn/(def*def+zero)
!---------mixing length---------------------!
         if (mix_length.eq.1) then
            mixl = karm*z(iz)/(1.+karm*z(iz)/lmax)
         elseif (mix_length.eq.2) then
            lambda_b = lambda_s/rich
            mixl = karm*z(iz)/(1.+karm*z(iz)/lmax*(1. + lmax/lambda_b))
         endif
!-------------------------------------------!
         ri(iz)=min(5.,rich)    ! ri(:) is used for output only       
	 if(rich.ge.0) then
            fri = max(b,(1.-5.*rich))**2.
            difk(iz)=mixl**2.*def*fri
	    dift(iz)=difk(iz)
!             write(0,*) z(iz),dift(iz),difk(iz)
	 endif
	 if(rich.lt.0) then
	    difk(iz)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
         endif
!-------------K-profile----------------------------------------!
          if(Fv.gt.0.and.0.5*(z(iz+1)+z(iz+2)).le.hbl) then
	    ws=(ust_s**3.+7.*karm*wstar**3.*z(iz)/hbl)**(1./3.) 
	    Pr=1+(Pr0-1)*exp(-al*(z(iz)-eps*hbl)**2./hbl**2.)
!            difk(iz)=karm*ws*z(iz)*(1.-eps2*z(iz)/hbl)**2.
!            dift(iz)= difk(iz)/Pr
          !  difk(iz) = max(difk(iz),difk2(iz))
          !  dift(iz) = max(dift(iz),dift2(iz))
          !  write(0,*) iz,difk(iz)
         else
!           dift(iz)=0.
!           difk(iz)=0.
         endif
          if(0.5*(z(iz+2)+z(iz+1)).le.hbl) then
           dift3(iz) = 0.85*karm*(dR*hbl*g/thv)**(1./3.)*
     :               (z(iz))**2.
     :              /(hbl)*(1.-(z(iz))
     :              /(hbl))**0.5
           difk3(iz) = 0.75*dift3(iz)
         else 
           dift3(iz)=0.
           difk3(iz)=0.
         endif
!         difk(iz)=difk(iz)+difk3(iz)
!         dift(iz)=dift(iz)+dift3(iz)
      enddo
!----------Km and Kh in the surface layer--------------!
      if(-tst_s.le.0) then
         dudz=(1.+5.*dzits)*ust_s/(karm*z_sl)
         dthdz=(1.+5.*dzits)*tst_s/(karm*z_sl)
      else
         dudz=(1.-16.*dzits)**(0.25)*ust_s/(karm*z_sl)
         dthdz=(1.-16.*dzits)**(0.25)*tst_s/(karm*z_sl)
      endif
      mixl=karm*z_sl/(1.+karm*z_sl/lmax)
      def=dudz
      xn=g*dthdz/th(1,2)
      rich=xn/(def*def+zero)     
      if(-tst_s.le.0) then
         difk(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
         dift(1)=difk(1)
      else
		write(0,*) 1.-16.*rich
		write(0,*) th(1,2),u(1,2),v(1,2)
         difk(1)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
         dift(1)=difk(1)*(min(3.,(1.-16.*rich)**0.25))	
      endif
!----------Km and Kh at half-levels-------------!
      do iz = 2,nz
         dift_hl(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl(iz)=0.5*(difk(iz-1)+difk(iz))
      enddo
!-------------implicit scheme---------------------------!
      if (implicit) then
         F(1:nz)=xthl(1:nz)     ! right-hand side
 !        do iz = 1,nz
 !           if (z(iz).le.hbl) then
 !           h3e(iz+1)=wthl_h*min(1.,(z(iz)/hbl)**7.)
 !           write(0,*) h3e(iz)
 !           endif
 !        enddo
 !        stop
 !        if (nstep.gt.2000) then
         do iz = 1, nz-1
            if (hbl.ge.0.5*(z(iz-1)+z(iz)).and.hbl.lt.
     :         0.5*(z(iz)+z(iz+1))) then
            F(iz) = xthl(iz) + dtl*(-wthl_h)
     :              /(0.5*(dz(iz+1)+dz(iz)))
            F(iz-1) = xthl(iz-1) + dtl*(wthl_h)
     :              /(0.5*(dz(iz)+dz(iz-1)))
            endif
         enddo
 !        endif
!         F(1:nz)=th(1:nz,1)     ! right-hand side
         call implicit_dif(F,xthl(1:nz),dift_hl,ust_s*tst_s,difunt2)
         F(1:nz)=u(1:nz,1)
         call implicit_dif(F,u(1:nz,1),difk_hl,cdm*u(1,2),difunu2)
         F(1:nz)=v(1:nz,1)
         call implicit_dif(F,v(1:nz,1),difk_hl,cdm*v(1,2),difunv2)
         F(1:nz)= xqt(1:nz)
!         if(nstep.gt.20000) then
!          do iz = 1, nz
!            if (hbl.gt.0.5*(z(iz-1)+z(iz)).and.hbl.lt.
!     :         0.5*(z(iz)+z(iz+1))) then
!            F(iz) = xqt(iz) + dtl*(-wqt_h)/(0.5*(dz(iz+1)+dz(iz)))
!            F(iz-1) = xqt(iz-1) + dtl*(wqt_h)/(0.5*(dz(iz)+dz(iz-1)))
!            endif
!         enddo
!         endif
         call implicit_dif(F,xqt(1:nz),dift_hl,ust_s*qst_s,difunqv2)
         if (ifwr.ne.0) then
            F(1:nz)=qc(1:nz,1)
            call implicit_dif(F,qc(1:nz,1),dift_hl,0.,difunqc2)
            F(1:nz)=qr(1:nz,1)
            call implicit_dif(F,qr(1:nz,1),dift_hl,0.,difunqr2)
         endif
!         if (ifmf.ne.0) then
!            F(1:nz)=qci(1:nz,1)
!            call implicit_dif(F,qci(1:nz,1),dift_hl,0.,difunqci2)
!            F(1:nz)=qsn(1:nz,1)
!            call implicit_dif(F,qsn(1:nz,1),dift_hl,0.,difunqsn2)
!         endif

!--------------FLUXES OF CONSERVATIVE VARIABLES-------------------!
         def13(1)=cdm*u(1,2) !surface flux of momentum
         def23(1)=cdm*v(1,2)
         fthl(1)=ust_s*tst_s ! surface flux of heat
         h3(1)=ust_s*tst_s
         fqt(1)=ust_s*qst_s !surface flux of moisture
         wq3(1)=ust_s*qst_s
         wq3c(1) =0.
         do iz = 2,nz
!            def13(iz) = difunu2(iz-1)*0.5*(dz(iz)+dz(iz-1))+def13(iz-1)
!            def23(iz) = difunv2(iz-1)*0.5*(dz(iz)+dz(iz-1))+def23(iz-1)
            fthl(iz) =  difunt2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + fthl(iz-1)
            fqt(iz) = difunqv2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + fqt(iz-1)
!            wq3c(iz) = difunqc2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + wq3c(iz-1)
!            t_hl = 0.5*(t(iz)+t(iz-1))
!            pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
!            fthl(iz) = h3(iz) - hlatcp*wq3c(iz)  !/pot_hl
!            fqt(iz) = wq3(iz) + wq3c(iz)
   !         write(0,*) z(iz),fthl(iz),hbl
         enddo
!---------------FLUXES OF NONCONSERVATIVE TH, Qv and Qc ----------!
!---------------following Deardorff, 1976-------------------------!
         do iz = 2,nz
            if(qc(iz-1,2).gt.0.and.qc(iz,2).gt.0)
     :          then
               t_hl = 0.5*(t(iz)+t(iz-1))
               pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
               alpha = hlat*0.5*(qsat(t(iz),p(iz,2))
     :                 +qsat(t(iz-1),p(iz-1,2)))/rv/t_hl/t_hl
              h3(iz) = (fthl(iz)+hlatcp*fqt(iz)/pot_hl)
     :                    /(1.+hlatcp*alpha)
               wq3(iz) = alpha*h3(iz)*pot_hl
                wq3c(iz) = fqt(iz) - wq3(iz)
            else
                h3(iz) = fthl(iz)
                wq3(iz) = fqt(iz)
                wq3c(iz) = 0.
            endif
         enddo

         do iz=1,nz-1
             difunt(iz)= (h3(iz+1)-h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
            if(qc(iz,2).eq.0) then
  !          difunt(iz) = difunt2(iz)
                difunqv(iz) = difunqv2(iz)
                difunqc(iz) = 0.
            else
               
 !              pot_hl = t(iz)/(0.5*(th(iz,2)+th(iz-1,2)))
!               alpha = hlat*qsat(t(iz),p(iz,2))/rv/t(iz)/t(iz)
!            difunt(iz) = (difunt2(iz)+hlatcp*difunqv2(iz))    !1./pot_hl*
!     :                    /(1.+hlatcp*alpha)
 !           difunt(iz)= (h3(iz+1)-h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
 !           difunqv(iz) = alpha*difunt(iz)
            difunqv(iz)=(wq3(iz+1)-wq3(iz))/(0.5*(dz(iz+1)+dz(iz)))
             difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
            endif
         enddo
       
         if (ifwr.ne.0) then
            do iz=1,nz-1
   !             if(qc(iz-1,2).eq.0) then
   !                 difunqc(iz) = 0.
   !             else
   !                 difunqc(iz) = (difunqv2(iz) - difunqv(iz))
   !             endif
               difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
              ! difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
            enddo
   !          do iz=2,nz-2
   !              if(qc(iz-1,2).ne.0.and.qc(iz,2).ne.0) then
   !                 difunqc(iz) = 1./5.*(difunqc(iz-2)+
   !  :                  difunqc(iz-1)+difunqc(iz)+
   !  :                  difunqc(iz+1)+  difunqc(iz+2))
              
   !              endif
 
   !         enddo
         endif

         difunu(1:nz)=difunu2
!         difunt(1:nz)=difunt2
         difunv(1:nz)=difunv2
 !        difunqv(1:nz)=difunqv2
!         difunqc(1:nz)=difunqc2
         difunqr(1:nz)=difunqr2
         difunqci(1:nz)=difunqci2
         difunqsn(1:nz)=difunqsn2
      
!------------------explicit scheme---------------------!
      else
         do iz=2,nz
            def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
            def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
            def13(iz)=def13(iz)*difk_hl(iz)   ! momentum flux
            def23(iz)=def23(iz)*difk_hl(iz)
         enddo
         def13(1)=cdm*u(1,2) !surface flux of momentum
         def23(1)=cdm*v(1,2)      
         do iz=1,nz-1
            difunu(iz)=(def13(iz+1)-def13(iz))/(0.5*(dz(iz+1)+dz(iz)))
            difunv(iz)=(def23(iz+1)-def23(iz))/(0.5*(dz(iz+1)+dz(iz)))
         enddo
         do iz=2,nz
            fthl(iz) = (xthl(iz)-xthl(iz-1))/dz(iz) 
            fthl(iz) = fthl(iz)*dift_hl(iz)    !  theta_l flux 
         !   if (nstep.gt.2000) then
            if(0.5*(z(iz)+z(iz+1)).le.hbl.and.
     :           0.5*(z(iz+1)+z(iz+2)).gt.hbl ) then
 !              h3e(iz+1)=wthl_h*min(1.,0.5*(z(iz+1)+z(iz))/hbl)**3.
 !              wq3e(iz+1)=wqt_h*min(1.,0.5*(z(iz+1)+z(iz))/hbl)**3.
                 h3e(iz+1) = wthl_h
                 wq3e(iz+1) = wqt_h
!                 write(0,*) wthl_h,fthl(iz),dift_hl(iz)
            endif
    !        endif
            fqt(iz) = (xqt(iz) - xqt(iz-1))/dz(iz)  
            fqt(iz) = fqt(iz)*dift_hl(iz)      ! q_total flux
            if (qc(iz,1).gt.0.and.qc(iz-1,1).gt.0) then
               t_hl = 0.5*(t(iz)+t(iz-1))
               pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
               alpha = hlat*0.5*(qsat(t(iz),p(iz,2))
     :                 +qsat(t(iz-1),p(iz-1,2)))/rv/t_hl/t_hl
              h3(iz) = ((fthl(iz)+h3e(iz))
     :                 +hlatcp*(fqt(iz)+wq3e(iz))/pot_hl)    !1./pot_hl*
     :                    /(1.+hlatcp*alpha)
               wq3(iz) = alpha*h3(iz)*pot_hl
                wq3c(iz) = (fqt(iz)+wq3e(iz)) - wq3(iz)
               else
                  h3(iz) = fthl(iz)+h3e(iz)
                  wq3(iz) = fqt(iz)+wq3e(iz)
                  wq3c(iz) = 0.
            endif
            if(0.5*(z(iz)+z(iz+1)).le.hbl.and.
     :           0.5*(z(iz+1)+z(iz+2)).gt.hbl ) then
               wq3c(iz+1) = 0.
            endif
            fthl(iz) = fthl(iz) + h3e(iz)
            fqt(iz) = fqt(iz) + wq3e(iz)
!            h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
!            h3(iz)=h3(iz)*dift_hl(iz)         ! heat flux
!            wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
!            wq3(iz)=wq3(iz)*dift_hl(iz)       ! moisture flux


         enddo
      
         if (ifwr.ne.0) then
            do iz=2,nz
 !              wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
 !              wq3c(iz)=wq3c(iz)*dift_hl(iz)  ! cloud water flux
               wq3r(iz)=(qr(iz,2)-qr(iz-1,2))/dz(iz)
               wq3r(iz)=wq3r(iz)*dift_hl(iz)  ! rain water flux
            enddo
         endif
    
         if (ifmf.ne.0) then
            do iz=2,nz
               wq3ci(iz)=(qci(iz,2)-qci(iz-1,2))/dz(iz)
               wq3ci(iz)=wq3ci(iz)*dift_hl(iz) ! cloud ice flux
               wq3sn(iz)=(qsn(iz,2)-qsn(iz-1,2))/dz(iz)
               wq3sn(iz)=wq3sn(iz)*dift_hl(iz) ! snow flux
            enddo
         endif
       
         h3(1)=ust_s*tst_s ! surface flux of heat
         wq3(1)=ust_s*qst_s !surface flux of moisture
      
         do iz=1,nz-1
            difunt(iz)=(h3(iz+1)    !+h3e(iz+1)
     :          -h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
            difunqv(iz)=(wq3(iz+1)   !+wq3e(iz+1)
     :          -wq3(iz))/(0.5*(dz(iz+1)+dz(iz)))
         enddo
       
         if (ifwr.ne.0) then
            do iz=1,nz-1
               difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
               difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
            enddo
         endif
      
         if (ifmf.ne.0) then
            do iz=1,nz-1
               difunqci(iz)=(wq3ci(iz+1)-wq3ci(iz))
     :                      /(0.5*(dz(iz+1)+dz(iz)))
               difunqsn(iz)=(wq3sn(iz+1)-wq3sn(iz))
     :                      /(0.5*(dz(iz+1)+dz(iz)))
            enddo
         endif

      endif
      
      
      
!      if(mod(nstep*dt,3600.).eq.0) then
!         def13(1)=cdm*u(1,2) !surface flux of momentum
!         def23(1)=cdm*v(1,2)
!         h3(1)=ust_s*tst_s ! surface flux of heat
!         wq3(1)=ust_s*qst_s !surface flux of moisture
!!         do iz = 2,nz
!            h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
!            h3(iz)=h3(iz)*dift_hl(iz)
!            wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
!            wq3(iz)=wq3(iz)*dift_hl(iz)
!            wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
!            wq3c(iz)=wq3c(iz)*dift_hl(iz)
!            def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
!            def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
!            def13(iz)=def13(iz)*difk_hl(iz)   ! momentum flux
!            def23(iz)=def23(iz)*difk_hl(iz)
!         enddo
!      endif
      


      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
      difunqv(nz)=0.
      difunqc(nz)=0.
      difunqr(nz)=0.
      difunqci(nz)=0.
      difunqsn(nz)=0.
      
      end
