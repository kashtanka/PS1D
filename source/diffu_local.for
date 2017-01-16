      subroutine diffu_local
      use alloc_1d
      implicit none
      integer iz
      real,external:: qsat
      real*8 karm,lmax,mixl,zero,xn,dudz,dthdz,b, rich,def
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     : difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz),
     :     difunqsn2(1:nz),difunqci2(1:nz),xthl(1:nz),
     :     xqt(1:nz),fri_m,fri_h
      real*8 t_hl, pot_hl, alpha, thv, wstar, w_m3
      real, parameter:: B2 = 5.
      real, parameter:: b_0=6.5
      real, parameter:: al = 3.
      real*8 gammah,tstar_f,gammaq,wsm,ws0,eps,eps2,ws,
     :       Pr,Pr0,fi_m,fi_h, fri, surf_flux
      integer lambda_b, lambda_s, mix_length
      real*8 tstar,ustar,dzeta
      dift_hl=0.
      difk_hl=0.
      h3e = 0.
      wq3e=0.
      eps2=1.
!%%%%%%%%% LOCAL CLOSURE %%%%%%%%%%%%%%%%%%%
!  shear  - vertical wind shear
!  mixl   - mixing length (calc. using Blackadar f-la)
!  xn     - square of Brunt-Vaisala frequency
!  ri     - Richardson number
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
      b=0.00000001 ! 0.005  ! fixed value for F(Ri) for strongly-stable stratification
      mix_length = 1 ! 1 - Blackadar original; 2 - Sorbjan's modification
!-----------CONSERVATIVE VARIABLES--------------------------------!
      do iz = 1,nz
         xthl(iz) = th(iz,2)-th(iz,2)/t(iz)*hlatcp*qc(iz,2)
         xqt(iz) = qv(iz,2) + qc(iz,2)
      enddo

!-------asymptotic mixing length for use in Blackadar formula-----!
      if(frac*Fv+(1.-frac)*Fv2.lt.0) then          ! stable stratification
         if(mix_length.eq.2) then ! Sorbjan/Blackadar version
            if (frac.lt.1) then
               lmax = 0.009*(frac*ust_s + (1. - frac)*ust_s2)/fcor
            else
               lmax = 0.009*ust_s/fcor
            endif
         elseif (mix_length .eq. 1) then ! practical version
            lmax=40.       
         endif
      else                      ! unstable stratification
         lmax=100. !0.15*hbl
      endif
      
      call richardson    ! Richardson number and wind shear

      do iz = 2,nz-1
!---------mixing length---------------------!
         if (mix_length.eq.1) then
            mixl = karm*z(iz)/(1.+karm*z(iz)/lmax)
         elseif (mix_length.eq.2) then    ! Sorbjan modification
            lambda_b = lambda_s/rich
            mixl = karm*z(iz)/(1.+karm*z(iz)/lmax*(1. + lmax/lambda_b))
         endif
!-------------------------------------------!       
	 if(ri(iz).ge.0) then
            fri = max(b,(1.-5.*ri(iz)))**2.
            fri = (1.+5.*ri(iz)+44.*ri(iz)**2.)**(-2.)
            fri_m = (1.+300.*ri(iz)**2.)**(-3./2.)
            fri_h = 1./0.9/(1.+250.*ri(iz))**(3./2.)
            difk(iz)=mixl**2.*shear(iz)*fri_m
	    dift(iz)=mixl**2.*shear(iz)*fri_h !difk(iz)
!             write(0,*) z(iz),dift(iz),difk(iz)
	 endif
	 if(ri(iz).lt.0) then
	    difk(iz)=mixl**2.*shear(iz)*(min(9.,sqrt(1.-16.*ri(iz))))
	    dift(iz)=difk(iz)*(min(3.,(1.-16.*ri(iz))**0.25))
         endif
      enddo
!----------Km and Kh in the surface layer--------------!
         
      dzeta = frac*dzits + (1-frac)*dzits2
      ustar = frac*ust_s + (1-frac)*ust_s2
      tstar = frac*tst_s + (1-frac)*tst_s2
      if(-(frac*tst_s+(1-frac)*tst_s2).le.0) then    
         dudz=(1.+5.*dzeta)*ustar/(karm*z_sl)
         dthdz=(1.+5.*dzeta)*tstar/(karm*z_sl)
      else
         dudz=(1.-16.*dzeta)**(0.25)*ustar/(karm*z_sl)
         dthdz=(1.-16.*dzeta)**(0.25)*tstar/(karm*z_sl)
      endif
      mixl=karm*z_sl/(1.+karm*z_sl/lmax)
      def=dudz
      xn=g*dthdz/th(1,2)
      rich=xn/(def*def+zero)     
      if(-(frac*tst_s+(1-frac)*tst_s2).le.0) then
         difk(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
         dift(1)=difk(1)
      else
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

!--------------SURFACE FLUXES-------------------!
         
         if (frac.lt.1) then
            def13(1) = u(1,2)*(frac*cdm + (1.-frac)*cdm2)
            def23(1) = v(1,2)*(frac*cdm + (1.-frac)*cdm2)
            h3(1) = frac*ust_s*tst_s + (1.-frac)*ust_s2*tst_s2
            fthl(1) = frac*ust_s*tst_s + (1.-frac)*ust_s2*tst_s2
            wq3(1) = frac*ust_s*qst_s + (1.-frac)*ust_s2*qst_s2
            fqt(1) = frac*ust_s*qst_s + (1.-frac)*ust_s2*qst_s2
         else
            def13(1)=cdm*u(1,2) !surface flux of momentum
            def23(1)=cdm*v(1,2)
            fthl(1)=ust_s*tst_s ! surface flux of heat
            h3(1)=ust_s*tst_s
            fqt(1)=ust_s*qst_s  !surface flux of moisture
            wq3(1) = ust_s*qst_s
         endif
         wq3c(1) =0.
!-------------PROGONKA--------------------------------!
         F(1:nz)=xthl(1:nz)     ! right-hand side
!         F(1:nz)=th(1:nz,1)     ! right-hand side
         call implicit_dif(F,xthl(1:nz),dift_hl,fthl(1),difunt2)
         F(1:nz)=u(1:nz,1)
         call implicit_dif(F,u(1:nz,1),difk_hl,def13(1),difunu2)
         F(1:nz)=v(1:nz,1)
         call implicit_dif(F,v(1:nz,1),difk_hl,def23(1),difunv2)
         F(1:nz)= xqt(1:nz)
         call implicit_dif(F,xqt(1:nz),dift_hl,fqt(1),difunqv2)
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
         
         do iz = 2,nz
            def13(iz) = difunu2(iz-1)*0.5*(dz(iz)+dz(iz-1))+def13(iz-1)
            def23(iz) = difunv2(iz-1)*0.5*(dz(iz)+dz(iz-1))+def23(iz-1)
            h3(iz) = difunt2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + h3(iz-1)
c            fthl(iz) =  difunt2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + fthl(iz-1)
c            fqt(iz) = difunqv2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + fqt(iz-1)
!            wq3c(iz) = difunqc2(iz-1)*0.5*(dz(iz)+dz(iz-1)) + wq3c(iz-1)
!            t_hl = 0.5*(t(iz)+t(iz-1))
!            pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
!            fthl(iz) = h3(iz) - hlatcp*wq3c(iz)  !/pot_hl
!            fqt(iz) = wq3(iz) + wq3c(iz)
   !         write(0,*) z(iz),fthl(iz),hbl
         enddo
!---------------FLUXES OF NONCONSERVATIVE TH, Qv and Qc ----------!
!---------------following Deardorff, 1976-------------------------!
c         do iz = 2,nz
c            if(qc(iz-1,2).gt.0.and.qc(iz,2).gt.0)
c     :           then
c               t_hl = 0.5*(t(iz)+t(iz-1))
c               pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
c               alpha = hlat*0.5*(qsat(t(iz),p(iz,2))
c     :              +qsat(t(iz-1),p(iz-1,2)))/rv/t_hl/t_hl
c               h3(iz) = (fthl(iz)+hlatcp*fqt(iz)/pot_hl)
c     :              /(1.+hlatcp*alpha)
c               wq3(iz) = alpha*h3(iz)*pot_hl
c               wq3c(iz) = fqt(iz) - wq3(iz)
c            else
c               h3(iz) = fthl(iz)
c               wq3(iz) = fqt(iz)
c               wq3c(iz) = 0.
c            endif
c         enddo

c         do iz=1,nz-1
c            difunt(iz)= (h3(iz+1)-h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
c            if(qc(iz,2).eq.0) then
!     difunt(iz) = difunt2(iz)
c               difunqv(iz) = difunqv2(iz)
c               difunqc(iz) = 0.
c            else
               
!              pot_hl = t(iz)/(0.5*(th(iz,2)+th(iz-1,2)))
!     alpha = hlat*qsat(t(iz),p(iz,2))/rv/t(iz)/t(iz)
!     difunt(iz) = (difunt2(iz)+hlatcp*difunqv2(iz))    !1./pot_hl*
!     :                    /(1.+hlatcp*alpha)
!           difunt(iz)= (h3(iz+1)-h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
!           difunqv(iz) = alpha*difunt(iz)
c               difunqv(iz)=(wq3(iz+1)-wq3(iz))/(0.5*(dz(iz+1)+dz(iz)))
c               difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
c            endif
c         enddo
       
c         if (ifwr.ne.0) then
c            do iz=1,nz-1
!             if(qc(iz-1,2).eq.0) then
!                 difunqc(iz) = 0.
!             else
!                 difunqc(iz) = (difunqv2(iz) - difunqv(iz))
!             endif
c               difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
! difunqr(iz)=(wq3r(iz+1)-wq3r(iz))/(0.5*(dz(iz+1)+dz(iz)))
c            enddo
!          do iz=2,nz-2
!              if(qc(iz-1,2).ne.0.and.qc(iz,2).ne.0) then
!                 difunqc(iz) = 1./5.*(difunqc(iz-2)+
!  :                  difunqc(iz-1)+difunqc(iz)+
!  :                  difunqc(iz+1)+  difunqc(iz+2))
            
!              endif
            
!         enddo
c         endif

         difunu(1:nz)=difunu2
         difunt(1:nz)=difunt2    !
         difunv(1:nz)=difunv2
         difunqv(1:nz)=difunqv2   !
         difunqc(1:nz)=difunqc2   !
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
