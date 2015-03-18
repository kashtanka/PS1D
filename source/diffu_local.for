      subroutine diffu_local
      use alloc_1d
      implicit none
      integer iz
      real,external:: qsat
      real*8 karm,def,lmax,mixl,rich,zero,xn,dudz,dthdz,b
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     : difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz),
     :     difunqsn2(1:nz),difunqci2(1:nz)
      dift_hl=0.
      difk_hl=0.
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
      zero=1.e-3    ! minumum vertical wind shear
      b=0.00 ! 0.005  ! fixed value for F(Ri) for strongly-stable stratification

!-------asymptotic mixing length for use in Blackadar formula-----!
      if(Fv.lt.0) then      ! stable stratification
         lmax=40.       
      else                  ! unstable stratification
         lmax=0.15*hbl
      endif

      do iz=2,nz-1
         def=dsqrt(((u(iz+1,2)-u(iz-1,2))
     :       /(z(iz+1)-z(iz-1)))**2.+
     :       ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.)
 !       def=dsqrt(((u(iz,2)-u(iz-1,2))/(z(iz)-z(iz-1)))**2.+
 !    :         ((v(iz,2)-v(iz-1,2))/(z(iz)-z(iz-1)))**2.)
         mixl=karm*z(iz)/(1.+karm*z(iz)/lmax)
         if(qif.ne.0) then 
             if(ifwr.eq.0.) then              ! virtual theta
                xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2))
     :             -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :             /(z(iz+1)-z(iz-1))/
     :             th(iz,2)/(1.+0.61*qv(iz,2))
             else   
                if(qv(iz,2)
     :             .ge.0.9*qsat(t(iz),p(iz,2))) then        
!            xn=g*(1+hlat*qsat(t(iz),p(iz,2))/r/t(iz))
!     :      *(1+0.622*hlat**2.*qsat(t(iz),p(iz,2))/cp/r/t(iz)**2)**(-1)
!     :      *((log(th(iz+1,2))
!     :      -log(th(iz-1,2)))/(z(iz+1)-z(iz-1))+
!     :       hlat/cp/t(iz)*(qsat(t(iz+1),p(iz+1,2))-
!     :       qsat(t(iz-1),p(iz-1,2))))/(z(iz+1)-z(iz-1))-
!     :       g*(qv(iz+1,2)+qc(iz+1,2)-qv(iz-1,2)-qc(iz-1,2))/
!     :       (z(iz+1)-z(iz-1))
                    xn=g*((th(iz+1,2)-hlat/cp*qc(iz+1,2))
     :                 -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
     :                 /(z(iz+1)-z(iz-1))/th(iz,2)
 !              xn=g*(th(iz+1,2)
 !    :      -th(iz-1,2))
 !    :      /(z(iz+1)-z(iz-1))/
 !    :     th(iz,2)
                else
                   xn=g*(th(iz+1,2)*(1.+0.61*qv(iz+1,2))
     :                -th(iz-1,2)*(1.+0.61*qv(iz-1,2)))
     :                /(z(iz+1)-z(iz-1))/
     :                th(iz,2)/(1.+0.61*qv(iz,2))
                endif
            endif
         else
            xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
         endif
         rich=xn/(def*def+zero)
         ri(iz)=min(5.,rich)    ! ri(:) is used for output only       
	 if(rich.ge.0) then
            difk(iz)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
	    dift(iz)=difk(iz)
!             write(0,*) z(iz),dift(iz),difk(iz)
	 endif
	 if(rich.lt.0) then
	    difk(iz)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift(iz)=difk(iz)*(min(3.,(1.-16.*rich)**0.25))
         endif
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
         F(1:nz)=th(1:nz,1)     ! right-hand side
         call implicit_dif(F,th(1:nz,1),dift_hl,ust_s*tst_s,difunt2)
         F(1:nz)=u(1:nz,1)
         call implicit_dif(F,u(1:nz,1),difk_hl,cdm*u(1,2),difunu2)
         F(1:nz)=v(1:nz,1)
         call implicit_dif(F,v(1:nz,1),difk_hl,cdm*v(1,2),difunv2)
         F(1:nz)=qv(1:nz,1)
         call implicit_dif(F,qv(1:nz,1),dift_hl,ust_s*qst_s,difunqv2)
         if (ifwr.ne.0) then
            F(1:nz)=qc(1:nz,1)
            call implicit_dif(F,qc(1:nz,1),dift_hl,0.,difunqc2)
            F(1:nz)=qr(1:nz,1)
            call implicit_dif(F,qr(1:nz,1),dift_hl,0.,difunqr2)
         endif
         if (ifmf.ne.0) then
            F(1:nz)=qci(1:nz,1)
            call implicit_dif(F,qci(1:nz,1),dift_hl,0.,difunqci2)
            F(1:nz)=qsn(1:nz,1)
            call implicit_dif(F,qsn(1:nz,1),dift_hl,0.,difunqsn2)
         endif
         difunu(1:nz)=difunu2
         difunt(1:nz)=difunt2
         difunv(1:nz)=difunv2
         difunqv(1:nz)=difunqv2
         difunqc(1:nz)=difunqc2
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
            h3(iz)=(th(iz,2)-th(iz-1,2))/dz(iz)
            h3(iz)=h3(iz)*dift_hl(iz)         ! heat flux
            wq3(iz)=(qv(iz,2)-qv(iz-1,2))/dz(iz)
            wq3(iz)=wq3(iz)*dift_hl(iz)       ! moisture flux
         enddo
      
         if (ifwr.ne.0) then
            do iz=2,nz
               wq3c(iz)=(qc(iz,2)-qc(iz-1,2))/dz(iz)
               wq3c(iz)=wq3c(iz)*dift_hl(iz)  ! cloud water flux
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
               difunqci(iz)=(wq3ci(iz+1)-wq3ci(iz))
     :                      /(0.5*(dz(iz+1)+dz(iz)))
               difunqsn(iz)=(wq3sn(iz+1)-wq3sn(iz))
     :                      /(0.5*(dz(iz+1)+dz(iz)))
            enddo
         endif

      endif
      


      difunu(nz)=0.
      difunv(nz)=0.
      difunt(nz)=0.
      difunqv(nz)=0.
      difunqc(nz)=0.
      difunqr(nz)=0.
      difunqci(nz)=0.
      difunqsn(nz)=0.
      
      end
