      subroutine diffu_INM
      use alloc_1d
      implicit none
      real*8 lmax,lmax_h,mixl,xn,def,rich,mixl_h,F_ml
      real*8 Fm,Fh,dudz,dthdz
      real*8 F(1:nz),dift_hl(1:nz),difk_hl(1:nz)
      real*8 difunt2(1:nz),difunu2(1:nz),difunv2(1:nz),
     :       difunqv2(1:nz),difunqc2(1:nz),difunqr2(1:nz),
     :       difunqci2(1:nz),difunqsn2(1:nz)
      integer iz    
      real, external:: qsat
      real, parameter:: karm = 0.4
      real, parameter:: d = 5.
      real, parameter:: zero = 1.e-8
      real, parameter:: b = 5.
      real, parameter:: c = b**2./sqrt(3.)
      real, parameter:: ricr = 2./(3.*d)    ! limit for Richardson number in stable stratification
      dift_hl=0.
      difk_hl=0.
!%%%%%%%%% INM LOCAL CLOSURE %%%%%%%%%%%%%%%%%%%
!  def    - vertical wind shear
!  mixl   - mixing length
!  mixl_h - mixing length for scalars
!  xn     - square of Brunt-Vaisala frequency
!  rich   - Richardson number
!  ricr   - limiter for Richardson number
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


!--------- asymptotic length scale lmax at z --> inf
      lmax = 160.                    ! for momentum
      lmax_h = sqrt(3.*d/2.)*lmax    ! for heat and scalars

      do iz=2,nz-1
        def=max(zero,dsqrt(((u(iz+1,2)-u(iz-1,2))
     :          /(z(iz+1)-z(iz-1)))**2.+
     :         ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.))
 !       def=dsqrt(((u(iz,2)-u(iz-1,2))/(z(iz)-z(iz-1)))**2.+
 !    :         ((v(iz,2)-v(iz-1,2))/(z(iz)-z(iz-1)))**2.)

        if(qif.ne.0) then 
           if(ifwr.eq.0.) then              ! virtual theta
              xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2)) 
     :           -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :           /(z(iz+1)-z(iz-1))/
     :           th(iz,2)/(1.+0.61*qv(iz,2))
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
     :              -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
     :              /(z(iz+1)-z(iz-1))/th(iz,2)
               else
                    xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2))
     :                 -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :                 /(z(iz+1)-z(iz-1))/
     :                 th(iz,2)/(1.+0.61*qv(iz,2))
               endif
            endif
         else
            xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
         endif
!--------Richardson number------------------!
        rich=min(ricr,xn/(def*def+zero))
        ri(iz)=rich
!-------------- mixing length---------------! 
        if (rich.gt.0.) then
           F_ml = max(0.,1. - z(iz-1)/hbl)    !  stable stratification
        else
           F_ml = 1.                        !  neutral and unstable
        endif
        mixl = karm*z(iz)/(1.+karm*z(iz)/lmax)*F_ml
        mixl_h = karm*z(iz)/(1.+karm*z(iz)/lmax_h)*F_ml
!---------------eddy diffusivities----------!
        if (rich.ge.0.) then
           Fm = 1./(1.+2.*b*rich/sqrt(1+d*rich))
           Fh = 1./(1.+3.*b*rich*sqrt(1+d*rich))
        else
           Fm = 1.-2.*b*rich/(1+c*(mixl/z(iz))**2.*sqrt(-rich))
           Fh = 1.-3.*b*rich/(1+c*(mixl_h/z(iz))**2.*sqrt(-rich))
        endif
           difk(iz) = mixl**2.*def*Fm
           dift(iz) = mixl_h**2.*def*Fh
!--------------limiters----------------------!
!        if (xn/g*th(iz,2).gt.7.e-3.and.z(iz).lt.hbl) then
!           dift(iz) = 0.
!        endif
      enddo
!---------------eddy diffusivities at lowest level------------! 
      if(Fv.le.0) then
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
      
      if(Fv.le.0) then
         difk(1)=mixl**2.*def*(max(b,(1.-5.*rich))**2.)
         dift(1)=difk(1)
      else
         difk(1)=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
         dift(1)=difk(1)*(min(3.,(1.-16.*rich)**0.25))	
      endif
!-----------eddy diffusivities at half-levels-----------------!
      do iz = 2,nz
         dift_hl(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl(iz)=0.5*(difk(iz-1)+difk(iz))
      enddo
!---------------------implicit scheme-------------------------!
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
!------------------explicit scheme---------------------!
         difunu(1:nz)=difunu2
         difunt(1:nz)=difunt2
         difunv(1:nz)=difunv2
         difunqv(1:nz)=difunqv2
         difunqc(1:nz)=difunqc2
         difunqr(1:nz)=difunqr2
         difunqci(1:nz)=difunqci2
         difunqsn(1:nz)=difunqsn2
      else
         write(0,*) 'explicit scheme for INM closure is not used'
         stop
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
