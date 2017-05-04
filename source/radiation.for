      subroutine radiation
      use alloc_1d
      use goddard_mod
      use ice_mod,only: Tsi
      implicit none
      integer iz
      real*8,parameter :: xkp = 85.
      real*8,parameter :: DIV = 3.75e-6
      real*8,parameter :: F0 = 70.
      real*8,parameter :: F1 = 22.
      real*8,parameter:: gq = 1.7e-6
      real*8 zi,Q1,Q0,Q2,R3,mdR,dR2,LWP,LWP2,qmax,dRmax
      real*8 Fplus
      real*8 p_lw(nz+1), t_lw(nz), qv_lw(nz)
      real*8 coolr(nz), coolr2(nz),lwdif2(nz),xqv(nz)
      real,external:: qsati

      
      
      if (rad_par.eq.1) then
  
!-----------------------------------------!
!-!    net longwave radiative flux for  !-!
!-!    Sc experiment at half-levels     !-!
!-----------------------------------------!
      do iz=1,nz-1
      Q0=Q0+xkp*ro(iz)*qc(iz,2)*0.5*(dz(iz)+dz(iz+1))
      enddo

      do iz=2,nz-1
!--------determine qmax for LWP calc--------------!
        if (z(iz).le.hbl) then
           LWP = LWP+ro(iz)*qc(iz,2)*1000.*dz(iz)
           if (qc(iz,2).gt.qmax) qmax = qc(iz,2)
        endif       
         if(qv(iz,2)+qc(iz,2).ge.8.e-3.and.
     :      qv(iz+1,2)+qc(iz+1,2).lt.8.e-3) then
            zi=0.5*(z(iz)+z(iz+1))
            !zi = hbl
            endif
            zi=max(100.,zi)
      enddo
      LWP2=ro(1)*qmax**2./(2.*gq)*1000.
!      zi=800.
      Q1=Q0
      rfl(0)=F0*exp(-Q1)+F1
      do iz=2,nz-1
         Q1=Q1-xkp*ro(iz)*qc(iz-1,2)*0.5*(dz(iz)+dz(iz+1))
         Q2=Q2+xkp*ro(iz)*qc(iz-1,2)*0.5*(dz(iz)+dz(iz+1))
         if(0.5*(z(iz)+z(iz-1)).ge.zi) then
         R3=0.5*(ro(iz)+ro(iz-1))
     :           *cp*DIV*((0.5*(z(iz-2)+z(iz-1))-zi)**(4./3.)/4.
     :           +zi*(0.5*(z(iz-1)+z(iz-2))-zi)**(1./3.))
         else
         R3=0.
         endif
         
         rfl(iz)=F0*exp(-Q1)+F1*exp(-Q2) +R3
      enddo
      dR=0.
      mdR=0.
      do iz=1,nz-1
        rad(iz)=-(rfl(iz+1)-rfl(iz))/(0.5*(dz(iz+1)+dz(iz)))/cp/ro(iz)
!         dR=(rfl(iz)-rfl(iz-1))/cp/ro(iz)
!         if (dR.gt.mdR) mdR=dR
         if (0.5*(z(iz)+z(iz-1)).lt.hbl.and.
     :            0.5*(z(iz)+z(iz+1)).gt.hbl) then
            dR=(rfl(iz)-rfl(iz-1))/cp/ro(iz)
            Fplus = rfl(iz)
         endif
      enddo
      dRmax = 0.
      do iz=1,nz-1
         if(qc(iz,2).gt.0) then
            if ((Fplus - rfl(iz))/cp/ro(iz).gt.dRmax) 
     :          dRmax = (Fplus - rfl(iz))/cp/ro(iz)
         endif
      enddo
      dRmax = 55/cp/ro(1)
       dR2 = Fplus*(1. - exp(-0.03*sqrt(LWP)))
       dR=dRmax
       
       elseif (rad_par.eq.2) then
!-------------SOME INPUT PARAMETERS-----------------_!
          ta_s = 0.5*(t(1)+Tsi) !Tsi !0.5*(Tsi + t(1)) ! surface air temperature, K
          oa(1:nz) = 5.e-7 ! ozone mixing ratio by mass, g/g
          co2 = 300.e-6 ! co2 mixing ratio by volume, pppv
          high = .true.  ! option (.true. slower but more accurate)
          trace = .false.  ! option (absorbtion in some minor bands)
          n2o = 0. ! n2o mixing ratio by volume, pppv
          ch4 = 0. ! ch4 mixing ratio by volume, pppv
          cfc11 = 0.  !
          cfc12 = 0.  !
          cfc22 = 0.  !
          vege = .false. ! If vege=.true., a vegetation layer is added
          fs(1,1) = 1. ! fractional cover of subgrid regions
          fs(1,2) = 0. !
          tsurfs(1,1) = Tsi !Tsi !Tsi ! land or ocean surface temperature
          write(0,*) 'tsurfs =', tsurfs(1,1)
          tsurfs(1,2) = 271.35 !271.35 ! land or ocean surface temperature
          eg(1,1:2,1:10) = 0.98 ! land or ocean surface emissivity
          ev(1,1:2,1:10) = 0.9 ! vegetation emissivity
          rvir(1,1:2,1:10) = 0.1 ! vegetation reflectivity
          overcast_lw = .true. ! 
          cldwater_lw = .true. !
          cwc(1,:,1) = qci(:,2)
          cwc(1,:,2) = qc(:,2)
          cwc(1,:,3) = qr(:,2)
          taucl(1,:,:) = 0.     ! cloud optical thickness
          fcld(1,1:nz) = 0. ! cloud fraction
          ict = nz !level index separating high and middle clouds
          icb = nz ! level index separating middle and low  clouds 
          aerosol = .false. ! option to include aerosols
!          na = 0  ! number of aerosol types 
          taual_lw(1,:,:,:) = 0 ! aerosol optical thickness
          ssaal_lw(1,:,:,:) = 0 ! aerosol single-scat albedo
          asyal_lw(1,:,:,:) = 0 ! aerosol asymmetry factor
 !------------flip atmospheric profiles-------!
          do iz = 1,nz
             xqv(iz) = 0.8*qsati(t(iz),p(iz,2))
          enddo
          do iz = 2,nz+1
             p_lw(iz) = phl(nz+2-iz)/100.
          enddo 
          p_lw(1) = 0.
          do iz = 2, nz
             t_lw(iz) = t(nz+1-iz)
             qv_lw(iz) = xqv(nz+1-iz) !qv(nz+1-iz,2)
          enddo 
          t_lw(1) = t_lw(2)
          qv_lw(1) = xqv(2)  !qv_lw(2)
          
          if(nstep.eq.1. or. mod(nstep,400).eq.0) then
         call irrad(1,nz,p_lw,t_lw,qv_lw,oa,
     &             ta_s,co2,high,trace,n2o,ch4,cfc11,cfc12,cfc22,
     &             vege,nsur,fs,tsurfs(1,1:2),
     &             eg,tsurfs(1,1:2),ev,rvir,
     &             overcast_lw,cldwater_lw,cwc,taucl,fcld,ict,icb,
     &             aerosol,na,taual_lw,ssaal_lw,asyal_lw, 
     &             flx_lw,flc_lw,dfdts,sfcem)
!         write(0,*) flx_lw(1,1:nz)
         write(0,*) 'sfc emmis', sfcem 
         LW = flx_lw(1,nz+1)
!         write(0,*) flx_lw(1,nz:1)
         do iz = 1,nz+1
            flx_lw2(iz) = flx_lw(1,nz+2-iz)
         enddo
         do iz = 1,nz
            coolr(iz) = (flx_lw(1,iz+1) - flx_lw(1,iz))*8.44/
     :                   (p_lw(iz+1)-p_lw(iz))
            coolr2(iz) = (flx_lw(1,iz+1) - flx_lw(1,iz))/
     :                    (ro(nz+1-iz)*cp)/
     :                   (z(nz+1-iz+1)-z(nz+1 -iz))*24*3600
            lwdif2(iz) = (flx_lw(1,iz+1) - flx_lw(1,iz))
     :                     /(ro(nz+1-iz)*cp)/
     :                   dz(nz+1-iz)
            ! if (iz.gt.200) then
            ! write(0,*)t_lw(iz),coolr2(iz),flx_lw(1,iz+1)
           !  endif
         enddo
         do iz = 1,nz
           if (z(iz)-z_sl.lt.1000) then
            lwdif(iz) =  lwdif2(nz+1-iz)
            else
               lwdif(iz) = 0.
            endif
         enddo
         endif
 !        do iz = 1,nz
 !           if (iz.le.200) then
 !           write(0,*)(flx_lw2(iz+1)-flx_lw2(iz))/(ro(nz+1-iz)*cp)*(3600
 !    :       *24)/dz(nz+1-iz),flx_lw2(iz+1),iz !(flx_lw(1,iz)-flx_lw(1,iz+1))/dz(iz)*dtl
 !           endif
 !        enddo
!         stop
      endif
      
      end
