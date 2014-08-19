      program nh1d
      use alloc_1d
      implicit real*8 (a-h,o-z)
      logical scrout
	logical difloc,nonloc_tm86,nonloc_ls96,nonloc_noh03,loc_inm,
     :          nonloc_lock
	character*80 fnmap,fnfor,fext,fnout,fngrd
      character*100 fileturb,fileprof,filemean,filewater
      integer hour
      real hour2,hour3,height
      real*8 distY_1,dp,LWP
      
      difloc=.false.
      nonloc_tm86=.false.
      nonloc_ls96=.false.
      nonloc_noh03=.false.
      nonloc_lock=.false.
      loc_inm=.false.  
!  read parameters      
      call readpa(fnmap,fnfor,fnout,itheta
     :   ,iwind,scrout,fngrd,difloc,nonloc_tm86,nonloc_ls96,
     :    nonloc_noh03)
     
!  open files for output     
      open(10,file='./results/exp/timeserie_GABLS.txt')
!      open(11,file='./results/exp/theta_crossection.txt')
!      open(12,file='./results/exp/u_crossection.txt')
!      open(13,file='./results/exp/v_crossection.txt')
!      open(14,file='./results/exp/mom_crossection.txt')
!      open(15,file='./results/exp/hflx_crossection.txt')
      open(16,file='./results/exp/means_series.txt')
      fileturb='./results/exp/turbulence00.txt'
      fileprof='./results/exp/profiles00.txt'
      filewater='./results/exp/microphys00.txt'
      filemean='./results/exp/mean00.txt'
      
      
!  allocate variables   
      call allocvar
      call vertical_grid
      call iniprofs


      
     
      l=1           ! only for the 1st time step, after that l=2 for leap-frog scheme
      dtl=l*dt
      phi=phi/180.*pi
      fcor=2.*omega*sin(phi)
      distY=0.
      
! *** start main time loop  ***
    
      do 9000 nstep=1,ntime
      if( nstep.ne.1) l=2
      dtl=l*dt
      
      call eq_state    ! p from hydrostatics, ro and t from equation of state

!---------------initial profile of specific humidity----------!      
!      if(nstep.eq.1) then
!      if (qif.ne.0) then
!      do iz=1,nz
!      qv(iz,1)=0. !0.4*qsat(t(iz),p(iz,2))
!      qv(iz,2)=qv(iz,1)
!      qv(iz,3)=qv(iz,1)
!      enddo
!      endif
!      endif
      !------------- surface fluxes-------------!
      if (ifhle.eq.1) then
        h=100.
        le=0.
        cdm=0.1
      endif
      if (tfix.eq.1) then
        call surf_layer_t
      endif 
!      write(0,*) 'AFTER SURFACE'
      !-----------------------------------------!      

      call bl_depth(nonloc_tm86,difloc,nonloc_lock)   !diagnostics of the ABL height
!      write(0,*)'AFTER BL_DEPTH'
      !----------turbulent diffusion------------!      
      if (difloc) then 
	    call diffu_local
	elseif (nonloc_tm86) then
	    call diffu_TM86
	elseif (nonloc_ls96) then
	    call diffu_LS96
	elseif (nonloc_noh03) then
	    call diffu_Noh03
        elseif (loc_inm) then
           call diffu_INM
        elseif (nonloc_lock) then
           call diffu_Lock
      endif
      !-----------------------------------------! 
      call entrainment           ! diagnostics of entrainment rate     
      call meanABL               ! calcultaion of mean ABL quantities
      !-----------------------------------------!    
      dy=ablv*dt                 ! distance in north-south direction
      if (dt*nstep.gt.3600.*icetime) then
        distY_1=distY
        distY=distY+abs(ablv*dt/1000.)
      endif  
      !-----------------------------------------!
    !  call baroclinity           ! calculating pressure from hydrostatic + thermal wind
      
      call radiation  
      call moment                ! integration of equations for u and v
      call thermo                ! integration of equation for theta
      
      if (qif.ne.0) then
        call humid                 ! integration of equation for qv
      endif
      
      if (ifwr.ne.0) then
        call microphysics
        call microphys_update
        call eq_state
        call saturation
        call satur_adjust
        do iz=1,nz
          if(qc(iz,3).lt.0) qc(iz,3)=0.
          if(qr(iz,3).lt.0) qr(iz,3)=0.
          if(qci(iz,3).lt.0) qci(iz,3)=0.
          if(qsn(iz,3).lt.0) qsn(iz,3)=0.
        enddo
      endif
      
     
      !-------------ASELYN TIME FILTER-----------------------!
       do iz=1,nz
         u(iz,2)=u(iz,2)+0.1*(u(iz,1)+u(iz,3)-2.*u(iz,2))
         v(iz,2)=v(iz,2)+0.1*(v(iz,1)+v(iz,3)-2.*v(iz,2))
         th(iz,2)=th(iz,2)+0.1*(th(iz,1)+th(iz,3)-2.*th(iz,2))
         qv(iz,2)=qv(iz,2)+0.1*(qv(iz,1)+qv(iz,3)-2.*qv(iz,2))
         qc(iz,2)=qc(iz,2)+0.1*(qc(iz,1)+qc(iz,3)-2.*qc(iz,2))
         qr(iz,2)=qr(iz,2)+0.1*(qr(iz,1)+qr(iz,3)-2.*qr(iz,2))
         qci(iz,2)=qci(iz,2)+0.1*(qci(iz,1)+qci(iz,3)-2.*qci(iz,2))
         qsn(iz,2)=qsn(iz,2)+0.1*(qsn(iz,1)+qsn(iz,3)-2.*qsn(iz,2))
       enddo
       
       
      !------------------------------------------------------!
       dp=(p(1,2)-p(1,1))/dy*1000.
       LWP=0.
       do iz=1,nz
         u(iz,1)=u(iz,2)
         u(iz,2)=u(iz,3)
         v(iz,1)=v(iz,2)
         v(iz,2)=v(iz,3)
         th(iz,1)=th(iz,2)
         th(iz,2)=th(iz,3)
         qv(iz,1)=qv(iz,2)
         qv(iz,2)=qv(iz,3)
         qc(iz,1)=qc(iz,2)
         qc(iz,2)=qc(iz,3)
         qr(iz,1)=qr(iz,2)
         qr(iz,2)=qr(iz,3)
         qci(iz,1)=qci(iz,2)
         qci(iz,2)=qci(iz,3)
         qsn(iz,1)=qsn(iz,2)
         qsn(iz,2)=qsn(iz,3)
         p(iz,1)=p(iz,2)
         mom(iz)=sqrt((def13(iz)+def13c(iz))**2.
     :    +(def23(iz)+def13c(iz))**2.)
         ht(iz)=-h3(iz)-h3c(iz)-h3e(iz)
         LWP=LWP+ro(iz)*qc(iz,2)*1000.*dz(iz)
       enddo
             

       !-------writing YZ cross-sections to output-----------!
!       if (distY_1.eq.0.and.distY.ne.0) then
!         call crossection(distY_1,distY,5.,nz,z_sl,z,th,11)
!         call crossection(distY_1,distY,5.,nz,z_sl,z,u,12)
!         call crossection(distY_1,distY,5.,nz,z_sl,z,v,13)
!         call crossection2(distY_1,distY,5.,nz,z_sl,z,mom,14)
!         call crossection2(distY_1,distY,5.,nz,z_sl,z,ht,15)
!       endif
!       if (distY.gt.0.and.distY.le.350) then
!       if (mod(distY_1,5.).gt.mod(distY,5.)) then 
!         write(0,*) distY-mod(distY,5.)
!         call crossection(distY_1,distY,5.,nz,z_sl,z,th,11)
!         call crossection(distY_1,distY,5.,nz,z_sl,z,u,12)
!         call crossection(distY_1,distY,5.,nz,z_sl,z,v,13)
!         call crossection2(distY_1,distY,5.,nz,z_sl,z,mom,14)
!         call crossection2(distY_1,distY,5.,nz,z_sl,z,ht,15)
!       endif
!       endif
       !------------------------------------------------------!
       
       !-----------writing time series------------------------!
      if(mod(nstep*dt,60.).eq.0) write(10,
     : '(f7.2,f11.3,3f9.2,f6.3,f9.2,f9.4,f10.7,2f10.2,
     :     f9.2,f13.2,f9.3)')
     : nstep*dt/60.,distY,u(1,1),v(1,1),
     : sqrt(u(1,1)**2+v(1,1)**2),ust_s,
     : hbl,-tst_s*ust_s,-qst_s*ust_s,h,le,LWP, th(1,1),qv(1,1)
     
       if(mod(nstep*dt,60.).eq.0) write(16,
     :  '(f7.2,f11.3,f10.2,5f12.6)')
     :  nstep*dt/60.,distY,mth,mqv*1000.,mqc*1000.,mqr*1000.,
     : mqci*1000.,mqsn*1000.
       !------------------------------------------------------!
       
       !-----------writing profiles---------------------------!
       if(mod(nstep*dt,3600.).eq.0)then
       
  !     if(distY.gt.0.and.mod(int(distY-dy),50).ne.0
  !   :          .and.mod(int(distY),50).eq.0)then
        hour=nstep*dt/3600
  
        if(hour.gt.9) write(fileturb(25:26),'(i2)')hour
         if(hour.le.9) write(fileturb(26:26),'(i1)')hour
         if(hour.gt.9) write(fileprof(23:24),'(i2)')hour
         if(hour.le.9) write(fileprof(24:24),'(i1)')hour
         if(hour.gt.9) write(filewater(24:25),'(i2)')hour
         if(hour.le.9) write(filewater(25:25),'(i1)')hour
    !       write(fileturb(38:48),'(i4)')int(distY)
    !      write(fileprof(36:46),'(i4)')int(distY)
    !      write(filewater(37:47),'(i4)')int(distY)
         open(20,file=fileprof)
         open(21,file=fileturb)
         open(22,file=filewater)
         do iz=0,nz
           if(iz.eq.0) height=z(0)
           if(iz.gt.0) height=z(iz)-z_sl
           
              write(20,'(f5.0,4 f10.2,f10.6,3f10.5)')
     :        height, u(iz,3),v(iz,3),
     :        sqrt(u(iz,3)**2.+v(iz,3)**2),th(iz,3),qv(iz,3),
     :        difk(iz),dift(iz),ri(iz) !,sqrt((ug+dpdy(iz)/fcor)**2+vgeos(iz)**2)
     
              write(22,'(f5.0,8f12.8,3f13.9)')
     :        height,qv(iz,3),qsat(t(iz),p(iz,2)),qsati(t(iz),p(iz,2)),
     :        qs(iz),qc(iz,2),qr(iz,2),qci(iz,2),qsn(iz,2)
     :        ,hlat/cp*(p00/p(iz,2))**akapa*condensat(iz)
     :        ,sublim(iz),difunt(iz)-
     :    hlat/cp*difunqv(iz)
         enddo
         do iz=1,nz
           write(21,'(f5.0,4f10.4,f13.8,f10.4)')
     :     0.5*(z(iz-1)+z(iz))-z_sl*0.5,
     :     ht(iz),mom(iz),def13(iz)+def13c(iz),def23(iz)+def23c(iz),
     :     wq3(iz),ht(iz)+0.61*th(iz,3)*wq3(iz)+0.61*qv(iz,3)*ht(iz)
         enddo
         close(20)
         close(21)
         close(22)
      endif
      !-------------------------------------------------------------!
      
      !---------------writing averaged profiles---------------------!
      hour2=nstep*dt/3600.
      hour3=hour
       if(hour2.ge.hour.and.hour2.lt.hour+0.5) then
       do iz=1,nz
       sh3(iz)=sh3(iz)+h3(iz)
       enddo
       endif
       
       if(mod(hour2,hour+0.5).eq.0) then
          if(hour.gt.9) write(filemean(32:33),'(i2)')hour
          if(hour.le.9) write(filemean(33:33),'(i1)')hour
          open(22,file=filemean)
          do iz=1,nz
            write(22,'(f5.0,f10.5)')0.5*(z(iz-1)+z(iz))-z_sl*0.5
     :            ,-sh3(iz)/(1800./dt)
          enddo
          sh3=0.
          close(22)
       endif
       !-----------------------------------------------------------!
           
9000  continue
      
      end program nh1d
      
