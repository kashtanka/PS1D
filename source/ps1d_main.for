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
      
!  read parameters from text file   
      call readpa(fnmap,fnfor,fnout,itheta
     :   ,iwind,scrout,fngrd,difloc,nonloc_tm86,nonloc_ls96,
     :    nonloc_noh03,loc_inm,nonloc_lock)
     
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
!---- allocate variables   
      call allocvar
      call vertical_grid
      call iniprofs         ! initial profs interpolated on model grid
      l=1           ! only for the 1st time step, after that l=2 for leap-frog scheme
      dtl=l*dt
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
      if (tfix.eq.1) then        ! 
        call surf_layer_t        ! Monin-Obukhov similarity
      endif 
      if (ifhle.eq.1) then       ! constant surface fluxes
         call const_fluxes
      endif
      !-----------------------------------------!      

      call bl_depth(nonloc_tm86,difloc,nonloc_lock)   !diagnostics of the ABL height
      !----------turbulent diffusion------------!      
      if (difloc) then 
	    call diffu_local      ! local closure
	elseif (nonloc_tm86) then
	    call diffu_TM86       ! Troen and Mahrt 
	elseif (nonloc_ls96) then
	    call diffu_LS96       ! Lupkes and Schlunzen
	elseif (nonloc_noh03) then
	    call diffu_Noh03      ! Noh et al
        elseif (loc_inm) then
           call diffu_INM         ! local closure from the INM model
        elseif (nonloc_lock) then
           call diffu_Lock        ! Lock
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
    !  call baroclinity         
      if (radif) then
         call radiation   ! rad fluxes and tendencies
      endif  
      if (vadv) then
         call vert_adv    ! vertical advection
      endif
      call moment      ! integration of equations for u and v
      call thermo      ! integration of equation for theta    
      if (qif.ne.0) then
         call humid     ! integration of equation for qv
      endif
    
     
      if (ifwr.ne.0) then          ! if there are clouds
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

!--------Aselyn time filter for the leap-frog scheme---------!
      call aselyn       
!------------------------------------------------------------!
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
!---------write output to files------!
       call output(fileturb,fileprof,filemean,filewater,hour,
     :             hour2,hour3,LWP)             
!------------------------------------!
           
9000  continue
      
      end program nh1d
      
