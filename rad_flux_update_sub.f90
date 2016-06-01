      MODULE WORK_OPTICS

      SAVE

!-----Input arrays into sorad---------------------
      real, allocatable :: pl_rad (:,:)
      real, allocatable :: ta_rad (:,:)
      real, allocatable :: wa_rad (:,:)
      real, allocatable :: zlayer (:,:)
      real, allocatable :: fcld   (:,:)
      real, allocatable :: zdel   (:,:)
      real, allocatable :: cwc    (:,:,:)
!-----Output arrays from sorad---------------------           
      real, allocatable :: flx     (:,:), flc     (:,:)
      real, allocatable :: flx_d   (:,:), flx_u   (:,:)
      real, allocatable :: flc_d   (:,:), flc_u   (:,:)
      real, allocatable :: fdiruv  (:)  , fdifuv  (:)
      real, allocatable :: fdirpar (:)  , fdifpar (:)
      real, allocatable :: fdirir  (:)  , fdifir  (:)
!-----Input arrays to irrad------------------------
      real, allocatable :: fs     (:,:), tg(:,:), tv(:,:)
      real, allocatable :: tsurfs (:,:)
      real, allocatable :: t_atm_surf(:)
!-----Output arrays from irrad---------------------                 
      real, allocatable :: flx_lw(:,:), flc_lw(:,:)
      real, allocatable :: dfdts (:,:), sfcem (:)
      real, allocatable :: cosz  (:)

      real, allocatable :: alb_y(:), emis_y(:)

      contains
      SUBROUTINE ALLOC_WORK_OPTICS(nyi,nye,ns)
      implicit none
      integer(4), intent(in) :: nyi, nye, ns

      allocate (pl_rad (nyi:nye,1:ns) ) 
      allocate (ta_rad (nyi:nye,1:ns-1) )
      allocate (wa_rad (nyi:nye,1:ns-1) )
      allocate (zlayer (nyi:nye,1:ns) )
      allocate (fcld   (nyi:nye,1:ns-1) )
      allocate (zdel   (nyi:nye,1:ns-1) )
      allocate (cwc    (nyi:nye,1:ns-1,1:3) )
      allocate (flx    (nyi:nye,1:ns), flc  (nyi:nye,1:ns) )
      allocate (flx_d  (nyi:nye,1:ns), flx_u(nyi:nye,1:ns) )
      allocate (flc_d  (nyi:nye,1:ns), flc_u(nyi:nye,1:ns) )
      allocate (fdiruv (nyi:nye), fdifuv  (nyi:nye) )
      allocate (fdirpar(nyi:nye), fdifpar (nyi:nye) )
      allocate (fdirir (nyi:nye), fdifir  (nyi:nye) )
      allocate (fs     (nyi:nye,1:2), tg(nyi:nye,1:2), tv(nyi:nye,1:2) )
      allocate (tsurfs (nyi:nye,1:2) )
      allocate (t_atm_surf(nyi:nye) )
      allocate (flx_lw (nyi:nye,1:ns), flc_lw(nyi:nye,1:ns) )
      allocate (dfdts  (nyi:nye,1:ns), sfcem(nyi:nye) )
      allocate (cosz   (nyi:nye) )
      allocate (alb_y  (nyi:nye), emis_y(nyi:nye))
      END SUBROUTINE ALLOC_WORK_OPTICS


      SUBROUTINE DEALLOC_WORK_OPTICS
      implicit none
      deallocate (pl_rad) 
      deallocate (ta_rad)
      deallocate (wa_rad)
      deallocate (zlayer)
      deallocate (fcld)
      deallocate (zdel)
      deallocate (cwc)
      deallocate (flx, flc)
      deallocate (flx_d, flx_u)
      deallocate (flc_d, flc_u)
      deallocate (fdiruv, fdifuv)
      deallocate (fdirpar, fdifpar)
      deallocate (fdirir, fdifir)
      deallocate (fs, tg, tv)
      deallocate (tsurfs)
      deallocate (t_atm_surf)
      deallocate (flx_lw, flc_lw)
      deallocate (dfdts, sfcem)
      deallocate (cosz)
      deallocate (alb_y, emis_y)
      END SUBROUTINE DEALLOC_WORK_OPTICS


      END MODULE WORK_OPTICS


      SUBROUTINE RAD_FLUX_UPDATE &
      & (nx, ny, ns, &
      & imonth,iday,ihour,iminu,iseco,ntime,xlatit, &
      & ifqv, ifqc, ifqr, shortwave, longwave, &
      & sigma1,sigma0,ds1,pp,ptop,phis,phi, &
      & pts,pt,qv,qc,qr,qci,qsn, &
      & hsuf,tsurf,alb,emis, &
      & Srad,Srad_surf,Lrad,Lrad_surf, &
      & SraddirUV_surf,SraddirPAR_surf,SraddirIR_surf, &
      & SraddifUV_surf,SraddifPAR_surf,SraddifIR_surf)
      
!     Subroutine RAD_FLUX_UPDATE updates
!     the shortwave and longwave radiation fluxes

      use OPTIC_PARAMETERS
      use ALLOC, only: &
      & g, &
      & p00, &
      & akapa
      use MPI_VARIABLES, only : &
      & nxs2i,nxs2e,nys2i,nys2e, &
      & nxi,nxe,nyi,nye, &
      & nx0i,nx1e,ny0i,ny1e, &
      & rank_comm3d, comm3d, coords 
!      & MPI_BARRIER
      use WORK_OPTICS

      implicit none

      include 'mpif.h'


!     Input variables
      integer(4), intent(in) :: nx, ny, ns
      integer(4), intent(in) :: imonth, iday, ihour, iminu, iseco, ntime 
      integer(4), intent(in) :: ifqv, ifqc, ifqr
      integer(4), intent(in) :: shortwave, longwave
      
      real(8), intent(in) :: xlatit
      
      real(8), intent(in) :: sigma1 (0:ns+1)
      real(8), intent(in) :: sigma0 (0:ns+1)
      real(8), intent(in) :: ds1    (0:ns+1)     
      real(8), intent(in) :: pp     (nxs2i:nxs2e,nys2i:nys2e,1:4)
      real(8), intent(in) :: ptop
      real(8), intent(in) :: phis   (nxs2i:nxs2e,nys2i:nys2e,0:ns+1)
      real(8), intent(in) :: phi    (nxs2i:nxs2e,nys2i:nys2e,0:ns+1)
      real(8), intent(in) :: pts    (nxs2i:nxs2e,nys2i:nys2e,0:ns+1)
      real(8), intent(in) :: pt     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: qv     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: qc     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: qr     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: qci    (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: qsn    (nxs2i:nxs2e,nys2i:nys2e,0:ns+1,1:3)
      real(8), intent(in) :: hsuf   (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(in) :: tsurf  (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(in) :: alb    (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(in) :: emis   (nxs2i:nxs2e,nys2i:nys2e)
      
!     Output variables      
      real(8), intent(out) :: Srad     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1) 
      real(8), intent(out) :: Lrad     (nxs2i:nxs2e,nys2i:nys2e,0:ns+1) 
      real(8), intent(out) :: Srad_surf(nxs2i:nxs2e,nys2i:nys2e) 
      real(8), intent(out) :: Lrad_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddirUV_surf (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddirPAR_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddirIR_surf (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddifUV_surf (nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddifPAR_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(out) :: SraddifIR_surf (nxs2i:nxs2e,nys2i:nys2e)

!     Local variables

      real(8), external :: SIN_SUN
            
      real :: t1_rad, t2_rad 
     
      integer(4) :: ict, icb
      integer(4) :: ix, iy, is ! Loop indices
      integer(4) :: nn
      
!     Switches for radiation parameterizations      
      logical, save :: cldwater = .true.
      logical, save :: overcast = .true.
      logical, save :: overcast_lw = .true.
      logical, save :: cldwater_lw = .true.
      logical, save :: high = .false.
      logical, save :: trace = .true.
      logical, save :: vege = .true.
      logical, save :: aerosol = .true.
      logical, save :: firstcall = .true.
             
!     The body of subroutine
      if (firstcall) then
        call ALLOC_OPTICS(nyi,nye,ns-1)
      endif

      call ALLOC_WORK_OPTICS(nyi,nye,ns)

!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d == 0) write(*,*) 'The 2 barrier passed'


      ix_cycle: do ix = nxi, nxe !1, nx 

        cosz(nyi:nye) = sngl(SIN_SUN(imonth,iday,ihour,iminu,iseco,xlatit))

        do is = 1, ns
          do iy = nyi, nye ! 1, ny 
!           The pressure pl_rad is expressed in hPa          
            pl_rad(iy,is) = sngl( 1.d-2*(sigma1(is)*pp(ix,iy,2)+ptop) )
!           The numbers of sigma-levels corresponding to 400 hPa and 700 hPa
!           are evaluated in the middle in meridional section of the domain (iy=ny/2)

!_commented for debug            if (is >= 2 .and. iy == int(0.5*(float(nyi)+float(nye)))) then 
!_commented for debug             if (pl_rad(iy,is)>400. .and. pl_rad(iy,is-1)<400.) ict = is 
!_commented for debug             if (pl_rad(iy,is)>700. .and. pl_rad(iy,is-1)<700.) icb = is 
!_commented for debug           endif

            ict = 5 ! Debugging
            icb = 10 ! debugging

!            zlayer(iy,is) = (phis(ix,iy,is)+phi(ix,iy,is))/g
            zlayer(iy,is) = sngl ( ( phis(ix,iy,is-1) + phi(ix,iy,is-1) + &
            & ds1(is-1)/(ds1(is-1) + ds1(is)) * &
            & (phis(ix,iy,is) + phi(ix,iy,is) - &
            & phis(ix,iy,is-1) - phi(ix,iy,is-1)) ) / g )

            if (is == ns) zlayer(iy,is) = sngl (hsuf(ix,iy))
!           Converting zlayer to km, as needed to input in SORAD
            zlayer(iy,is) = sngl ( zlayer(iy,is)*1.d-3 )
          enddo  
        enddo  

        do is = 1, ns-1 ! 1, ns-1  
          do iy = nyi, nye  ! 1, ny
!            t1_rad = (pts(ix,iy,is)+pt(ix,iy,is,2))/ & 
!            & (p00/(sigma1(is)*pp(ix,iy,2)+ptop))**akapa  
!            t2_rad = (pts(ix,iy,is+1)+pt(ix,iy,is+1,2))/ &
!            & (p00/(sigma1(is+1)*pp(ix,iy,2)+ptop))**akapa 
!            ta_rad(iy,is) = 0.5d0*(t1_rad + t2_rad)  
            ta_rad(iy,is) = sngl ( (pts(ix,iy,is)+pt(ix,iy,is,2))/ & 
            & (p00/(sigma0(is)*pp(ix,iy,2)+ptop))**akapa )
            wa_rad(iy,is) = sngl (qv(ix,iy,is,2))
            zdel(iy,is)   = zlayer(iy,is) - zlayer(iy,is+1) 
!--------Notes:
!--------fcld is not currently calculated in atmospheric model---------------
!        fcld is cloud amount (fraction)
            fcld(iy,is)   = 0.
            cwc(iy,is,1)  = qci(ix,iy,is,2) ! Ice crystals, snow is not included
            cwc(iy,is,2)  = qc(ix,iy,is,2) ! Cloud droplets
            cwc(iy,is,3)  = qr(ix,iy,is,2) ! Rain droplets
            if (ifqc == 1) & 
            & cwc(iy,is,2) = sngl ( qc(ix,iy,is,2) )
            if (ifqr == 1) &
            & cwc(iy,is,3) = sngl ( qr(ix,iy,is,2) )
          enddo  
        enddo 

        do iy = nyi, nye !1, ny
          tsurfs(iy,1:2) = sngl ( tsurf(ix,iy) )
          t_atm_surf(iy) = sngl ( 0.5d0*(tsurfs(iy,1) + ta_rad(iy,ns-1)) )
          alb_y(iy) = sngl (alb(ix,iy))
          emis_y(iy) = sngl (emis(ix,iy))
        enddo

        call OPTICPARSET(nyi,nye,nyi,nye,ns-1,alb_y,emis_y,zlayer)

! start debugging 

!        if (ix == 20 .and. coords(1) == 0 .and. coords(2) == 0) then
!            write(*,'(a10, 2i)') '1var', nye, ns
!            print*,'2var', crel
!            write(*,'(a10, 2f20.15)') '2.5var',zlayer(1:2,20)
!            write(*,'(a10, 2f20.15)') '3var', pl_rad(1:2,20)
!            write(*,'(a10, 2f20.15)') '4var', ta_rad(1:2,20)
!            write(*,'(a10, 2f20.15)') '5var', wa_rad(1:2,20) 
!            write(*,'(a10, 2f20.15)') '6var', oa(1:2,20)
!            write(*,'(a10, f20.15)') '7var', co2_sorad 
!            write(*,'(a10, 2f20.15)') '8var', zdel(1:2,20) 
!            print*,'9var', overcast
!            print*,'10var', cldwater 
!            write(*,'(a10, 2f20.15)') '11var', cwc(1:2,20,1)
!            write(*,'(a10, 2f20.15)') '12var', taucld(1:2,20,1)
!            write(*,'(a10, 2f20.15)') '13var', reff(1:2,20,1)
!            write(*,'(a10, 2f20.15)') '14var', fcld(1:2,20) 
!            write(*,'(a10, 2i)') '15var', ict, icb 
!            write(*,'(a10, 2f20.15)') '16var', taual(1:2,20,1)
!            write(*,'(a10, 2f20.15)') '17var', ssaal(1:2,20,1)
!            write(*,'(a10, 2f20.15)') '18var', asyal(1:2,20,1) 
!            write(*,'(a10, 2f20.15)') '19var', cosz(1:2)  
!            write(*,'(a10, 2f20.15)') '20var', rsuvbm(1:2)
!            write(*,'(a10, 2f20.15)') '21var', rsuvdf(1:2)
!            write(*,'(a10, 2f20.15)') '22var', rsirbm(1:2)
!            write(*,'(a10, 2f20.15)') '23var', rsirdf(1:2)  
!            write(*,'(a10, 2f20.15)') '24var', flx(1:2,20)
!            write(*,'(a10, 2f20.15)') '25var', flc(1:2,20) 
!            write(*,'(a10, 2f20.15)') '26var', fdiruv(1:2)
!            write(*,'(a10, 2f20.15)') '27var', fdifuv(1:2)
!            write(*,'(a10, 2f20.15)') '28var', fdirpar(1:2) 
!            write(*,'(a10, 2f20.15)') '29var', fdifpar(1:2)
!            write(*,'(a10, 2f20.15)') '30var', fdirir(1:2)
!            write(*,'(a10, 2f20.15)') '31var', fdifir(1:2)   
!            write(*,'(a10, 2f20.15)') '32var', flx_d(1:2,20)
!            write(*,'(a10, 2f20.15)') '33var', flx_u(1:2,20)
!            write(*,'(a10, 2f20.15)') '34var', flc_d(1:2,20)
!            write(*,'(a10, 2f20.15)') '35var', flc_u(1:2,20)
!            write(*,'(a10, 2f20.15)') '36var', cosz(1:2)
!            write(*,'(a10, 2f20.15)') '42var', tsurfs(1:2,1)
!            write(*,'(a10, 2f20.15)') '43var', alb_y(1:2)
!            write(*,'(a10, 2f20.15)') '44var', emis_y(1:2)
!            STOP
!        endif
        

!        oa = 3.E-04
!        taucld = 0.001
!        reff = 0.01
!        taual = 0.001
!        ssaal = 0.0001
!        asyal = 0.0001
!        rsuvbm = 0.5
!        rsuvdf = 0.5
!        rsirbm = 0.5
!        rsirdf = 0.5

! end debugging

!        call MPI_BARRIER(comm3d,nn)
!        if (rank_comm3d==0) write(*,*) 'The 3 barrier passed', ix

        if (cosz(nyi) > 0) then   
          call SORAD (nye-nyi+1,ns-1,&
          & crel, &
          & pl_rad(nyi,1),ta_rad(nyi,1),wa_rad(nyi,1), &
          & oa(nyi,1),co2_sorad, &
          & zdel(nyi,1), &
          & overcast,cldwater, &
          & cwc(nyi,1,1),taucld(nyi,1,1),reff(nyi,1,1),fcld(nyi,1), &
          & ict,icb, &
          & taual(nyi,1,1),ssaal(nyi,1,1),asyal(nyi,1,1), &
          & cosz(nyi), & 
          & rsuvbm(nyi),rsuvdf(nyi),rsirbm(nyi),rsirdf(nyi), & 
          & flx(nyi,1),flc(nyi,1), &
          & fdiruv(nyi),fdifuv(nyi),fdirpar(nyi), &
          & fdifpar(nyi),fdirir(nyi),fdifir(nyi), &  
          & flx_d(nyi,1),flx_u(nyi,1),flc_d(nyi,1),flc_u(nyi,1))

!          call MPI_BARRIER(comm3d,nn)
!          if (rank_comm3d==0) write(*,*) 'The 4 barrier passed', ix

          if (shortwave == 1) then 
            do is = 1, ns 
              do iy = nyi, nye
                Srad(ix,iy,is) = s0*flx(iy,is)*amax1(cosz(iy),0.)
              enddo
            enddo  
          else 
            Srad(ix,nyi:nye,1:ns) = 0. 
          endif  

!         Srad_surf is the downward solar radiation at the horizontal earth surface
!         (positive downward)
          do iy = nyi, nye !1, ny
            Srad_surf(ix,iy) = &
            & flx(iy,ns)*s0*amax1(cosz(iy),0.)/(1.d0-alb(ix,iy))
            SraddirUV_surf(ix,iy) = &
            & fdiruv(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsuvbm(iy))
            SraddirPAR_surf(ix,iy) = &
            & fdirpar(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsuvbm(iy))
            SraddirIR_surf(ix,iy) = &
            & fdirir(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsirbm(iy))
            SraddifUV_surf(ix,iy) = &
            & fdifuv(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsuvdf(iy))
            SraddifPAR_surf(ix,iy) = &
            & fdifpar(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsuvdf(iy))
            SraddifIR_surf(ix,iy) = &
            & fdifir(iy)*s0*amax1(cosz(iy),0.)/(1.d0-rsirdf(iy))
          enddo
        else  
          Srad(ix,nyi:nye,1:ns) = 0.
          Srad_surf(ix,nyi:nye) = 0.
          SraddirUV_surf(ix,nyi:nye) = 0.
          SraddirPAR_surf(ix,nyi:nye) = 0.
          SraddirIR_surf(ix,nyi:nye) = 0.
          SraddifUV_surf(ix,nyi:nye) = 0.
          SraddifPAR_surf(ix,nyi:nye) = 0.
          SraddifIR_surf(ix,nyi:nye) = 0.
        endif

!        if (ntime == 720*3 .and. ix == (nx+1)/2 .and. &
!           & nyi <= (ny+1)/2 .and. nye >= (ny+1)/2) then
!           write(*,*) 'Radiation ', 'rank ', rank_comm3d, &
!           & 'm,d,h,m,s', imonth, iday, ihour, iminu, iseco, &
!           & 'Srad', Srad_surf(ix,(ny+1)/2), &
!           & 'cosz', cosz((ny+1)/2), 'pres', pl_rad((ny+1)/2,:), &
!           & 'temp', ta_rad((ny+1)/2,:), 'hum', wa_rad((ny+1)/2,:)
!           STOP
!        endif

! Debug output of the array
!      write(*,*) 'xycoords = ', nxi ,nxe, nyi, nye
!      call ARRGATHERWRI_MPI(Srad(nxs2i,nys2i,0), & !pt(nxs2i,nys2i,0,3), &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns+1, &
!      & 0,nx+1,0,ny+1,0,ns+1,nx+1,ny+1,ns+1)
!      STOP



        fs(nyi:nye,1) = 0.5  
        fs(nyi:nye,2) = 0.5 
!         print*, 'nsur, na', nsur, na
!        call MPI_BARRIER(comm3d,nn)
!        if (rank_comm3d==0) write(*,*) 'The 5 barrier passed', ix

        call IRRAD(nye-nyi+1,ns-1,pl_rad,ta_rad,wa_rad,oa, & 
        &          t_atm_surf, &
        &          co2,high,trace,n2o,ch4,cfc11,cfc12,cfc22, &
        &          vege,nsur,fs,tsurfs, &
        &          eg,tsurfs,ev,rvir, &
        &          overcast_lw,cldwater_lw,cwc,taucl_lw,fcld,ict,icb, &
        &          aerosol,na,taual_lw,ssaal_lw,asyal_lw, &
        &          flx_lw,flc_lw,dfdts,sfcem)

!        call MPI_BARRIER(comm3d,nn)
!        if (rank_comm3d==0) write(*,*) 'The 6 barrier passed', ix

!       Lrad_surf is the atmospheric radiation at the earth surface
!       (positive downward)
        do iy = nyi, nye !1, ny
          Lrad_surf(ix,iy) = flx_lw(iy,ns) - sfcem(iy)
        enddo
        if (longwave==1) then
          do is = 1, ns
            do iy = nyi, nye !1, ny
              Lrad(ix,iy,is) = flx_lw(iy,is)
            enddo
          enddo
        else  
          Lrad(ix,nyi:nye,1:ns) = 0.  
        endif 
      enddo ix_cycle 

!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 7 barrier passed'
     
!      deallocate (pl_rad)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-1 barrier passed'
!      deallocate (ta_rad)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-2 barrier passed'
!      deallocate (wa_rad)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-3 barrier passed'
!      deallocate (zlayer)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-4 barrier passed'
!      deallocate (fcld)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-5 barrier passed'
!      deallocate (zdel)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-6 barrier passed'
!      deallocate (cwc)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-7 barrier passed'
!      deallocate (flx)
!      deallocate (flc)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-8 barrier passed'
!      deallocate (flx_d)
!      deallocate (flx_u)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-9 barrier passed'
!      deallocate (flc_d)
!      deallocate (flc_u)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-10 barrier passed'
!      deallocate (fdiruv)
!      deallocate (fdifuv)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-11 barrier passed'
!      deallocate (fdirpar)
!      deallocate (fdifpar)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-12 barrier passed'
!      deallocate (fdirir)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-12-1 barrier passed'
!      deallocate (fdifir)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-13 barrier passed'
!      deallocate (fs) 
!      deallocate (tg)
!      deallocate (tv)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-14 barrier passed'
!      deallocate (tsurfs)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-15 barrier passed'
!      deallocate (t_atm_surf)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-16 barrier passed'
!      deallocate (flx_lw) 
!      deallocate (flc_lw)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-17 barrier passed'
!      deallocate (dfdts)
!      deallocate (sfcem)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-18 barrier passed'
!      deallocate (cosz)
!      call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-19 barrier passed'
!      deallocate (alb_y)
!      deallocate (emis_y)
!       call MPI_BARRIER(comm3d,nn)
!      if (rank_comm3d==0) write(*,*) 'The 8-20 barrier passed'  
     

      call DEALLOC_WORK_OPTICS

      if (firstcall) firstcall = .false.
      END SUBROUTINE RAD_FLUX_UPDATE
