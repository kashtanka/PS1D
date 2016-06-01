#include "debug.inc"

      SUBROUTINE THERMO_MPI &
     & (timeho,ptm2,wrk, &
     & nstep,nx,ny,ns,nx1,ny1,ns1, &
     & nx2,ny2,nx1qv,ny1qv,ns1qv, &
     & nx1qc,ny1qc,ns1qc,nx1qr,ny1qr,ns1qr, &
     & s,dtxs,dtys,dxys,ds0,ds08a,ds02,ds12,ds1, &
     & sigma0,sigma1,taudra,hdamp,ugeos,vgeos,tauspo,tauspoy, &
     & dtxy,dxdt,dydt,dt,dtl,dx,dy,dx2,dy2,fcor,deltaz, &
     & u,v,w,wsig,pt,pts,dptsdtsig,pp,pp10,pp01,ptop, &
     & phis,phi,difunt,difuntcg,dift001, &
     & cond,evap,vini,hmfrz,sacrw,iacr,sbercw,ibercw,sacrr,smlt,imlt, &
     & qv,qc,qr,qci,qsn,rhoasurf,cdcoef,filt,vdampt, &
     & ptbx,ptby,ptcc,ptbxe,ptbye, &
     & ubxe, ubye, vbxe, vbye, & 
     & Radheat,Srad,Lrad,Srad_surf,Lrad_surf, &
     & SraddirUV_surf,SraddirPAR_surf,SraddirIR_surf, &
     & SraddifUV_surf,SraddifPAR_surf,SraddifIR_surf, &
     & hk1, &
     & xlatit,uvsuf,hsuf,h,tsurf,alb,emis,xlake,tslake, &
     & isqv,ifqc,ifqr,ifqi, &
     & radpar,shortwave,longwave,nradcall,ifsoil,ifhle, &
     & iobptx,iobpty,nchn1,idrmax,numsmo, &
     & nxsponge1,nxsponge2,nysponge1,nysponge2, &
     & ntypesponge, & 
     & imonth,iday,ihour,iminu,iseco, &
     & tfct,prt,ptbout,raylei,dohsmo,dovsmt,impldiff,tsdif,refsadjust)
     
!-----------------------------------------------------------------------
!     time integration of potential temperature equation
!-----------------------------------------------------------------------
!

      use ALLOC, only: &
      & hlatcp, hsubcp, hfuscp, &
      & p00, &
      & akapa, &
      & cp, &
      & r, g

      use MPI_VARIABLES, only : &
      & nxs2i,nxs2e,nys2i,nys2e, &
      & nxi,nxe,nyi,nye, &
      & nx0i,nx1e,nx2e,ny0i,ny1e,ny2e, &
      & nxsh,nysh,nxsh2,nysh2,nssh, &
      & comm3d,size_MPI,ndim,coords,dims, &
      & prplcxy,prplc, &
      & parallel, &
      & req, request_send_2d, request_recv_2d, &
      & request_send_3d, request_recv_3d, &
      & isperiodic, &
      & xlbound, xrbound, ylbound, yrbound

      use ARRNDEXCH, only : &
      & ARR3DEXCH, & 
      & ARR2DEXCH

      use TIMING_MOD, only : &
      & pr_time, proc_time

      use MODCONTR, only : &
      & lbcs

      implicit none

      include 'mpif.h'
      
!     Input/output variables      

      real(8), intent(inout) :: ptm2(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: wrk(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: pt(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(inout) :: pp(nxs2i:nxs2e,nys2i:nys2e,1:4)
      real(8), intent(inout) :: s(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: pts(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: dptsdtsig(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: u(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(in) :: v(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(in) :: w(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(inout) :: wsig(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(inout) :: phis(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: phi(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: difunt(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: difuntcg(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: dift001(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: cond(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: evap(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: vini(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: hmfrz(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: sacrw(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: iacr(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: sbercw(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: ibercw(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: sacrr(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: smlt(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: imlt(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: pp10(nxs2i:nxs2e,nys2i:nys2e,1:4)
      real(8), intent(inout) :: pp01(nxs2i:nxs2e,nys2i:nys2e,1:4)
      real(8), intent(inout) :: qv(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(inout) :: qc(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(inout) :: qr(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(in) :: qci(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
      real(8), intent(in) :: qsn(nxs2i:nxs2e,nys2i:nys2e,0:ns1,1:3)
   
      real(8), intent(inout) :: ptbx(1:2,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: ptby(nxs2i:nxs2e,1:2,0:ns1)
      real(8), intent(inout) :: ptcc(1:2,1:2,0:ns1)
      real(8), intent(in) :: ptbxe(1:2,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: ptbye(nxs2i:nxs2e,1:2,0:ns1)
      real(8), intent(in) :: ubxe(1:2,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: ubye(nxs2i:nxs2e,1:2,0:ns1)
      real(8), intent(in) :: vbxe(1:2,nys2i:nys2e,0:ns1)
      real(8), intent(in) :: vbye(nxs2i:nxs2e,1:2,0:ns1)
      
      real(8), intent(inout) :: uvsuf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: hsuf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: h(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: tsurf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: alb(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: emis(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: xlake(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: tslake(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: deltaz(nxs2i:nxs2e,nys2i:nys2e)
      
      real(8), intent(inout) :: Radheat(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: Srad(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: Lrad(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
      real(8), intent(inout) :: Srad_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: Lrad_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddirUV_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddirPAR_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddirIR_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddifUV_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddifPAR_surf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: SraddifIR_surf(nxs2i:nxs2e,nys2i:nys2e)
      
      real(8), intent(inout) :: hk1(nxs2i:nxs2e,nys2i:nys2e)
      
      real(8), intent(inout) :: dtxs(0:ns1)
      real(8), intent(inout) :: dtys(0:ns1)
      real(8), intent(inout) :: dxys(0:ns1)
      real(8), intent(in) :: ds0(0:ns1)
      real(8), intent(in) :: ds02(0:ns1)
      real(8), intent(in) :: ds08a(0:ns1)
      real(8), intent(in) :: ds12(0:ns1)
      real(8), intent(in) :: ds1(0:ns1)
      real(8), intent(in) :: sigma0(0:ns1)
      real(8), intent(in) :: sigma1(0:ns1)
      real(8), intent(inout) :: taudra(0:ns)
      real(8), intent(inout) :: hdamp(0:ns1)
      real(8), intent(in) :: ugeos(0:ns), vgeos(0:ns)
      real(8), intent(in) :: tauspo(nxs2i:nxs2e)
      real(8), intent(in) :: tauspoy(nys2i:nys2e)
      
      real(8), intent(in) :: dtxy
      real(8), intent(in) :: dxdt
      real(8), intent(in) :: dydt
      real(8), intent(in) :: dt
      real(8), intent(in) :: dtl
            
      real(8), intent(in) :: timeho
      real(8), intent(in) :: ptop
      real(8), intent(in) :: dx, dx2
      real(8), intent(in) :: dy, dy2
      
      real(8), intent(in) :: fcor

      real(8), intent(inout) :: filt
      real(8), intent(inout) :: vdampt
      real(8), intent(inout) :: xlatit
      real(8), intent(in) :: rhoasurf(nxs2i:nxs2e,nys2i:nys2e)
      real(8), intent(inout) :: cdcoef
      real(8), intent(in) :: tsdif
      real(8), intent(in) :: refsadjust
      
      integer(4), intent(in) :: nstep            
      integer(4), intent(in) :: ns, ns1
      integer(4), intent(in) :: nx, nx1, nx2
      integer(4), intent(in) :: ny, ny1, ny2
      integer(4), intent(in) :: nx1qv
      integer(4), intent(in) :: ny1qv
      integer(4), intent(in) :: ns1qv
      integer(4), intent(in) :: nx1qc
      integer(4), intent(in) :: ny1qc
      integer(4), intent(in) :: ns1qc
      integer(4), intent(in) :: nx1qr
      integer(4), intent(in) :: ny1qr
      integer(4), intent(in) :: ns1qr
      
      integer(4), intent(in) :: imonth
      integer(4), intent(in) :: iday
      integer(4), intent(in) :: ihour
      integer(4), intent(in) :: iminu
      integer(4), intent(in) :: iseco
            
      integer(4), intent(in) :: isqv, ifqc, ifqr, ifqi
      integer(4), intent(in) :: radpar
      integer(4), intent(in) :: shortwave
      integer(4), intent(in) :: longwave
      integer(4), intent(inout) :: nradcall
      integer(4), intent(in) :: ifsoil
      integer(4), intent(inout) :: ifhle
      integer(4), intent(in) :: iobptx
      integer(4), intent(in) :: iobpty
      integer(4), intent(in) :: nchn1
      integer(4), intent(inout) :: idrmax
      integer(4), intent(inout) :: numsmo
      integer(4), intent(in) :: nxsponge1, nxsponge2
      integer(4), intent(in) :: nysponge1, nysponge2
      integer(4), intent(in) :: ntypesponge
      
      logical, intent(inout) :: tfct
      logical, intent(inout) :: prt
      logical, intent(inout) :: ptbout
      logical, intent(inout) :: raylei
      logical, intent(inout) :: dohsmo
      logical, intent(inout) :: dovsmt
      logical, intent(in) :: impldiff
            
!     Local variables      
      real(8), allocatable :: work1(:,:,:)
      real(8), allocatable :: work2(:,:,:)
      real(8), allocatable :: work3(:,:,:)
      real(8), allocatable :: work4(:,:,:)
      real(8), allocatable :: work5(:,:,:)
      real(8), allocatable :: work6(:,:,:)
      real(8), allocatable :: work7(:,:,:)
      real(8), allocatable :: work8(:,:,:)

      real(8), allocatable :: difs001(:,:,:)
      real(8), allocatable :: st(:,:,:)
      real(8), allocatable :: acoef(:,:,:)
      real(8), allocatable :: bcoef(:,:,:)
      real(8), allocatable :: ccoef(:,:,:)
      real(8), allocatable :: dcoef(:,:,:)
      real(8), allocatable :: work_arrsectget(:,:,:)

      real(8) :: p
      real(8) :: q
      real(8) :: constt
      real(8), allocatable :: tavx(:,:,:), tavy(:,:,:), tavs(:,:,:)
      real(8) :: condens, transform
      real(8) :: temperature
      real(8) :: ctbulk
      real(8) :: tsfunc
       
      integer(4) :: i, j, k, indunst
      integer(4) :: ix, iy, is
      integer(4), parameter :: nodata = -999
      integer(4), parameter :: neighborhood = 2

      character(len=20) :: work_char
       
      logical :: firstcall = .true.

      logical, parameter :: sim = .true. 
      logical, parameter :: nao = .false.
      logical, parameter :: debug_output = .true.

!     External functions
      real(8), external :: TSOIL
      real(8), external :: TLAKE

      if (impldiff) then
        allocate (difs001(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate (st(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate (acoef(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate (bcoef(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate (ccoef(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate (dcoef(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
      endif

      allocate (tavx(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
      allocate (tavy(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
      allocate (tavs(nxs2i:nxs2e,nys2i:nys2e,0:ns1))

!     MPI-exchanges section

      if (parallel) then
!       3-D arrays exchanges
!        pr_time = MPI_WTIME() ! Timing

!        prplc(1) = 1
!        call ARR3DEXCH(u(nxs2i,nys2i,0,2), &
!        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        &  nx0i,nx2e,ny0i,ny2e,0,ns1, &
!        &  nxsh,nysh,0, &
!        &  size_MPI,comm3d,req,ndim,prplc)
!        prplc(3) = 1
!        call ARR3DEXCH(v(nxs2i,nys2i,0,2), &
!        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        &  nx0i,nx2e,ny0i,ny2e,0,ns1, &
!        &  nxsh,nysh,0, &
!        &  size_MPI,comm3d,req,ndim,prplc)
        prplc(1:4) = 1
        call ARR3DEXCH(pt(nxs2i,nys2i,0,2), &
        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
        &  nx0i,nx2e,ny0i,ny2e,0,ns1, &
        &  nxsh,nysh,0, &
        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
        &  size_MPI,comm3d,req,ndim,prplc, &
        &  request_send_3d,request_recv_3d)

!     2-D arrays exchanges
        if (nstep == 1) then
          prplcxy(1) = 1
          call ARR2DEXCH(pp10(nxs2i,nys2i,2), &
          &  nxs2i,nxs2e,nys2i,nys2e, &
          &  nx0i,nx2e,ny0i,ny2e, &
          &  nxsh,nysh, &
          &  comm3d,size_MPI,ndim,prplcxy, &
          &  request_send_2d, request_recv_2d)
          prplcxy(3) = 1
          call ARR2DEXCH(pp01(nxs2i,nys2i,2), &
          &  nxs2i,nxs2e,nys2i,nys2e, &
          &  nx0i,nx2e,ny0i,ny2e, &
          &  nxsh,nysh, &
          &  comm3d,size_MPI,ndim,prplcxy, &
          &  request_send_2d, request_recv_2d)
        endif

!        proc_time(1,2) = proc_time(1,2) + MPI_WTIME() - pr_time ! Timing
      endif




! NOTE: in current MPI-version the fct3d algorithm is not implemented
!       i.e. tfct must be set to .false.

! Start of an unparallelized section TFCT
      if (tfct) then
        call OPTWARNSTOP('tfct==.true.')
        do k = 2, 3
          do is = 0, ns
            do iy = 1, ny1
              do ix = 1, nx1
                pt(ix,iy,is,k)=pt(ix,iy,is,k)*pp(ix,iy,k)
              enddo
            enddo
          enddo
        enddo
        do is = 0, ns
          do iy = 1, ny1
            do ix = 1, nx1
              pt(ix,iy,is,1)=pt(ix,iy,is,1)*pp(ix,iy,1)
            enddo
          enddo
        enddo
        allocate(work1(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work2(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work3(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work4(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work5(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work6(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work7(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
        allocate(work8(nxs2i:nxs2e,nys2i:nys2e,0:ns1))
!         call fct3d(pt(0,0,0,3),pt(0,0,0,2),ptm2, &
!         & u(0,0,0,3),u(0,0,0,2), &
!         & v(0,0,0,3),v(0,0,0,2), &
!         & wsig(0,0,0,3),wsig(0,0,0,2), &
!         & nx,ny,ns, &
!         & dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0, &
!         & work1,work2,work3,work4,work5,work6,work7,work8)
        deallocate(work1)
        deallocate(work2)
        deallocate(work3)
        deallocate(work4)
        deallocate(work5)
        deallocate(work6)
        deallocate(work7)
        deallocate(work8)

        do k = 2, 3
          do is = 0, ns
            do iy = 1, ny1
              do ix = 1, nx1
                pt(ix,iy,is,k)=pt(ix,iy,is,k)/pp(ix,iy,k)
              enddo
            enddo
          enddo
        enddo
        do is = 0, ns
          do iy = 1, ny1
            do ix = 1, nx1
              pt(ix,iy,is,1)=pt(ix,iy,is,1)/pp(ix,iy,1)
            enddo
          enddo
        enddo
        do is = 1, ns-1
          do iy = 2, ny
            do ix = 2, nx
              q=(s(ix,iy,is)+s(ix,iy,is+1))*(pts(ix,iy,is+1) - &
              &  pts(ix,iy,is-1))*(w(ix,iy,is+1,2)+w(ix,iy,is,2))/ds08a(is)
              p=sigma0(is)*pp(ix,iy,2)+ptop
              constt=hlatcp*(p00/p)**akapa
              if(ifqc.ne.0) then
                pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl * &
                & (pp(ix,iy,2)*(q+difunt(ix,iy,is) + &
                &  constt*(cond(ix,iy,is)-evap(ix,iy,is))))/pp(ix,iy,3)
              else
                pt(ix,iy,is,3)=pt(ix,iy,is,3) + dtl * &
                & (pp(ix,iy,2)*(q+difunt(ix,iy,is)))/pp(ix,iy,3)
              endif
            enddo
          enddo
        enddo
! End of an unparallelized section TFCT


      else


        if (impldiff) then

          do is = 1, ns
            do iy = max(nyi,2), nye
              do ix = max(nxi,2), nxe
                difs001(ix,iy,is) = s(ix,iy,is) * dift001(ix,iy,is) / ds0(is)
              enddo
            enddo
          enddo
          do is = 1, ns-1
            do iy = max(nyi,2), nye
              do ix = max(nxi,2), nxe
                st(ix,iy,is) = ( s(ix,iy,is+1) + s(ix,iy,is) ) / ds12(is)
              enddo
            enddo
          enddo

!          call THERMOADV

          do is = 1, ns-1 ! It is asumed that there is no MPI decomposition of the domain in vertical
            do iy = max(nyi,2), nye   !2, ny
              do ix = max(nxi,2), nxe !2, nx
                q = (s(ix,iy,is) + s(ix,iy,is+1)) &
                & * (pts(ix,iy,is+1) - pts(ix,iy,is-1)) &
                & * ((w(ix,iy,is+1,2) + w(ix,iy,is,2)) - &
                & fcor/g*(vgeos(is)*(u(ix,iy,is,2) + u(ix-1,iy,is,2) ) - &
                & ugeos(is)*(v(ix,iy,is,2) + v(ix,iy-1,is,2) )))/ ds08a(is)
                tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
                & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
                & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
                & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
                tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
                & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
                & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
                & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
                tavs(ix,iy,is) = (wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2)) &
                & -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2))) &
                & /ds12(is)
                call MICROPHYS_TERMS
                acoef(ix,iy,is) = - 0.5*dtl*st(ix,iy,is)*difs001(ix,iy,is)*pp(ix,iy,2)
                bcoef(ix,iy,is) = - 0.5*dtl*st(ix,iy,is)*difs001(ix,iy,is+1)*pp(ix,iy,2)
                ccoef(ix,iy,is) = acoef(ix,iy,is) + bcoef(ix,iy,is) - pp(ix,iy,3)
                dcoef(ix,iy,is) = - pt(ix,iy,is,1)*pp(ix,iy,1) + &
                & dtl*(tavx(ix,iy,is) + tavy(ix,iy,is) + &
                & pp(ix,iy,2)*(refsadjust*dptsdtsig(ix,iy,is) + &
                & tavs(ix,iy,is) - q - condens - transform  - difuntcg(ix,iy,is) - &
                & st(ix,iy,is) * &
                & (0.5*difs001(ix,iy,is+1)*pt(ix,iy,is+1,1) - &
                &  0.5*(difs001(ix,iy,is+1) + difs001(ix,iy,is))*pt(ix,iy,is,1) + &
                &  0.5*difs001(ix,iy,is)*pt(ix,iy,is-1,1) + &
                &  tsdif * &
                &  (difs001(ix,iy,is+1)*pts(ix,iy,is+1) - &
                &  (difs001(ix,iy,is+1) + difs001(ix,iy,is))*pts(ix,iy,is) + &
                &  difs001(ix,iy,is)*pts(ix,iy,is-1)) ) ))
!                if (nstep >= 700 .and. is == 1 .and. ix == (nx-1)/2 .and. &
!                  & iy == (ny-1)/2 ) then
!                  print*, ix,iy,is,st(ix,iy,is),difs001(ix,iy,is:is+1),&
!                  & pp(ix,iy,1:3),pt(ix,iy,is-1:is+1,1),tavx(ix,iy,is),tavy(ix,iy,is),&
!                  & tavs(ix,iy,is), q, condens, transform, difuntcg(ix,iy,is)
!                  read*
!                endif
              enddo
            enddo
          enddo

          do iy = max(nyi,2), nye
            do ix = max(nxi,2), nxe
!              acoef(ix,iy,ns) = 1.
!              ccoef(ix,iy,ns) = 1.
!              dcoef(ix,iy,ns) = 0.
              acoef(ix,iy,ns) = 0.5*dift001(ix,iy,ns)/deltaz(ix,iy)
              ccoef(ix,iy,ns) = 0.5*dift001(ix,iy,ns)/deltaz(ix,iy)
              p = sigma0(ns-1)*pp(ix,iy,2) + ptop
              dcoef(ix,iy,ns) = h(ix,iy)*(p00/p)**akapa/(cp*rhoasurf(ix,iy))
              ccoef(ix,iy,0) = 1.
              bcoef(ix,iy,0) = 1.
              dcoef(ix,iy,0) = 0.
            enddo
          enddo

          call PROGONKAS2D &
          & (acoef, bcoef, ccoef, dcoef, pt(nxs2i,nys2i,0,3), &
          & nxs2i, nxs2e, nys2i, nys2e, 0, ns1, &
          & max(nxi,2), nxe, max(nyi,2), nye, 0, ns)

!                do is = ns-1,0,-1
!                do iy = max(nyi,2), nye
!                do ix = max(nxi,2), nxe
!                if (abs(pt(ix,iy,is,3)) >= 50.) then
!                  print*, 'thermo:max exceeded', 
!                  & pt(ix,iy,is,3),ix,iy,is, &
!                  & st(ix,iy,is-1:is+1), &
!                  & difs001(ix,iy,is-1:is+1), &
!                  & pp(ix,iy,1:2), &
!                  & pt(ix,iy,is,1:2), &
!                  & tavx(ix,iy,is), &
!                  & tavy(ix,iy,is), &
!                  & tavs(ix,iy,is), &
!                  & difuntcg(ix,iy,is), &
!                  & h(ix,iy), deltaz(ix,iy)
!                  read*
!                endif
!                enddo
!                enddo
!                enddo

        else

!         call THERMOADV

          do is = 1, ns-1 ! It is asumed that there is no MPI decomposition of the domain in vertical
            do iy = max(nyi,2), nye   !2, ny
              do ix = max(nxi,2), nxe !2, nx
                q = (s(ix,iy,is) + s(ix,iy,is+1)) &
                & * (pts(ix,iy,is+1) - pts(ix,iy,is-1)) &
                & * ((w(ix,iy,is+1,2) + w(ix,iy,is,2)) - &
                & fcor/g*(vgeos(is)*(u(ix,iy,is,2) + u(ix-1,iy,is,2) ) - &
                & ugeos(is)*(v(ix,iy,is,2) + v(ix,iy-1,is,2) )))/ ds08a(is)
                tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
                & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
                & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
                & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
                tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
                & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
                & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
                & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
                tavs(ix,iy,is) = (wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2)) &
                & -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2))) &
                & /ds12(is)
                call MICROPHYS_TERMS
                pt(ix,iy,is,3)=(pt(ix,iy,is,1)*pp(ix,iy,1) &
                & -dtl*(tavx(ix,iy,is)+tavy(ix,iy,is)+pp(ix,iy,2) &
                & *(refsadjust*dptsdtsig(ix,iy,is)+tavs(ix,iy,is)-q-difunt(ix,iy,is) &
                & -condens-transform)))/pp(ix,iy,3)
              enddo
            enddo
          enddo

        endif

      endif


!        write (work_char,*) nstep
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt1'//work_char)

!     Debug output of the array
!      if (nstep == 1) then
!        call ARRGATHERWRI_MPI(phi(nxs2i,nys2i,0), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'phi')
!        call ARRGATHERWRI_MPI(phis(nxs2i,nys2i,0), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'phis')
!        call ARRGATHERWRI_MPI(alb(nxs2i,nys2i), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,0, &
!        & 0,nx1,0,ny1,0,0,nx1,ny1,ns1,'alb')
!        call ARRGATHERWRI_MPI(emis(nxs2i,nys2i), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,0, &
!        & 0,nx1,0,ny1,0,0,nx1,ny1,ns1,'emis')


!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt1')
        !STOP
!      endif
! Debug output of the array
!      call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), & !pt(nxs2i,nys2i,0,3), &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & 1,nx+1,1,ny+1,0,ns,nx1,ny1,ns1)
!      STOP

     
      
!------Heating due to shortwave and longwave radiation absorption-------------         
      ifradpar: if (radpar == 3) then
!------Clirad subroutines SORAD and IRRAD are called once   ------
!------per nradcall time steps,                             ------
!------in between Radheat is not updated -------------------------   
        radcall: if (mod(nstep,nradcall) == 0 .or. nstep == 1) then
          call RAD_FLUX_UPDATE &
          & (nx, ny, ns, &
          & imonth,iday,ihour,iminu,iseco,nstep,xlatit, &
          & isqv, ifqc, ifqr, shortwave, longwave, &
          & sigma1,sigma0,ds1,pp,ptop,phis,phi, &
          & pts,pt,qv,qc,qr,qci,qsn, &
          & hsuf,tsurf,alb,emis, &
          & Srad,Srad_surf,Lrad,Lrad_surf, &
          & SraddirUV_surf,SraddirPAR_surf,SraddirIR_surf, &
          & SraddifUV_surf,SraddifPAR_surf,SraddifIR_surf)
        endif radcall

        call EXTRAH_MPI &
        & (Srad_surf, nxs2i, nxs2e, nys2i, nys2e, 0, 0, &
        & 0, nx1, 0, ny1, 0, 0)
        call EXTRAH_MPI &
        & (Lrad_surf, nxs2i, nxs2e, nys2i, nys2e, 0, 0, &
        & 0, nx1, 0, ny1, 0, 0)

!      write(*,*) 'After radiation'
!      call MPI_BARRIER(comm3d,ix)
!      write (*,*) 'The barrier 9 passed'
!      STOP
              
!      Positive radiation fluxes are assumed to be downward

        do is = 1, ns-1
          do iy = max(nyi,2), nye   !2, ny
            do ix = max(nxi,2), nxe !2, nx
              p = sigma0(is)*pp(ix,iy,2) + ptop
              temperature = (pts(ix,iy,is) + pt(ix,iy,is,2))/(p00/p)**akapa
              Radheat(ix,iy,is) = - 0.5*(s(ix,iy,is) + s(ix,iy,is+1)) * &
              & (Srad(ix,iy,is+1)-Srad(ix,iy,is) + &
              &  Lrad(ix,iy,is+1)-Lrad(ix,iy,is)) / &
              &  ds1(is)/cp/p*temperature*r*(p00/p)**akapa
              pt(ix,iy,is,3) = pt(ix,iy,is,3) + dtl * &
              & pp(ix,iy,2)*Radheat(ix,iy,is)/pp(ix,iy,3)
            enddo
          enddo
        enddo
     
      endif ifradpar

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nxi,nx1e,nyi,ny1e,1,ns-1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'pt exceeded min or max after progonka', i, j, k, &
!        & 'pt_1:3',pt(i,j,k,1:3), &
!        & 'difs_k:k+1',difs001(i,j,k:k+1), &
!        & 'dift001_k:k+1', dift001(i,j,k:k+1), &
!        & 's_k:k+1', s(i,j,k:k+1), &
!        & 'pp_1:3',pp(i,j,1:3), 'tavx_tend', dtl*tavx(i,j,k), &
!        & 'tavy_tend', dtl*tavy(i,j,k), 'tavs_tend', &
!        & dtl*pp(i,j,2)*tavs(i,j,k), 'pts_k-1:k+1', pts(i,j,k-1:k+1), &
!        & 'w2_k:k+1', w(i,j,k:k+1,2), 'Radheat_tend', Radheat(i,j,k)
!      endif

!      write (work_char,*) nstep
!      call ARRGATHERWRI_MPI(Srad(nxs2i,nys2i,0), &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'Srad'//work_char)
!      write (work_char,*) nstep
!      call ARRGATHERWRI_MPI(Lrad(nxs2i,nys2i,0), &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'Lrad'//work_char)
!      STOP

!        write (work_char,*) nstep
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt2'//work_char)


!     Debug output of the array
!      if (nstep == 1) then
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt2')
!        call ARRGATHERWRI_MPI(Srad(nxs2i,nys2i,0), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'Srad')
!        call ARRGATHERWRI_MPI(Lrad(nxs2i,nys2i,0), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'Lrad')

        !STOP
!      endif

!      call MPI_BARRIER(comm3d,ix)
!      write(*,*) 'flag1'

! surface heat flux (bulk parametrization):

      if (ifsoil /= 0 .or. ifhle == 1) then
        if (.not. impldiff) then
          do iy = max(nyi,2), nye    !2,ny
            do ix = max(nxi,2), nxe  !2,nx
              p = sigma0(ns-1)*pp(ix,iy,2) + ptop
              pt(ix,iy,ns-1,3) = pt(ix,iy,ns-1,3) &
              & +dtl*h(ix,iy)*(p00/p)**akapa &
              & /(2.*deltaz(ix,iy)*cp*rhoasurf(ix,iy)) ! 2.* - bug corrected
            enddo
          enddo
        endif
      elseif (cdcoef > 0.) then
        ctbulk = cdcoef*dtl
        do iy = max(nyi,2), nye   !2,ny
          do ix = max(nxi,2), nxe !2,nx
            p = sigma0(ns-1)*pp(ix,iy,2) + ptop
            if (xlake(ix,iy).eq.0.) then
              tsfunc = TSOIL(ix,iy,timeho)
            else
              tslake(ix,iy) = TLAKE(ix,iy,timeho)
              tsfunc = tslake(ix,iy)
            endif
            pt(ix,iy,ns-1,3) = pt(ix,iy,ns-1,3)+ctbulk*(p00/p)**akapa &
            & *uvsuf(ix,iy)*(tsfunc-pts(ix,iy,ns-1) &
            & -pt(ix,iy,ns-1,2))/deltaz(ix,iy)
          enddo
        enddo
      endif

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nxi,nx1e,nyi,ny1e,1,ns-1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'pt exceeded min or max surface forcing', i, j, k, &
!        & 'pt_1:3',pt(i,j,k,1:3)
!      endif

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nx0i,nx1e,ny0i,ny1e,0,ns1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'theta exceeded min or max after progonka', i, j, k, &
!        & pt(i,j,k,1:3), pt(i,j,k-1:k+1,1), &
!        & difs001(i,j,k), dift001(i,j,k), s(i,j,k:k+1), pts(i,j,k-1:k+1), &
!        & w(i,j,k:k+1,2), pp(i,j,1:3), h(i,j)
!      endif

!     Debug output of the array
!      if (nstep == 1) then
!        call ARRGATHERWRI_MPI(h(nxs2i,nys2i), &
!       & nxs2i,nxs2e,nys2i,nys2e,0,0, &
!        & 0,nx1,0,ny1,0,0,nx1,ny1,ns1,'h')
!        call ARRGATHERWRI_MPI(deltaz(nxs2i,nys2i), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,0, &
!        & 0,nx1,0,ny1,0,0,nx1,ny1,ns1,'deltaz')

!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt3')
        !STOP
!      endif

!      call MPI_BARRIER(comm3d,ix)
!      write(*,*) 'flag2'


! horizontal boundary conditions:

!      call PTBC_MPI(nx,ny,ns,pt(nxs2i,nys2i,0,3),ptbx,ptby,ptcc, &
!      & iobptx,iobpty,ptbout,prt,nchn1)
      call PTBC2_MPI(pt,pp,ptbx,ptby,ptcc, &
      & ptbxe,ptbye,ubxe,vbxe,ubye,vbye, &
      & pts,phi,s,ds0, &
      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
      & nx,ny,ns,iobptx,iobpty,ptbout,prt,nchn1)

!      call MPI_BARRIER(comm3d,ix)
!      write(*,*) 'flag3'

!       write(*,*) 'After ptbc'
!       STOP

! vertical boundary conditions:

      do iy = nyi, ny1e   !1,ny1
        do ix = nxi, nx1e !1,nx1
          pt(ix,iy,ns,3) = pt(ix,iy,ns-1,3)
          pt(ix,iy,0,3) = pt(ix,iy,1,3)
        enddo
      enddo

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nxi,nx1e,nyi,ny1e,1,ns-1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'pt exceeded min or max after b.c.s', i, j, k, &
!        & 'pt_1:3',pt(i,j,k,1:3)
!      endif

!        write (work_char,*) nstep
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt3'//work_char)


! Raleigh damping
      if (raylei) then
        do is = 1, idrmax
          do iy = nyi, ny1e    !1, ny1
            do ix = nxi, nx1e  !1, nx1
              pt(ix,iy,is,3) = pt(ix,iy,is,3) - dtl/taudra(is)*pt(ix,iy,is,2)
            enddo
          enddo
        enddo
      endif


!     Debug output of the array
!      if (nstep == 1) then
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt4')
!        STOP
!      endif
! smoothing:

      if (dohsmo .and. (mod(nstep,numsmo).eq.0) ) then
!        call MPI_BARRIER(comm3d,ix)
!        write(*,*) 'flag4-01'
        if (parallel) then
!          pr_time = MPI_WTIME() ! Timing

          prplc(1:4) = 1
          call ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
          &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
          &  nxi,nx1e,nyi,ny1e,1,ns-1, &
          &  nxsh2,nysh2,0, &
          &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
          &  size_MPI,comm3d,req,ndim,prplc, &
          &  request_send_3d,request_recv_3d)
          ! Exchanges for 4-order filtering at boundaries for periodic l.b.c.s
          if (iobptx == lbcs%per) then
            prplc(1) = 1
            call ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  nxi,min(nxe,nx-1),nyi,ny1e,1,ns, &
            &  nxsh2,0,0, &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  size_MPI,comm3d,req,ndim,prplc, &
            &  request_send_3d,request_recv_3d,1)
            prplc(2) = 1
            call ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  max(nxi,3),nx1e,nyi,ny1e,1,ns, &
            &  nxsh2,0,0, &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  size_MPI,comm3d,req,ndim,prplc, &
            &  request_send_3d,request_recv_3d,1)
          endif
          if (iobpty == lbcs%per) then
            prplc(3) = 1
            call ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  nxi,nx1e,nyi,min(nye,ny-1),1,ns, &
            &  0,nysh2,0, &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  size_MPI,comm3d,req,ndim,prplc, &
            &  request_send_3d,request_recv_3d,1)
            prplc(4) = 1
            call ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  nxi,nx1e,max(nyi,3),ny1e,1,ns, &
            &  0,nysh2,0, &
            &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
            &  size_MPI,comm3d,req,ndim,prplc, &
            &  request_send_3d,request_recv_3d,1)
          endif
        else
          if (iobptx == lbcs%per) then
            do is = 1, ns
              do iy = nyi, ny1e
                pt(  -1:0   ,iy,is,3) = pt(nx-2:nx-1,iy,is,3)
                pt(nx+2:nx+3,iy,is,3) = pt(   3:4   ,iy,is,3)
              enddo
            enddo
          endif
          if (iobpty == lbcs%per) then
            do is = 1, ns
              do ix = nxi, nx1e
                pt(ix,  -1:0   ,is,3) = pt(ix,ny-2:ny-1,is,3)
                pt(ix,ny+2:ny+3,is,3) = pt(ix,   3:4   ,is,3)
              enddo
            enddo
          endif
!          proc_time(2,2) = proc_time(2,2) + MPI_WTIME() - pr_time ! Timing
        endif
!        call MPI_BARRIER(comm3d,ix)
!        write(*,*) 'flag4-0'
!        call HSMOOT_PAR(pt(0,0,0,3),hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1,hdamp)
        call HSMOOT_MPI(pt(nxs2i,nys2i,0,3),hk1, &
        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
        & nxi,nx1e,nyi,ny1e,1,ns-1,hdamp, &
        & iobptx == lbcs%per, iobpty == lbcs%per)
      endif

!     Debug output of the array
!      if (nstep == 1) then
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt')
!       STOP
!      endif


      if (dovsmt .and. (mod(nstep,numsmo).eq.0) ) then
!        call VSMOOT_PAR(pt(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1,vdampt)
        call VSMOOT_MPI(pt(nxs2i,nys2i,0,3),wrk,nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
        & nxi,nx1e,nyi,ny1e,0,ns,vdampt) ! 1, ns-1
      endif


      if (.not.tfct) then
        ! Check the procedure arguments
!        call ASELIN_PAR(pt(0,0,0,3),pt(0,0,0,2),ptm2,0,nx1,0,ny1,0,ns,filt)
        call ASELIN_MPI(pt(nxs2i,nys2i,0,3),pt(nxs2i,nys2i,0,2),pt(nxs2i,nys2i,0,1), &
        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
        & nxi,nx1e,nyi,ny1e,0,ns,filt, &
        & xlbound, xrbound, ylbound, yrbound)
      endif

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nxi,nx1e,nyi,ny1e,1,ns-1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'pt exceeded min or max after smoothing', i, j, k, &
!        & 'pt_1:3',pt(i,j,k,1:3)
!      endif

!        write (work_char,*) nstep
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt4'//work_char)


!      call MPI_BARRIER(comm3d,ix)
!      write(*,*) 'flag4'

!     Debug output of the array
!      call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,2), &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & 0,nx+1,0,ny+1,0,ns,nx1,ny1,ns1)
!      STOP

      if (parallel) then
        ! Note this exchange is a subject for opimization since
        ! it performs much more exchanges than those needed by RADBCH_MPI

!        pr_time = MPI_WTIME()

!        prplc(1:4) = 1       
!        call  ARR3DEXCH(pt(nxs2i,nys2i,0,3), &
!        &  nxs2i, nxs2e, nys2i, nys2e, 0, ns1, &
!        &  nxi, nx1e, nyi, ny1e, 1, ns-1, &
!        &  nxsh, nysh, 0, &
!        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        &  size_MPI, comm3d, req, ndim, prplc, &
!        &  request_send_3d,request_recv_3d)
!        prplc(1:4) = 1
!        call  ARR3DEXCH(pt(nxs2i,nys2i,0,2), &
!        &  nxs2i, nxs2e, nys2i, nys2e, 0, ns1, &
!        &  nxi, nx1e, nyi, ny1e, 1, ns-1, &
!        &  nxsh, nysh, 0, &
!        &  nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        &  size_MPI, comm3d, req, ndim, prplc, &
!        &  request_send_3d,request_recv_3d)

        call RADBCHEXCH(pt(nxs2i,nys2i,0,2), &
        &  nxs2i, nxs2e, nys2i, nys2e, 0, ns1, &
        &  nxi, nx1e, nyi, ny1e, 1, ns-1, 2)
        call RADBCHEXCH(pt(nxs2i,nys2i,0,3), &
        &  nxs2i, nxs2e, nys2i, nys2e, 0, ns1, &
        &  nxi, nx1e, nyi, ny1e, 1, ns-1, 2)
!        proc_time(3,2) = proc_time(3,2) + MPI_WTIME() - pr_time

      endif

!      call radbch(ptcc,ptbx,ptby,pt(0,0,0,3),pt(0,0,0,2),ptm2, &
!      & nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)

       call RADBCH_MPI &
       & (ptcc,ptbx,ptby, &
       & pt(nxs2i,nys2i,0,3),pt(nxs2i,nys2i,0,2),pt(nxs2i,nys2i,0,1), &
       & s, &
       & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
       & 1,ns-1,nxi,nx1e,nyi,ny1e,4)

!        write (work_char,*) nstep
!        call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
!        & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!        & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt5'//work_char)

       if (nxsponge2 > 0 .and. nxsponge2 >= nxi) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nxsponge2:nxsponge2,ny0i:ny1e,0:ns))
         call SPONGE_EXCH(nxsponge2,nxsponge2,ny0i,ny1e,0,ns, &
         & nxsponge2,nxs2e,1)
         endif

         if (ntypesponge == 2) then
         do is = 0, ns
           do iy = ny0i, ny1e !0, ny1
             do ix = nx0i, min(nxe,nxsponge2)
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - &
               & (pt(ix,iy,is,2) - work_arrsectget(nxsponge2,iy,is))*dtl/tauspo(ix)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         elseif (ntypesponge == 1) then
         do is = 0, ns
           do iy = ny0i, ny1e !0, ny1
             do ix = nx0i, min(nxe,nxsponge2)
               pt(ix,iy,is,3) = pt(ix,iy,is,3)*(1.-dtl/tauspo(ix))
             enddo
           enddo
         enddo

         elseif (ntypesponge == 3) then
         do is = 0, ns
           do iy = ny0i, ny1e !0, ny1
             do ix = nx0i, min(nxe,nxsponge2)
               if (ubxe(1,iy,is) > 0.) then
                 q = ptbxe(1,iy,is)
               else
                 q = work_arrsectget(nxsponge2,iy,is)
               endif
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - (pt(ix,iy,is,2) - q)*dtl/tauspo(ix)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)
         endif

         !pt(ix,iy,is,3) = pt(ix,iy,is,3)*(1.-dtl/tauspo(ix))
       elseif (nxsponge2 > 0) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nxsponge2:nxsponge2,ny0i:ny1e,0:ns))
         call SPONGE_EXCH(nxsponge2,nxsponge2,ny0i,ny1e,0,ns, &
         & nxsponge2,nxs2e,2)
         deallocate (work_arrsectget)
         endif

       endif

       if (nxsponge1 > 0 .and. nx1-nxsponge1 <= nxe) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx1-nxsponge1:nx1-nxsponge1,ny0i:ny1e,0:ns))
         call SPONGE_EXCH(nx1-nxsponge1,nx1-nxsponge1,ny0i,ny1e,0,ns, &
         & nxs2i,nx1-nxsponge1,1)
         endif

         if (ntypesponge == 2) then
         do is = 0, ns
           do iy = ny0i, ny1e 
             do ix = max(nx0i,nx1-nxsponge1), nx1e
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - &
               & (pt(ix,iy,is,2) - work_arrsectget(nx1-nxsponge1,iy,is))*dtl/tauspo(ix)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         elseif (ntypesponge == 1) then
         do is = 0, ns
           do iy = ny0i, ny1e 
             do ix = max(nx0i,nx1-nxsponge1), nx1e
               pt(ix,iy,is,3) = pt(ix,iy,is,3)*(1.-dtl/tauspo(ix))
             enddo
           enddo
         enddo


         elseif (ntypesponge == 3) then
         do is = 0, ns
           do iy = ny0i, ny1e 
             do ix = max(nx0i,nx1-nxsponge1), nx1e
               if (ubxe(2,iy,is) < 0.) then
                 q = ptbxe(2,iy,is)
               else
                 q = work_arrsectget(nx1-nxsponge1,iy,is)
               endif
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - (pt(ix,iy,is,2) - q)*dtl/tauspo(ix)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         endif

       elseif (nxsponge1 > 0) then
         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx1-nxsponge1:nx1-nxsponge1,ny0i:ny1e,0:ns))
         call SPONGE_EXCH(nx1-nxsponge1,nx1-nxsponge1,ny0i,ny1e,0,ns, &
         & nxs2i,nx1-nxsponge1,2)
         deallocate (work_arrsectget)
         endif
       endif

       if (nysponge2 > 0 .and. nysponge2 >= nyi) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx0i:nx1e,nysponge2:nysponge2,0:ns))
         call SPONGE_EXCH(nx0i,nx1e,nysponge2,nysponge2,0,ns, &
         & nysponge2,nys2e,1)
         endif

         if (ntypesponge == 2) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = ny0i, min(nye,nysponge2)
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - &
               & (pt(ix,iy,is,2) - work_arrsectget(ix,nysponge2,is))*dtl/tauspoy(iy)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         elseif (ntypesponge == 1) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = ny0i, min(nye,nysponge2)
               pt(ix,iy,is,3) = pt(ix,iy,is,3)*(1.-dtl/tauspoy(iy))
             enddo
           enddo
         enddo

         elseif (ntypesponge == 3) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = ny0i, min(nye,nysponge2)
               if (vbye(ix,1,is) > 0.) then
                 q = ptbye(ix,1,is)
               else
                 q = work_arrsectget(ix,nysponge2,is)
               endif
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - (pt(ix,iy,is,2) - q)*dtl/tauspoy(iy)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         endif

       elseif (nysponge2 > 0) then
         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx0i:nx1e,nysponge2:nysponge2,0:ns))
         call SPONGE_EXCH(nx0i,nx1e,nysponge2,nysponge2,0,ns, &
         & nysponge2,nys2e,2)
         deallocate (work_arrsectget)
         endif
       endif

       if (nysponge1 > 0 .and. ny1-nysponge1 <= nye) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx0i:nx1e,ny1-nysponge1:ny1-nysponge1,0:ns))
         call SPONGE_EXCH(nx0i,nx1e,ny1-nysponge1,ny1-nysponge1,0,ns, &
         & nys2i,ny1-nysponge1,1)
         endif

         if (ntypesponge == 2) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = max(ny0i,ny1-nysponge1), ny1e
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - &
               & (pt(ix,iy,is,2) - work_arrsectget(ix,ny1-nysponge1,is))*dtl/tauspoy(iy)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         elseif (ntypesponge == 1) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = max(ny0i,ny1-nysponge1), ny1e
               pt(ix,iy,is,3) = pt(ix,iy,is,3)*(1.-dtl/tauspoy(iy))
             enddo
           enddo
         enddo


         elseif (ntypesponge == 3) then
         do is = 0, ns
           do ix = nx0i, nx1e !0, nx1
             do iy = max(ny0i,ny1-nysponge1), ny1e
               if (vbye(ix,2,is) < 0.) then
                 q = ptbye(ix,2,is)
               else
                 q = work_arrsectget(ix,ny1-nysponge1,is)
               endif
               pt(ix,iy,is,3) = pt(ix,iy,is,3) - (pt(ix,iy,is,2) - q)*dtl/tauspoy(iy)
             enddo
           enddo
         enddo
         deallocate (work_arrsectget)

         endif

       elseif (nysponge1 > 0) then

         if (ntypesponge == 2 .or. ntypesponge == 3) then
         allocate (work_arrsectget(nx0i:nx1e,ny1-nysponge1:ny1-nysponge1,0:ns))
         call SPONGE_EXCH(nx0i,nx1e,ny1-nysponge1,ny1-nysponge1,0,ns, &
         & nys2i,ny1-nysponge1,2)
         deallocate (work_arrsectget)
         endif

       endif

       if (iobptx /= lbcs%per .and. iobpty /= lbcs%per) then
         call EXTRAH_MPI &
         & (pt(nxs2i,nys2i,0,3),nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
         & 0,nx2,0,ny2,1,ns-1)
       else
         if (iobptx /= lbcs%per) then
           call EXTRAX_MPI &
           & (pt(nxs2i,nys2i,0,3),nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
           & 0,nx2,0,ny2,1,ns-1)
         endif
         if (iobpty /= lbcs%per) then
           call EXTRAY_MPI &
           & (pt(nxs2i,nys2i,0,3),nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
           & 0,nx2,0,ny2,1,ns-1)
         endif
         ! For output only
         if (xlbound .and. ylbound) pt(0,0,1:ns-1,3) = pt(1,1,1:ns-1,3)
       endif

!      call CHECKMINMAX &
!      & (pt(nxs2i,nys2i,0,3),-5.d0,5.d0, &
!      & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
!      & nxi,nx1e,nyi,ny1e,1,ns-1, &
!      & i,j,k,indunst)
!      if (indunst > 0 .and. &
!        & i /= nodata .and. j /= nodata .and. k /= nodata) then
!        write(*,*) 'pt exceeded min or max after sponge', i, j, k, &
!        & 'pt_1:3',pt(i,j,k,1:3)
!      endif

! vertical boundary conditions:

!      do iy=1,ny1
!      do ix=1,nx1
!        pt(ix,iy,ns,3)=pt(ix,iy,ns-1,3)
!        pt(ix,iy,0,3)=pt(ix,iy,1,3)
!      enddo
!      enddo

!     Debug output of the array
!      if (mod(nstep,360) == 0) then ! .or. nstep >= 6970) then ! .or. nstep >= 2660) then
!       if (debug_output) then
#ifdef debug1       
         write (work_char,'(i10.10)') nstep
         call ARRGATHERWRI_MPI(pt(nxs2i,nys2i,0,3), &
         & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
         & 0,nx1,0,ny1,0,ns1,nx1,ny1,ns1,'pt3_thermo'//work_char)
#endif         
!       endif 
!        STOP
!      endif

      if (impldiff) then
        deallocate (difs001)
        deallocate (st)
        deallocate (acoef)
        deallocate (bcoef)
        deallocate (ccoef)
        deallocate (dcoef)
      endif

      deallocate (tavx,tavy,tavs)

      return
      contains

      SUBROUTINE SPONGE_EXCH(i0,i1,j0,j1,k0,k1,lcoord,rcoord,proc)

!     This surboutine performs exchanges needed for sponge
!     at lateral boundaries

      implicit none

      integer(4), intent(in) :: i0,i1,j0,j1,k0,k1
      integer(4), intent(in) :: lcoord, rcoord
      integer(4), intent(in) :: proc

      if (proc == 1) then
        if (parallel) then
          call ARRSECTGET &
          & (pt(nxs2i,nys2i,0,2), &
          & work_arrsectget(i0,j0,k0), &
          & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
          & nx0i,nx1e,ny0i,ny1e,0,ns1, &
          & i0,i1,j0,j1,k0,k1, &
          & nxsh, nysh, nssh, &
          & comm3d, size_MPI, coords, dims, &
          & neighborhood, isperiodic)
          if (lcoord <= rcoord) then
            work_arrsectget(i0:i1,j0:j1,k0:k1) = pt(i0:i1,j0:j1,k0:k1,2)
          endif
        else
          work_arrsectget(i0:i1,j0:j1,k0:k1) = pt(i0:i1,j0:j1,k0:k1,2)
        endif
      elseif (proc == 2) then
        if (parallel) then
          call ARRSECTGET &
          & (pt(nxs2i,nys2i,0,2), &
          & work_arrsectget(i0,j0,k0), &
          & nxs2i,nxs2e,nys2i,nys2e,0,ns1, &
          & nx0i,nx1e,ny0i,ny1e,0,ns1, &
          & nodata,i1,j0,j1,k0,k1, &
          & nxsh, nysh, nssh, &
          & comm3d, size_MPI, coords, dims, &
          & neighborhood, isperiodic)
        endif
      endif

      END SUBROUTINE SPONGE_EXCH


      SUBROUTINE MICROPHYS_TERMS

      p = sigma0(is)*pp(ix,iy,2)+ptop
      if (ifqc .ne. 0 .and. ifqi .eq. 0) then
        condens = hlatcp*(p00/p)**akapa*(cond(ix,iy,is) &
        & - evap(ix,iy,is))
      endif
      if(ifqc.ne.0 .and. ifqi.ne.0) then
        condens = hlatcp*(p00/p)**akapa*(-evap(ix,iy,is))
      endif
      if (ifqc == 0 ) then
        condens = 0.
      endif
      if(ifqi /= 0) then
        transform = &
        & hsubcp*(p00/p)**akapa*(vini(ix,iy,is)) + &
        & hfuscp*(p00/p)**akapa* &
        & (hmfrz(ix,iy,is) + sacrw(ix,iy,is) + iacr(ix,iy,is) + &
        & sbercw(ix,iy,is) + ibercw(ix,iy,is) + sacrr(ix,iy,is) - &
        & smlt(ix,iy,is) - imlt(ix,iy,is))
      else
        transform = 0.
      endif

      END SUBROUTINE MICROPHYS_TERMS


      SUBROUTINE THERMOADV

      use MPI_VARIABLES, only : &
      & xlbound, xrbound, &
      & ylbound, yrbound

      implicit none

!     Subroutine calculates scalar advection terms using 
!     centered differences at C-Arakawa grid and forward differences at
!     lateral boundaries where the flow is directed outwards the model domain

!      real(8), intent(out) :: tavx_(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
!      real(8), intent(out) :: tavy_(nxs2i:nxs2e,nys2i:nys2e,0:ns1)
!      real(8), intent(out) :: tavs_(nxs2i:nxs2e,nys2i:nys2e,0:ns1)

      do is = 1, ns-1 ! It is asumed that there is no MPI decomposition of the domain in vertical
        do iy = max(nyi,3), min(nye,ny-1)   !2, ny
          do ix = max(nxi,3), min(nxe,nx-1) !2, nx
            tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
            & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
            & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
            & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
            tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
            & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
            & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
            & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
          enddo
        enddo
      enddo


      if (xlbound) then

        ix = 2
        do is = 1, ns-1
          do iy = max(nyi,2), nye !1, ny1
            if (iy /= 2 .and. iy /= ny) then
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
              & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
              & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
              & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
            endif
            if (u(ix,iy,is,3) > 0.) then
              tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
              & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
              & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
              & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
            else
              tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
              & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
              & -(u(ix-1,iy,is,2)*pp10(ix-1,iy,2) + &
              &   u(ix,iy,is,2)*pp10(ix,iy,2))*pt(ix,iy,is,2))/dx
            endif
          enddo
        enddo

      endif

      if (xrbound) then

        ix = nx
        do is = 1, ns-1
          do iy = max(nyi,2), nye !1, ny1
            if (iy /= 2 .and. iy /= ny) then
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
              & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
              & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
              & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
            endif
            if (u(ix,iy,is,3) < 0.) then
              tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
              & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
              & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
              & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
            else
              tavx(ix,iy,is) = ((u(ix-1,iy,is,2)*pp10(ix-1,iy,2) + &
              & u(ix,iy,is,2)*pp10(ix,iy,2))*pt(ix,iy,is,2) &
              & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
              & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx
            endif
          enddo
        enddo

      endif

      if (ylbound) then

        iy = 2
        do is = 1, ns-1
          do ix = max(nxi,2), nxe !1, ny1
            if (ix /= 2 .and. ix /= nx) then
              tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
              & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
              & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
              & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
            endif
            if (v(ix,iy,is,3) > 0.) then
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
              & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
              & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
              & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
            else
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
              & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
              & -(v(ix,iy-1,is,2)*pp01(ix,iy-1,2) + &
              &   v(ix,iy,is,2)*pp01(ix,iy,2))*pt(ix,iy,is,2))/dy
            endif
          enddo
        enddo

      endif

      if (yrbound) then

        iy = ny
        do is = 1, ns-1
          do ix = max(nxi,2), nxe !1, ny1
            if (ix /= 2 .and. ix /= nx) then
              tavx(ix,iy,is) = (u(ix,iy,is,2)*pp10(ix,iy,2) &
              & *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2)) &
              & -u(ix-1,iy,is,2)*pp10(ix-1,iy,2) &
              & *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
            endif
            if (v(ix,iy,is,3) < 0.) then
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) &
              & *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2)) &
              & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
              & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
            else
              tavy(ix,iy,is) = (v(ix,iy,is,2)*pp01(ix,iy,2) + &
              & v(ix,iy-1,is,2)*pp01(ix,iy-1,2))*pt(ix,iy,is,2) &
              & -v(ix,iy-1,is,2)*pp01(ix,iy-1,2) &
              & *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2))/dy
            endif
          enddo
        enddo

      endif

      do is = 1, ns-1 ! It is asumed that there is no MPI decomposition of the domain in vertical
        do iy = max(nyi,2), nye  !2, ny
          do ix = max(nxi,2), nxe  !2, nx
            tavs(ix,iy,is) = (wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2)) &
            & -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2))) &
            & /ds12(is)
          enddo
        enddo
      enddo

      END SUBROUTINE THERMOADV


      END SUBROUTINE THERMO_MPI
