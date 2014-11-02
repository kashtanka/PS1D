      module alloc_1d
      
      integer:: ns,ndat,npro,nz,ntime
      integer:: nstep
      integer:: ifhle,tfix  
      real(kind(0.d0)):: dt,dtl
      integer::
     :  ioreft,iorefq,iorefu,iorefv,ndth,ndus,ndvs,ndqvs
      integer :: iodif
      real(kind(0.d0)):: pts0,pa,fcor
      real(kind(0.d0)):: t0ini,ug,vg,psref,dzpro
      real(kind(0.d0)):: cdm,ust_s,tst_s,qst_s,dzits,ts,hbl,Fv,h,le,zct
      real(kind(0.d0)):: zc2,dR
      real(kind(0.d0)):: r,cp,g,akapa,qif,p00,omega,pi,entrt,hlat,ifwr,
     : ifmf 
      real(kind(0.d0)):: z_sl,ztop,distY,distYice,icetime
      real(kind(0.d0)):: ablv,dy,p0,dpdy_d,phi,ptop,bl_dpdy
      real(kind(0.d0)):: cond_heat,mth,mqv,mqc,mqr,mqci,mqsn
      parameter(r=287.05,cp=1005.,g=9.8066,akapa=r/cp,p00=1.e5,
     :          omega=7.2921e-5,pi=3.141593,hlat=2.501e6,
     :          hsub=2.837e6,hfus=3.336e5)
      
      real(kind(0.d0)),allocatable,dimension(:)::
     :   thdat,usdat,vsdat,qvsdat,thl
     :   ,zthdat,zusdat,zvsdat,zqvsdat,psdat,ptdat,pressdat
     :   ,tedat,dz,th0,h3,u0,v0,sh3,qv0,h3c,h3e,ug_bar,ht,mom,wq3,
     :   t,ro,qc0,qr0,dpdy,vgeos,qci0,qsn0,condensat,sublim,w,
     :   vat,rad,rfl,vaq,vaqc,vau,vav,wq3_c,wq3e,ri
      real(kind(0.d0)),allocatable,dimension(:,:)::
     : u,v,th,qv,p
     
      real(kind(0.d0)),allocatable,dimension(:,:)::
     : qc,qr,qci,qsn
      real(kind(0.d0)),allocatable,dimension(:)::
     : difunu,difunv,difunt,def13,def23,difk,dift,z,
     : difk2,dift2,difunqv,difunqc,difunqr,difunqci,difunqsn,
     : def13c,def23c,dift3,difk3
      real(kind(0.d0)) wth_h,wth_h2,we
!-------------MICROPHYSICS VARS AND CONSTANTS------------!   
      real(kind(0.d0)):: prec
      real(kind(0.d0)),allocatable,dimension(:)::
     : cond, evap, auto, col, divrain, vrain,wq3c,wq3r,
     : wq3ci,wq3sn,divsnow
      real(kind(0.d0)):: rv,hlatcp,xk4,qco,xnor,rhol,xk1,dv,
     :  xkt,xniu,rho0,xa,xb,sch,gm48,gm38,gm29
      parameter(rv=461.51,hlatcp=hlat/cp,xk4=5.e-4,qco=2e-4,
     : xnor=0.8e7,nso=3.e6,dsnow=0.25,csnow=11.72, brain=0.8,
     : Thom=233.15, Tfrz=273.15)           
      parameter(rhol=1.e3,xk1=1.506,dv=2.26e-5,xkt=2.43e-2
     : ,xniu=1.51e-5,rho0=1.23,xa=842.,xb=0.8,sch=xniu/dv
     : ,gm48=17.837870,gm38=4.694174,gm29=1.827355)
      real(kind(0.d0)),allocatable,dimension(:)::
     : qs,vini,sberci,sbercw,ibercw,iacr,raci,saci,sagg,
     : sacrw,sacrwr,hmfrz,smlt,vsnow, imlt,sacrr,vdepi,
     : vdeps
      
      end module alloc_1d
