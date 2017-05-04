      subroutine bl_depth(nonloc_tm86,difloc,nonloc_lock)
      use alloc_1d
      implicit none
      logical nonloc_tm86,difloc,nonloc_lock
      real*8 dthdz2,dthdz1,hbl1,hbl2,hbl3,th_h,th_m,hbl0,maxdt
      real*8 ws_m,wstar,thv,w_m3,thl_ct,hbl4
      integer iz,mz,i,nlct
      real, parameter:: b_0 = 6.5
      real, parameter:: b_t = 46.
      real,parameter:: B=5.
      real, parameter:: A=4.5
      real, parameter:: karm = 0.4
      real*8 dthdz(1:nz)
      real surf_flux
      dthdz=0.
      
      nlct=0
      
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
c----------------------------------------------------------
!      do iz=2,nz    !CRASHES WHEN THERE IS NO MOISTURE
!         if(th(iz,2)-hlat/cp*qc(iz,2).gt.th(1,2)+1.) then
!            hbl4=z(iz)-0.49*dz(iz)
!            goto 20
!         endif
!      enddo
 
      
! 20   do iz =3,nz-1
!         dthdz2= ((th(iz+1,2)-hlat/cp*qc(iz+1,2))
!     :  -(th(iz,2)-hlat/cp*qc(iz,2)))/(z(iz+1)-z(iz))
!      if(dthdz2.gt.0.003) then
!         write(0,*) 'dth',th(iz+1,1)-th(iz,1)
!     :   ,(th(iz+1,2)-hlat/cp*qc(iz+1,2))
!     :    -(th(iz,2)-hlat/cp*qc(iz,2))
!      hbl1 = max(50.,z(iz+1)-(dthdz2-0.003)/
!     :	(dthdz2)
!     :    *(z(iz+1)-z(iz)))
!      mz=int(iz/2.)
c------------------------------------------------------------      
 !     hbl1=z(iz)+0.1           !easy option, but not smooth
 !     goto 30
 !     endif
 !     enddo
!     30    continue
c----------------------------------------------------------      
      if (seaice.eq.1.and.frac.lt.1) then
         surf_flux = frac*ust_s*tst_s + (1.-frac)*ust_s2*tst_s2
      else
         surf_flux = ust_s*tst_s
      endif
      if (surf_flux.lt.0) then
         do iz=1,nz-1
            if(h3(iz)+h3c(iz)+h3e(iz).ge.0) then
               hbl2=z(iz)
               goto 40
            endif
         enddo
      else
         do iz=1,nz-1
            if(sqrt(def13(iz)**2.+def23(iz)**2.)
     :         .le.0.05*(ust_s**2.))then
               hbl2=z(iz)/0.95
               goto 40
            endif
         enddo
      endif
      
 40   continue
      
!      do i=1,3
!      hbl=hbl3
!      if(i.eq.1)hbl=hbl1
      hbl = hbl2
      
      wstar=(g/thv*Fv*hbl)**(1./3.)
      ws_m=(ust_s**3.+7.*karm*wstar**3.*0.5)**(1./3.)
      w_m3=wstar**3.+B*ust_s**3.
      !wth_h=-A*w_m3/hbl
      th_m=b_t*abs(wth_h)/ws_m !b_0*abs(Fv)/ws_m !
!      th_h=th(mz,2)+th_m !th(1,2)+th_m !
      
!      do iz = mz,nz-1
!      if(th(iz,2)-hlat/cp*qc(iz,2).lt.th_h
!     : .and.th(iz+1,2)-hlat/cp*qc(iz+1,2).ge.th_h) then
!      hbl3=z(iz)+(z(iz+1)-z(iz))*(th_h-th(iz,2)+hlat/cp*qc(iz,2))
!     : /(th(iz+1,2)-hlat/cp*qc(iz+1,2)-th(iz,2)+hlat/cp*qc(iz,2))
!      endif
!      enddo
      
!      enddo

!     diagnostics of cloud top height and 
!     bottom level of top-down diffusion
      if (ifwr.ne.0) then
      if (nonloc_lock) then
         do iz = 1,nz
            if (qc(iz,2).gt.0.and.qc(iz+1,2).eq.0) then
               zct=0.5*(z(iz)+z(iz+1))
               nlct=iz
            endif
         enddo
         do iz=1,nz
            thl(iz)=th(iz,2)+hlat/cp*qv(iz,2)
         enddo
         thl_ct=thl(nlct)
         do iz=nlct,1,-1
            if(thl(iz).lt.thl_ct.and.thl(iz+1).ge.thl_ct) 
     :      zc2=0.5*(z(iz)+z(iz+1))
         enddo
            
      endif
      endif

      
      
      
     
 50   hbl=hbl2                  !hbl4   !max(hbl3,hbl1)
      write(0,*) hbl
 !     write(0,*) hbl,zi_rec

 !     if (nstep.eq.1) then
 !        call reconstruction 
 !        zi1 = zi_rec
 !        zi2 = zi_rec
 !     endif
      
  !    call zi_prognostic
  !    call reconstruction
  !    hbl = zi_rec
      end

      subroutine reconstruction
      use alloc_1d
      implicit none
      real*8 a,b,c,d,deltaz,gfa,gml,sumth,thm,xnum
      real*8 th2(1:nz)
      integer kt,iz
   !   write(0,*) 'hbl=',hbl,wth_h
      sumth=0.
      xnum =0.
      do iz = 2,nz-1
         if (nstep.eq.1.and.z(iz).gt.hbl.and.z(iz-1).lt.hbl) kt = iz-1
         if (nstep.ne.1.and.0.5*(z(iz)+z(iz+1))
     :      .ge.hbl.and.0.5*(z(iz-1)+z(iz)).le.hbl) then
           kt = iz-1
!           WRITE(0,*) 'kt = ',kt
         endif
!         if (z(iz).gt.hbl.and.z(iz-1).lt.hbl) kt = iz-1
         th2(iz)=th(iz,2)-hlat/cp*qc(iz,2)*th(iz,2)/t(iz)
         if (z(iz).lt.hbl) then
            sumth = sumth+th2(iz)
            xnum=xnum+1
         endif
      enddo
      thm=sumth/xnum
      gfa = (th2(kt+3)-th2(kt+2))/dz(kt+2)
      gml = (th2(kt) - th2(kt-1))/dz(kt)
      if (nstep.eq.1.and.gfa.eq.0.and.gml.eq.0) then
         zi_rec = (0.5*(z(kt+1)+z(kt+2))*(th2(kt+2)-th2(kt+1))+
     :            0.5*(z(kt)+z(kt+1))*(th2(kt+1)-th2(kt)))/
     :             (th2(kt+2)-th2(kt))
      else
         a = 0.5*(gfa - gml)
         b = - (th2(kt+2)-gfa*(z(kt+2)-0.5*(z(kt+1)+z(kt+2))))+
     :    (th2(kt)+gml*(0.5*(z(kt+1)+z(kt+2))-z(kt)))
         c = (0.5*(z(kt+1)+z(kt+2))-0.5*(z(kt)+z(kt+1)))*
     :    (th2(kt+1) - (th2(kt) + gml*(z(kt+1) - z(kt))))
         d = b**2. - 4.*a*c
         deltaz = 0.5*(-b-sqrt(d))/a
         if(nstep.eq.1) then
            zi_rec = 0.5*(z(kt)+z(kt+1))+0.1
        endif
      endif
!      write(0,*) 'zi_rec =',zi_rec,z(kt),z(kt+1),z(kt+2)
!      write(0,*) 'Discriminant =', d
!      if(nstep.eq.145) stop
      
  !    write(0,*) th2(kt),th2(kt+1),th2(kt+2)
      
      thzi = th(kt,2)+(zi_rec - z(kt))*(th(kt,2)-th(kt-1,2))
     :      /(dz(kt))
      qczi = qc(kt,2)+(zi_rec - z(kt))*(qc(kt,2)-qc(kt-1,2))
     :      /dz(kt)
      qvzi = qv(kt,2)+(zi_rec - z(kt))*(qv(kt,2)-qv(kt-1,2))
     :      /dz(kt)
      delta_th = th(kt+2,2)-(z(kt+2)-zi_rec)*(th(kt+3,2)-th(kt+2,2))
     :      /dz(kt+2) - thzi
      delta_qc = qc(kt+2,2)-(z(kt+2)-zi_rec)*(qc(kt+3,2)-qc(kt+2,2))
     :      /dz(kt+2) - qczi
      delta_qv = qv(kt+2,2)-(z(kt+2)-zi_rec)*(qv(kt+3,2)-qv(kt+2,2))
     :      /dz(kt+2) - qvzi
   !   if (hbl.ge.900) then
   !         write(0,*) hbl, delta_th,delta_qc,delta_qv
           ! stop
  !       endif
      end

      subroutine zi_prognostic
      use alloc_1d
      implicit none
      real*8 w_ls,zi
      integer iz
      
      do iz = 2, nz-1
  !       if(z(iz).gt.hbl.and.z(iz-1).lt.hbl) then
         if (0.5*(z(iz)+z(iz+1))
     :      .ge.hbl.and.0.5*(z(iz-1)+z(iz)).lt.hbl) then
            w_ls =  0.5*(w(iz)+w(iz+1))
 !           delta_th = 9.
 !           delta_qv = 8.e-3
 !           delta_qc = 1.e-3
            w_e = wth_h/(delta_th+0.61*th(iz,2)*delta_qv)
            !w_e = -1.5e-3
         endif
      enddo
      ! hbl = zi_old
      ! if (Fv.gt.0) then
!      write(0,*) 'w_e+w_ls', -w_e+w_ls
         zi = zi1 + dtl*(-w_e + w_ls)
         
         hbl = zi1
         zi_rec=zi1
         zi1 = zi2
         zi2 = zi
      ! else
      !   zi = hbl
      ! endif

      end
