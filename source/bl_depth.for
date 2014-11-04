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
      dthdz=0.
      
      nlct=0
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
     

      do iz=2,nz
         if(th(iz,2)-hlat/cp*qc(iz,2).gt.th(1,2)+0.5) then
            hbl4=z(iz)-5.
            !write(0,*) z(iz),th(iz,2),th(1,2)
            goto 20
         endif
      enddo
      
 20   do iz =3,nz-1
         dthdz2= ((th(iz+1,2)-hlat/cp*qc(iz+1,2))
     :  -(th(iz,2)-hlat/cp*qc(iz,2)))/(z(iz+1)-z(iz))
      if(dthdz2.gt.0.003) then
!         write(0,*) 'dth',th(iz+1,1)-th(iz,1)
!     :   ,(th(iz+1,2)-hlat/cp*qc(iz+1,2))
!     :    -(th(iz,2)-hlat/cp*qc(iz,2))
      hbl1 = max(50.,z(iz+1)-(dthdz2-0.003)/
     :	(dthdz2)
     :    *(z(iz+1)-z(iz)))
      mz=int(iz/2.)
      
      hbl1=z(iz)+0.1           !easy option, but not smooth
      goto 30
      endif
      enddo
30    continue
      do iz=1,nz-1
      if (tst_s*ust_s.lt.0) then
      if(h3(iz)+h3c(iz)+h3e(iz).ge.0) then
            hbl2=z(iz)
            goto 40
         endif
         else
         if(sqrt(def13(iz)**2.+def23(iz)**2.).le.0.05*(ust_s**2.))then
      hbl2=z(iz)/0.95
      goto 40
      endif
      endif
   
      enddo
      
 40   continue
      
      do i=1,3
      hbl=hbl3
      if(i.eq.1)hbl=hbl1
      
      wstar=(g/thv*Fv*hbl)**(1./3.)
      ws_m=(ust_s**3.+7.*karm*wstar**3.*0.5)**(1./3.)
      w_m3=wstar**3.+B*ust_s**3.
      wth_h=-A*w_m3/hbl
      th_m=b_t*abs(wth_h)/ws_m !b_0*abs(Fv)/ws_m !
      th_h=th(mz,2)+th_m !th(1,2)+th_m !
      
      do iz = mz,nz-1
      if(th(iz,2)-hlat/cp*qc(iz,2).lt.th_h
     : .and.th(iz+1,2)-hlat/cp*qc(iz+1,2).ge.th_h) then
      hbl3=z(iz)+(z(iz+1)-z(iz))*(th_h-th(iz,2)+hlat/cp*qc(iz,2))
     : /(th(iz+1,2)-hlat/cp*qc(iz+1,2)-th(iz,2)+hlat/cp*qc(iz,2))
      endif
      enddo
      
      enddo

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

      
      
      
     
50    hbl=hbl4   !max(hbl3,hbl1)
      write(0,*) 'hbl=',hbl1,hbl4
      end
