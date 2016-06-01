      subroutine output(fileturb,fileprof,filemean,filewater,fileice,
     :                  hour,
     :                  hour2,hour3,LWP)
      use alloc_1d
      use ice_mod, only:
     :      dzi,Ti,dzs,Tsn,snow_h,nl,Tsi,nls
      implicit none
      character*100 fileturb,fileprof,filemean,filewater,fileice
      integer hour
      real hour2,hour3,height,LWP
      real thvflux(1:nz)
      real xxa,tflux,xt,xp,xxb,qsflux,xBf,dqsdz,dqsdz2,dqsdz3
      real alpha
      integer iz
      real,external:: qsat,qsati

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
     : '(f10.1,f11.3,3f9.2,f6.3,f9.2,f9.4,f13.8,2f10.2,
     :     f11.2,f13.2,f9.3,f10.1,f10.2,f10.6,f10.6,f10.2,2f10.6
     :     f10.2)') !,2f10.6,f10.3,f9.3,f11.7,f10.1)')
!-----------1,2,3,4--------------------------!
     : nstep*dt/60.,distY,u(1,1),v(1,1),
!------------5,6,7---------------------------!
     : sqrt(u(1,1)**2+v(1,1)**2),ust_s,hbl,
!------------8,9,10,11,12,13-----------------!
     : -tst_s*ust_s,-qst_s*ust_s,h,le,LWP, th(1,1),
!-----------14,15,16,17,18,19,20-------------!
     : qv(1,1),zi_rec,delta_th,delta_qv,w_e,zi2,wthl_h,
!-----------21,22-----------------------------!
     : wqt_h,Tsi
     : !wth_h,wth_h2,we,the,tende
     
       if(mod(nstep*dt,60.).eq.0) write(16,
     :  '(f10.1,f11.3,f10.2,5f12.6)')
     :  nstep*dt/60.,distY,mth,mqv*1000.,mqc*1000.,mqr*1000.,
     : mqci*1000.,mqsn*1000.
       !------------------------------------------------------!
       
       !-----------writing profiles---------------------------!
       if(mod(nstep*dt,3600.).eq.0)then
       
        hour=nstep*dt/3600
  
         if(hour.le.9) then 
            write(fileturb(27:27),'(i1)')hour
            write(fileprof(25:25),'(i1)')hour
            write(filewater(26:26),'(i1)')hour
            write(fileice(20:20),'(i1)')hour
         elseif (hour.gt.9.and.hour.le.99) then
            write(fileturb(26:27),'(i2)')hour
            write(fileprof(24:25),'(i2)')hour
            write(filewater(25:26),'(i2)')hour
            write(fileice(19:20),'(i2)')hour
         elseif (hour.gt.99) then
            write(fileturb(25:27),'(i3)')hour
            write(fileprof(23:25),'(i3)')hour
            write(filewater(24:26),'(i3)')hour
            write(fileice(18:20),'(i3)')hour
         endif 

         open(20,file=fileprof) 
         open(21,file=fileturb)
         open(22,file=filewater)
         if(seaice.eq.1) open(23,file=fileice)
         do iz=0,nz
           if(iz.eq.0) height=z(0)
           if(iz.gt.0) height=z(iz)-z_sl
           if (iz.gt.2.and.iz.lt.nz-1) then
              alpha = qsat(t(iz),p(iz,2))*hlat/rv/t(iz)/t(iz)
              dqsdz = alpha*(t(iz+1) - t(iz-1))/(z(iz+1) - z(iz-1))
              dqsdz3 = alpha*(th(iz+1,3)-th(iz-1,3))/(z(iz+1) - z(iz-1))
     :                 *(p(iz,2)/p00)**akapa - alpha*g/cp
              dqsdz2 = (qsat(t(iz+1),p(iz+1,2)) -
     :                  qsat(t(iz-1),p(iz-1,2)))/(z(iz+1)-z(iz-1))
           else 
              dqsdz = 0
              dqsdz2 = 0
              dqsdz3 = 0
           endif
              write(20,'(f7.2,5f10.4,f10.6,3f14.5,f12.6,2f12.5,3f15.9,
     :                   f12.6,f12.4)')
     :        height, u(iz,3),v(iz,3),
     :        sqrt(u(iz,3)**2.+v(iz,3)**2),th(iz,3),
     :        th(iz,3)-hlatcp*(p00/p(iz,2))**akapa*qc(iz,3), qv(iz,3),
     :        difk(iz),dift(iz),ri(iz),rfl(iz),vat(iz),rad(iz),
     :        difunt(iz),
     :        dqsdz2,dqsdz3,t(iz),p(iz,2) !,sqrt((ug+dpdy(iz)/fcor)**2+vgeos(iz)**2)
!------------------------ MICROPHYSICS---------------------------------!
              write(22,'(f7.2,8f12.8,3f13.9)')
     :        height,qv(iz,3),qsat(t(iz),p(iz,2)),qsati(t(iz),p(iz,2)),
     :        qs(iz),qc(iz,2),qr(iz,2),qci(iz,2),qsn(iz,2)
     :        ,hlat/cp*(p00/p(iz,2))**akapa*condensat(iz)
     :        ,sublim(iz),difunt(iz)-
     :    hlat/cp*difunqv(iz)
         enddo
!------------------------TURBULENCE------------------------------------!
!---------buoyancy flux diagnostics after Deardorff 1975---------------!
         do iz = 2,nz
            xt = 0.5*(t(iz)+t(iz-1))
            xp = 0.5*(p(iz,2)+p(iz-1,2))
            xxb = hlat/xt/rv*qsat(xt,xp)
            xxa = 1. + hlatcp/xt*xxb
            tflux = ((xp/p00)**akapa*HF2(iz) - hlatcp*
     :              (wq3(iz)+wq3c(iz)))/xxa
            qsflux = xxb*(1/th(iz,2)*HF2(iz)-hlatcp/t(iz)*
     :              (wq3(iz)+wq3c(iz)))/xxa
            if (0.5*(qc(iz,3)+qc(iz-1,3)).gt.0) then
               thvflux(iz) = (p00/xp)**akapa*(1. + 0.61*qsat(xt,xp) -
     :             0.5*(qc(iz,3)+qc(iz-1,3)))*tflux + 
     :             0.5*(th(iz,3)+th(iz-1,3))*(1.61*qsflux + 
     :             (wq3(iz)+wq3c(iz)))
            else
               thvflux(iz) =  (1.+0.61*(0.5*(qv(iz,3)+qv(iz-1,3))))*
     :           ht(iz) - 0.61 * 0.5*(th(iz,3) + th(iz-1,3))*wq3(iz)
!     :          HF(iz) - 0.61 * 0.5*(th(iz,3) + th(iz-1,3))*wq3(iz)
            endif
         enddo
         do iz=2,nz
             xp = 0.5*(p(iz,2)+p(iz-1,2))

           xBf = (xp/p00)**akapa*
     :           (1.+0.61*0.5*(qv(iz,3)+qv(iz-1,3)-qc(iz,3)-qc(iz-1,3)))
     :  *(-h3(iz)) - 0.5*(th(iz,3)+th(iz-1,3))*(0.61*wq3(iz) - wq3c(iz))
           write(21,'(f7.2,4f10.4,f13.8,f10.4,2f10.4,2f13.8,f10.4,
     :                 2f13.8)') 
     :     0.5*(z(iz-1)+z(iz))-z_sl*0.5,
     :     ht(iz),mom(iz),def13(iz)+def13c(iz),def23(iz)+def23c(iz),
     :     wq3(iz),xBf !ht(iz)-0.61*th(iz,3)*wq3(iz)+th(iz,3)*wq3c(iz)
     :     ,HF(iz),HF2(iz),wq3(iz)+wq3c(iz),wq3c(iz),thvflux(iz),
     :     fthl(iz),fqt(iz)
         enddo
!----------------------------------------------------------------------!
!---------------------ICE----------------------------------------------!
         if (seaice.eq.1) then
         if (snow_h.gt.0) then
            do iz = 1,nls
               write(23,'(f10.3,f10.2)') iz*dzs-0.5*dzs, Tsn(iz)
            enddo
            do iz = 1,nl
               write(23,'(f10.3,f10.2)') snow_h+iz*dzi-0.5*dzi, Ti(iz)
            enddo
         else   
            do iz = 1,nl
               write(23,'(f10.3,f10.2)') iz*dzi-0.5*dzi, Ti(iz)
            enddo
         endif
         close(23)
         endif
         close(20)
         close(21)
         close(22)
      endif
      !-------------------------------------------------------------!
      
      !---------------writing averaged profiles---------------------!
!      hour2=nstep*dt/3600.
!      hour3=hour
!       if(hour2.ge.hour.and.hour2.lt.hour+0.5) then
!       do iz=1,nz
!       sh3(iz)=sh3(iz)+h3(iz)
!       enddo
!       endif
       
!       if(mod(hour2,hour+0.5).eq.0) then
!          if(hour.gt.9) write(filemean(32:33),'(i2)')hour
!          if(hour.le.9) write(filemean(33:33),'(i1)')hour
!          open(22,file=filemean)
!          do iz=1,nz
 !           write(22,'(f5.0,f10.5)')0.5*(z(iz-1)+z(iz))-z_sl*0.5
 !    :            ,-sh3(iz)/(1800./dt)
 !         enddo
 !         sh3=0.
 !         close(22)
 !      endif

      end
