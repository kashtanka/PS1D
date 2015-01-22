      subroutine output(fileturb,fileprof,filemean,filewater,hour,
     :                  hour2,hour3,LWP)
      use alloc_1d
      implicit none
      character*100 fileturb,fileprof,filemean,filewater
      integer hour
      real hour2,hour3,height,LWP
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
     :     f11.2,f13.2,f9.3)') !,2f10.6,f10.3,f9.3,f11.7)')
     : nstep*dt/60.,distY,u(1,1),v(1,1),
     : sqrt(u(1,1)**2+v(1,1)**2),ust_s,
     : hbl,-tst_s*ust_s,-qst_s*ust_s,h,le,LWP, th(1,1),qv(1,1)
     : !wth_h,wth_h2,we,the,tende
     
       if(mod(nstep*dt,60.).eq.0) write(16,
     :  '(f10.1,f11.3,f10.2,5f12.6)')
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
           
              write(20,'(f5.0,4 f10.2,f10.6,3f10.5,f12.6)')
     :        height, u(iz,3),v(iz,3),
     :        sqrt(u(iz,3)**2.+v(iz,3)**2),th(iz,3),qv(iz,3),
     :        difk(iz),dift(iz),ri(iz),rfl(iz) !,sqrt((ug+dpdy(iz)/fcor)**2+vgeos(iz)**2)
     
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

      end
