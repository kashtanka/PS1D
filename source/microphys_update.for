      subroutine microphys_update
!-------------------------------------------------!
! prognostic variables are updated: 
! tendencies related with microphysics are added
! in frames of the time integration scheme
!-------------------------------------------------!
      use alloc_1d, only:
! prognostic variables:      
     : qv,qc,qr,qci,qsn,th,
! diagnostic variables:
     : p,
! constants
     : hlat,cp,p00,akapa,hfus,
! model parameters:
     : dtl,nz,ifwr,ifmf,
! microphysical tendencies:
     : evap,auto,col,divrain,
     : vdepi,vdeps,vini,imlt, 
     : hmfrz,sacrw,sacrwr,sbercw,
     : ibercw,smlt,sacrr,iacr,sberci,
     : sagg,saci,raci,divsnow

     
      implicit none
      integer iz
      
      if (ifwr.ne.0.and.ifmf.eq.0) then
!-------------------------------------------------!
! Tendencies related with liquid phase are added
! to prognostic variables    
!-------------------------------------------------!  
        do iz=1,nz
          qv(iz,3)=qv(iz,3)+dtl*
     :      (evap(iz)) 
          qc(iz,3)=qc(iz,3) +dtl*
     :      (-auto(iz)-col(iz))           
          qr(iz,3)=qr(iz,3)+dtl*
     :      (auto(iz)+col(iz)-evap(iz)+divrain(iz))
          th(iz,3)=th(iz,3)+dtl*(hlat/cp*(p00/p(iz,2))**akapa
     :      *(-evap(iz)))
        enddo
      endif
      
      if (ifmf.ne.0) then
!------------------------------------------------------------!
! Tendencies related with all phases are added.
!
! In the current version vdepi=0; vdeps=0.; vini=0  because
! ice depositional growth (vdepi) and
! ice initiation (vini) are taken into account in
! the saturation adjustment scheme. The corresponding subroutines are
! either absent or not invoked. Snow depositional growth (vdeps)
! is neglected for simplicity
!------------------------------------------------------------! 
        do iz=1,nz
          qv(iz,3)=qv(iz,3)+dtl*
     :      (evap(iz)-vdepi(iz)- 
     :      vdeps(iz)-vini(iz))
          qc(iz,3)=qc(iz,3) +dtl*
     :      (-auto(iz)-col(iz)+imlt(iz) 
     :      -hmfrz(iz)-sacrw(iz)-sacrwr(iz)
     :      -sbercw(iz)-ibercw(iz))
          qr(iz,3)=qr(iz,3)+dtl*
     :      (auto(iz)+col(iz)-evap(iz)+
     :      smlt(iz)-sacrr(iz)+sacrwr(iz)
     :      -iacr(iz)+divrain(iz))
          qci(iz,3)=qci(iz,3)+dtl*
     :      (vini(iz)+vdepi(iz)+hmfrz(iz)
     :      +ibercw(iz)-imlt(iz)-sberci(iz)
     :      -sagg(iz)-saci(iz)-raci(iz))
          qsn(iz,3)=qsn(iz,3)+dtl*
     :      (vdeps(iz)+sbercw(iz)+sberci(iz)+
     :      sacrw(iz)+sagg(iz)+saci(iz)+
     :      raci(iz)+iacr(iz)+sacrr(iz)-
     :      smlt(iz)+divsnow(iz))
          th(iz,3)=th(iz,3)+dtl*
     :      (hlat/cp*(p00/p(iz,2))**akapa*
     :      (-evap(iz))+hfus/cp*(p00/p(iz,2))**akapa*(
     :       hmfrz(iz)+sacrw(iz)+iacr(iz)+         
     :       sbercw(iz)+ibercw(iz)+sacrr(iz)-
     :       smlt(iz)-imlt(iz)))
        enddo
      endif
        
      end
