      subroutine meanABL
      use alloc_1d
      implicit none
      integer iz
      real sv,height,nbelow
      real*8 latentheat,sumheat,stheta,sqv,sqc,sqci,sqsn,sqr
      sv=0.
      nbelow=0
      do iz=1,nz
          !height=z(iz) !-z_sl
         if(z(iz).le.hbl) then
            sv=sv+v(iz,2)
            nbelow=nbelow+1
         endif
      enddo
c-------zatychka--------------
      nbelow = max(nbelow,1.)
c-----------------------------      
      ablv=sv/nbelow
       sv=0.
       sumheat=0.
       stheta=0.
       sqv=0.
       sqc=0.
       sqci=0.
       sqsn=0.
       nbelow=0.
       do iz=1,nz
          if(z(iz).le.hbl) then
             latentheat=ro(iz)*(hlat*(p00/p(iz,2))**akapa*
     :      (-evap(iz)+condensat(iz))+hfus*(p00/p(iz,2))**akapa*(
     :       sublim(iz)+hmfrz(iz)+sacrw(iz)+iacr(iz)+         
     :       sbercw(iz)+ibercw(iz)+sacrr(iz)-
     :       smlt(iz)-imlt(iz)))
             sumheat=sumheat+latentheat
             !sv=sv+1./ro(iz)*(p(iz,2)-p(iz,1))/dy
             nbelow=nbelow+1
             stheta=stheta+th(iz,3)
             sqv=sqv+qv(iz,2)
             sqc=sqc+qc(iz,2)
             sqr=sqr+qr(iz,2)
             sqci=sqci+qci(iz,2)
             sqsn=sqsn+qsn(iz,2)
!             write(0,*) hlat/cp*(p00/p(iz,2))**akapa*condensat(iz),
!     :       hlat/cp*(p00/p(iz,2))**akapa
          endif
       enddo
c-------zatychka--------------
      nbelow = max(nbelow,1.)
c-----------------------------
       bl_dpdy=sv/nbelow 
       cond_heat=sumheat/nbelow
       mth=stheta/nbelow
       mqv=sqv/nbelow
       mqc=sqc/nbelow
       mqr=sqr/nbelow
       mqci=sqci/nbelow
       mqsn=sqsn/nbelow
      end
