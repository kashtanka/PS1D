      subroutine thermo
      use alloc_1d
      implicit none
      integer iz
      real,parameter:: Kr = 2./86400.
      
      if(qif.ne.0) then
        do iz=1,nz
!           write(0,*) th(iz,3), difunt(iz)
          th(iz,3)=th(iz,1)+dtl*(difunt(iz)) !+condensat(iz)*hlatcp
 !    :               *(p00/p(iz,2))**akapa)
 !    :     -vat(iz)+rad(iz))
        enddo
      else
        do iz=1,nz
          th(iz,3)=th(iz,1) + dtl*( difunt(iz)  )
!     :       - Kr*exp(-(z(iz)-z_sl)/600.) )   ! LONGWAVE RADIATIVE COOLING AFTER VIHMA ET AL. 2003 
!-vat(iz)+rad(iz))
          !write(0,*) th(iz,3)
        enddo
      endif
      if (iftf) then
         do iz = 1,nz
            th(iz,3)=th(iz,3)+dtl*tfc(iz)
         enddo
      endif
!      th(nz,3)=th0(nz)
      end
