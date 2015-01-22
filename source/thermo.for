      subroutine thermo
      use alloc_1d
      implicit none
      integer iz
      
      if(qif.ne.0) then
        do iz=1,nz
          th(iz,3)=th(iz,1)+dtl*(difunt(iz)
     :    -vat(iz)+rad(iz))
        enddo
      else
        do iz=1,nz
          th(iz,3)=th(iz,1)+dtl*(difunt(iz)-vat(iz)+rad(iz))
        enddo
      endif
      if (iftf) then
         do iz = 1,nz
            th(iz,3)=th(iz,3)+dtl*tfc(iz)
         enddo
      endif
!      th(nz,3)=th0(nz)
      end
