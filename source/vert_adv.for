      subroutine vert_adv
      use alloc_1d
      implicit none
      integer iz

      do iz=1,nz-1                        ! vertical advection
         vat(iz)=w(iz)*(th(iz+1,1)-
     :   th(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vau(iz)=w(iz)*(u(iz+1,1)-
     :   u(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vav(iz)=w(iz)*(v(iz+1,1)-
     :   v(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vaq(iz)=w(iz)*(qv(iz+1,1)-
     :   qv(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vaqc(iz)=w(iz)*(qc(iz+1,1)-
     :   qc(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
      enddo
      vat(nz)=0.
      end

