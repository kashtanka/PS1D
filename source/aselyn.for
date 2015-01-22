      subroutine aselyn
      use alloc_1d
      implicit none
      integer iz
!------------------------ASELYN TIME FILTER-----------------------!
       do iz=1,nz
         u(iz,2)=u(iz,2)+0.1*(u(iz,1)+u(iz,3)-2.*u(iz,2))
         v(iz,2)=v(iz,2)+0.1*(v(iz,1)+v(iz,3)-2.*v(iz,2))
         th(iz,2)=th(iz,2)+0.1*(th(iz,1)+th(iz,3)-2.*th(iz,2))
         qv(iz,2)=qv(iz,2)+0.1*(qv(iz,1)+qv(iz,3)-2.*qv(iz,2))
         qc(iz,2)=qc(iz,2)+0.1*(qc(iz,1)+qc(iz,3)-2.*qc(iz,2))
         qr(iz,2)=qr(iz,2)+0.1*(qr(iz,1)+qr(iz,3)-2.*qr(iz,2))
         qci(iz,2)=qci(iz,2)+0.1*(qci(iz,1)+qci(iz,3)-2.*qci(iz,2))
         qsn(iz,2)=qsn(iz,2)+0.1*(qsn(iz,1)+qsn(iz,3)-2.*qsn(iz,2))
       enddo
      end
