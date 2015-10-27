      subroutine humid
      use alloc_1d
      implicit none
      integer iz

!      vaq=0.
!      vaqc=0.
      if (ifwr.ne.0) then
        do iz=1,nz
          qv(iz,3)=qv(iz,1)+dtl*
     :     (difunqv(iz) -vaq(iz)) !-condensat(iz)) !-vaq(iz))
          qc(iz,3)=qc(iz,1) +dtl*
     :     (difunqc(iz) -vaqc(iz)) !+condensat(iz)) !-vaqc(iz))
          qr(iz,3)=qr(iz,1)+dtl*
     :     (difunqr(iz))
        if (ifmf.ne.0) then
          qci(iz,3)=qci(iz,1)+dtl*
     :     (difunqci(iz))          
          qsn(iz,3)=qsn(iz,1)+dtl*
     :     (difunqsn(iz))
        endif
        enddo
      else
        do iz=1,nz
          qv(iz,3)=qv(iz,1)+dtl*difunqv(iz)
        enddo
      endif

!      qv(nz,3)=qv0(nz)
!      qc(nz,3)=qc0(nz)
!      qr(nz,3)=qr0(nz)
!      qci(nz,3)=qci0(nz)
!      qsn(nz,3)=qsn0(nz)
      
      end
