      subroutine humid
      use alloc_1d
      implicit none
      integer iz
       do iz=1,nz-1
         vaq(iz)=w(iz)*(qv(iz+1,1)-
     :   qv(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vaqc(iz)=w(iz)*(qc(iz+1,1)-
     :   qc(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
!         write(0,*)vat(iz),w(iz)
      enddo
!      vaq=0.
!      vaqc=0.
      if (ifwr.ne.0) then
        do iz=1,nz
          qv(iz,3)=qv(iz,1)+dtl*
     :     (difunqv(iz)-vaq(iz)) !-cond(iz)+evap(iz))
    !      write(0,*) qv(iz,3),cond(iz)
          qc(iz,3)=qc(iz,1) +dtl*
     :     (difunqc(iz)-vaqc(iz)) !+cond(iz)-auto(iz)-col(iz))
          qr(iz,3)=qr(iz,1)+dtl*
     :     (difunqr(iz)) !+auto(iz)+col(iz)-evap(iz)+divrain(iz))
        
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
