      subroutine thermo
      use alloc_1d
      implicit none
      integer iz

      do iz=1,nz-1
         vat(iz)=w(iz)*(th(iz+1,1)-
     :   th(iz,1))/
     :   (0.5*(dz(iz)+dz(iz+1)))
 !        write(0,*)vat(iz),w(iz)
      enddo
         vat(nz)=0.
         rad(nz)=0. 
!         vat=0.
!         rad=0.
      
      if(qif.ne.0) then
        do iz=1,nz
          th(iz,3)=th(iz,1)+dtl*(difunt(iz)
     :    -vat(iz)+rad(iz))
 !         write(0,*) z(iz),w(iz)
        enddo
      else
        do iz=1,nz
          th(iz,3)=th(iz,1)+dtl*(difunt(iz)-vat(iz)+rad(iz))
        enddo
      endif
!      th(nz,3)=th0(nz)      
      end
