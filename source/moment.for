      subroutine moment
      use alloc_1d
      implicit none
      integer iz
      real dpdz0,dpdz
      dpdz0=1.5e-3
      dpdy=0.

      do iz=1,nz

      u(iz,3)=u(iz,1)+dtl*(fcor*(v(iz,2)-vg)+difunu(iz)
     :        -vau(iz)) 

 !-vgeos(iz))+difunu(iz)) 
!      t=th(iz,2)*(p(iz,2)/p00)**akapa	
!      ro=p(iz,2)/(r*t)      
!      if (dt*nstep.gt.3600.*icetime
!     : .and.z(iz).le.hbl+200) then
!      dpdz=(dpdz0*(hbl-z(iz))/hbl)
!      dpdy(iz)=1./ro(iz)*(p(iz,2)-p(iz,1))/dy
      !write(0,*) z(iz),ug-1./ro/fcor*dpdz, ug_bar(iz)
 !     else
      dpdz=0.
 !     endif
            
      v(iz,3)=v(iz,1)+dtl*(-fcor*(u(iz,2)-ug)+difunv(iz)
     :        -vav(iz))
!     :  -dpdy(iz)) !-1./ro*dpdz)  !  !-1./ro*dpdz)
      !write(0,*)iz,ug-1./fcor*dpdy
      
      enddo
      u(nz,3)=ug
      v(nz,3)=vg !vgeos(nz)  !vg      
     
      end
