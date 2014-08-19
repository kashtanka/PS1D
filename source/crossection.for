      subroutine crossection(d1,d2,ystep,nz,z1,hz,var,nch)
      implicit none
      integer iz
      real*8 d1,d2,ystep
      real*8 dist,z1
      integer nz,nch
      real*8 var(0:nz,1:3),hz(0:nz),wrvar(0:nz)
      
      
      dist=d2-mod(d2,5.)
      
      do iz=0,nz
        
        wrvar(iz)=((dist-d1)*var(iz,1)+(d2-dist)*var(iz,2))/(d2-d1)
      enddo
      write(nch,'(<nz>f12.5)')wrvar(0:nz)
      end
      
       subroutine crossection2(d1,d2,ystep,nz,z1,hz,var,nch)
      implicit none
      integer iz
      real*8 d1,d2,ystep
      real*8 dist,z1
      integer nz,nch
      real*8 var(0:nz),hz(0:nz),wrvar(0:nz)
      
      
      dist=d2-mod(d2,5.)
      
      do iz=0,nz
        
        wrvar(iz)=((dist-d1)*var(iz)+(d2-dist)*var(iz))/(d2-d1)
      enddo
      write(nch,'(<nz>f12.5)')wrvar(0:nz)
      end