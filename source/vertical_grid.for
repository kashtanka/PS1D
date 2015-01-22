      subroutine vertical_grid
      use alloc_1d
      integer iz,gridz
      real*8 step
      step=ztop/nz
      do iz=1,nz
        z(iz)=iz*(step)
    !    z(iz)=z(iz-1)+(10.+50.*(iz/nz))
      enddo
      z_sl=0.5*z(1)
      do iz=1,nz
         dz(iz)=z(iz)-z(iz-1)
      enddo
      end
