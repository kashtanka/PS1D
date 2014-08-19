      subroutine entrainment
      use alloc_1d
      implicit none
      real*8 flx_min,flx
      integer iz
      flx_min=0.
      do iz=1,nz
      flx=-h3(iz)-h3c(iz)-h3e(iz)
      if(flx.lt.flx_min) flx_min=flx
      enddo
      if(-tst_s*ust_s.gt.0) then
      entrt=abs(flx_min/-tst_s*ust_s)
      else
      entrt=0.
      endif
      
      end