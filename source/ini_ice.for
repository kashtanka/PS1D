      subroutine ini_ice
      use ice_mod
      implicit none
      integer iz
      
      dzi = ice_h/nl   ! thickness of each level
      allocate(Ti(nl))
      allocate(fc(nl+1))
      allocate(dfc(nl))
!----------snow------------!
      if(snow_h.gt.0) then
         dzs = snow_h/nls
         allocate(Tsn(nls))
         allocate(fcs(nls+1))
         allocate(dfcs(nls))
      else
         dzs = 0.
      endif 

      Ki = 2.2     ! thermal conductivity of ice (Wm-1K-1)
      Ks = 0.21    ! thermal conductivity of snow (Wm-1K-1)
      roi = 916.   ! sea ice density (kgm-3)
      ci = 2100.   ! sea ice heat capacity
      rosn = 290.    ! snow density (kgm-3)
      Tb = 271.15  ! sea ice bottom temperature

      do iz = 1, nl
         Ti(iz) = 245. + iz*dzi*(Tb-245.)/ice_h
      enddo

      if (snow_h.gt.0) then
         do iz = 1,nls
            Tsn(iz) = 245.
         enddo
      endif

      end
