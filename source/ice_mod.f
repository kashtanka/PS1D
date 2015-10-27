      module ice_mod
      
      real Ki,Ks 
      real Tsi, T_int
      real Tb
      integer nl,nls
      real ice_h, dzi, snow_h, dzs
      real, allocatable:: Ti(:),fc(:),dfc(:)
      real, allocatable:: Tsn(:),fcs(:),dfcs(:)
      real, parameter:: eps = 0.98
      real, parameter:: sig = 5.67e-8
      real roi,ci,rosn

      end
