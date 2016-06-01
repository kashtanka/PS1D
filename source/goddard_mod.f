      module goddard_mod
!---------input vars-------------!
      real ta_s
      real co2,n2o,ch4,cfc11,cfc12,cfc22
      logical vege,high,trace,overcast_lw,cldwater_lw,aerosol
      integer nsur,ict,icb,na
      real,allocatable:: oa(:),fs(:,:),tsurfs(:,:),eg(:,:,:),ev(:,:,:)
      real,allocatable:: rvir(:,:,:),cwc(:,:,:),taucl(:,:,:)
      real,allocatable:: taual_lw(:,:,:,:),ssaal_lw(:,:,:,:)
      real,allocatable:: asyal_lw(:,:,:,:)
      real, allocatable:: fcld(:,:)
!--------output vars-------------!
      real,allocatable:: flx_lw(:,:),flc_lw(:,:),dfdts(:,:),sfcem(:)
      end
