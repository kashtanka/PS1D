      subroutine alloc_goddard
      use alloc_1d
      use goddard_mod
      implicit none
      nsur = 2
      na = 1
c---------INPUT VARS---------c

      allocate(oa(nz))
      allocate(fs(1,nsur))
      allocate(tsurfs(1,nsur))
      allocate(eg(1,nsur,9),ev(1,nsur,9))
      allocate(rvir(1,nsur,9))
      allocate(cwc(1,nz,3))
      allocate(taucl(1,nz,3))
      allocate(taual_lw(1,nz,10,na))
      allocate(ssaal_lw(1,nz,10,na))
      allocate(asyal_lw(1,nz,10,na))
      allocate(fcld(1,nz))

c---------OUTPUT VARS--------c
      allocate(flx_lw(1,nz+1))
      allocate(flc_lw(1,nz+1))
      allocate(dfdts(1,nz+1))
      allocate(sfcem(1))
      end
