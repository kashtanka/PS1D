      subroutine const_fluxes
      use alloc_1d
      implicit none
      real*8 uvs

!-------DYCOMSII surface forcing------------------!
      h=15.
      le=115.
      ust_s=0.25
!-------------------------------------------------!
      uvs=sqrt(u(1,2)**2.+v(1,2)**2.)
      tst_s=-h/ro(1)/cp/ust_s
      qst_s=-le/ro(1)/hlat/ust_s
      cdm=ust_s**2/uvs
      dzits=z_sl/(ust_s**2*th(1,2))*(g*0.4*tst_s)
      Fv=-ust_s*tst_s-0.61*t(1)*ust_s*qst_s
      end
