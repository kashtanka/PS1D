      subroutine satur_adjust
      use alloc_1d
      implicit none
      real*8 cond2
      real, external:: qsat,esat
      real*8 qsatur,rasrv
      integer iz
      
      
      do iz=1,nz
        rasrv=287.05/461.51
        qsatur=qsat(t(iz),p(iz,2))
        if(qv(iz,3).gt.qsatur.or.qc(iz,3).gt.0.) then
          cond2=-min(qc(iz,3)
     :    ,(qsatur-qv(iz,3))/   !hlat**2.*0.001/(cp*rv*250.**2))
     :    (1.+4097.93*hlat*qsatur/(cp*(t(iz)-35.86)**2)/(1.-
     :    (1.-rasrv)*esat(t(iz))/p(iz,2))))/dtl
        else
          cond2=0.
        endif
        qv(iz,3)=qv(iz,3)+dtl*(-cond2)
        qc(iz,3)=qc(iz,3) +dtl*cond2
      enddo
      
      end