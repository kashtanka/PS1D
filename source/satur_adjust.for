      subroutine satur_adjust
      use alloc_1d, only:
! prongnostic variables:
     : qv,qc,qci,th,
! diagnostic variable:
     : t,p,qs,
! model parameters:
     : nz,dz,dtl,ifmf,ifwr,hbl,
! constants:
     : hlat,hsub,cp,p00,Thom,Tfrz,akapa,
     : condensat,sublim,z
      implicit none
      real*8 cond2,subl2
      real, external:: qsat,esat,qsati
      real*8 qsatur,rasrv,temp
      real*8 DEP,CND,A1,A2,A3,r1,r2,dq,xMo,xNice
      integer iz
      cond2=0.
      subl2=0.
      xMo=1.e-12
!--------------------------------------------------------------------!
! In the warm rain scheme oversaturation/subsaturation is removed in one step 
! following the Rutledge and Hobbs (1983)
!--------------------------------------------------------------------!
      if (ifwr.ne.0.and.ifmf.eq.0) then
        do iz=1,nz
          rasrv=287.05/461.51
         ! temp=th(iz,2)*(p00/p(iz,2))**akapa
          qsatur=qsat(t(iz),p(iz,2))
          if(qv(iz,3).gt.qsatur.or.qc(iz,3).gt.0.) then
            cond2=-min(qc(iz,3)
     :      ,(qsatur-qv(iz,3))/  
     :      (1.+4097.93*hlat*qsatur/(cp*(t(iz)-35.86)**2)/(1.-
     :      (1.-rasrv)*esat(t(iz))/p(iz,2))))/dtl
          else
            cond2=0.
          endif
!          if (z(iz).gt.hbl) then
!             cond2=0.
!             endif
          qv(iz,3)=qv(iz,3)+dtl*(-cond2)
          qc(iz,3)=qc(iz,3) +dtl*cond2
          th(iz,3)=th(iz,3)+dtl*hlat/cp !*(p00/p(iz,2))**akapa
     :    *(cond2)
 !         write(0,*)z(iz),hlat/cp*cond2
 !         write(0,*) z(iz),p(iz,2)
          condensat(iz)=cond2
        enddo
      endif
!--------------------------------------------------------------------!
! In the mixed-phase scheme the super- / subsaturation is removed
! following the Tao et al. (1989) scheme
!--------------------------------------------------------------------!
      if (ifmf.ne.0) then
        do iz = 1,nz
          
	    if (t(iz).gt.Thom.and.t(iz).le.Tfrz) then
            CND=(t(iz)-Thom)/(Tfrz-Thom)
	      DEP=(Tfrz-t(iz))/(Tfrz-Thom)
	     elseif (t(iz).gt.Tfrz) then
	      CND=1.
	      DEP=0.
	     elseif (t(iz).le.Thom) then
	      CND=0.
	      DEP=1.
	     endif
	     A1=(237.3*17.27*(p(iz,2)/p00)**akapa)/(t(iz)-35.5)**2
	     A2=(265.5*21.88*(p(iz,2)/p00)**akapa)/(t(iz)-7.5)**2
	     A3=(hlat*CND+hsub*DEP)/(cp*(p(iz,2)/p00)**akapa)
	     r1=qv(iz,3)-qs(iz) 
	     if (t(iz).le.Tfrz) then
	       xNice=min(10.**5,0.01*exp(0.6*(273.15-t(iz))))
	       if(qci(iz,3).gt.0.or.qc(iz,3).gt.0.) then
	       r2=(A1*qc(iz,3)*qsat(t(iz),p(iz,2))+
     :	      A2*qci(iz,3)*qsati(t(iz),p(iz,2)))
     :        /(qc(iz,3)+qci(iz,3))
              else
              r2=CND*A1*qsat(t(iz),p(iz,2))+
     :        DEP*A2*qsati(t(iz),p(iz,2))         
              endif
           else
             r2=A1*qs(iz)
           endif
           dq=-r1/(1+r2*A3)
           if(r1.gt.0) then
             cond2=-dq*CND/dtl
             subl2=-dq*DEP/dtl
           else
             if (qc(iz,3).gt.0) then
             cond2=max(-dq*CND,-qc(iz,3))/dtl
             else
             cond2=0
             endif
             if (qci(iz,3).gt.0) then
             subl2=max(-dq*DEP,-qci(iz,3))/dtl
             else
             subl2=0
             endif
           endif
           qv(iz,3)=qv(iz,3)+dtl*(-cond2-subl2)
           qc(iz,3)=qc(iz,3)+dtl*cond2
           qci(iz,3)=qci(iz,3)+dtl*subl2
           th(iz,3)=th(iz,3)+dtl*(hlat/cp*(p00/p(iz,2))**akapa
     :     *(cond2)+hsub/cp*(p00/p(iz,2))**akapa
     :     *(subl2))
           condensat(iz)=cond2
           sublim(iz)=subl2
          
           
	     
        enddo
      endif
      
      end
