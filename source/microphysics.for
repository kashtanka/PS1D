      subroutine microphysics
!-------------------------------------!
! Here all the tendencies related with 
! cloud microphysics are calculated.
! Saturation adjustment is done later.
!-------------------------------------!
      use alloc_1d, only:
! model parameters:
     : ifwr,ifmf
      implicit none
      
      if (ifwr.ne.0) then
!-------------------------------------!
! tendencies related with liquid phase microphysics
! are calucalted
!-------------------------------------!
!        call evaporate
!        call convert
!        call collect
!        call rain
      endif
      
      if(ifmf.ne.0) then
!-------------------------------------!
! tendencies related with mixed/solid phase 
! microphysics are calucalted.
! Inclusion of Bergeron process
! described in Lin et al 1983 is under question
!-------------------------------------!
        call saggregation
        call accretcwsn
!        call bergeron
        call accretcisn
        call accret3comp
        call homofreeze
        call accretsr
        call snowmelt
        call icemelt
        call snowprec
      endif  
      end
      
      
      subroutine evaporate
!
! Evaluation of the evaporation of rainwater (rain water<->vapour)
!
! evap = a rather complex expression
!
! after Emanuel (1994) and Ogura and Takahashi (1973)
!
      use alloc_1d, only:
! prognostic variables:
     : qv,qr,
! diagnostic variables:
     : t,p,ro,
! tendencies:
     : evap,
! model parameters:
     : nz,
! constants:
     : pi,dv,xnor,hlat,xkt,sch,
     : gm29,xa,xniu,rho0,xb,rhol,
     : rv
      implicit none
      real*8 qsatur,xlr
      real, external:: qsat
      integer iz
      
      do iz=1,nz
            if(qr(iz,2).gt.0) then
            xlr=(ro(iz)*qr(iz,2)/(pi*rhol*xnor))**0.25
	        qsatur=qsat(t(iz),p(iz,2))
	      evap(iz)=2.*pi*dv*xnor*(qsatur-qv(iz,2))
     :      /(1.+dv*hlat*hlat*qsatur*ro(iz)/(xkt*rv*t(iz)**2.))*
     :        (0.78*xlr*xlr+0.31*sch**(1./3.)*gm29*sqrt(xa/xniu)*
     :        (rho0/ro(iz))**0.25*xlr**((xb+5.)/2.))
        else
          evap(iz)=0.
        endif
      enddo
      
      
      return
      end
      
      subroutine condens
c
c Evaluation of the condensation/evaporation (vapour<->cloud water)
c
c cond=d qv/dt=(qsat-qv)/(dt*(1+lv**2*qsat/cp*r*t))
c
c after Yau and Austin (1979), Wilhelmson and Klemp (1978)
c
      use alloc_1d, only:
! prongnostic variables:
     : qv,qc,
! diagnostic variables:
     : t,p,
! tendencies:
     : cond,
! model parameters:
     : nz,dtl,
! constants:
     : hlat,cp
      implicit none
      real, external:: qsat,esat
      real*8 qsatur,rasrv
      integer iz

      do iz=1,nz
            rasrv=287.05/461.51
            qsatur=qsat(t(iz),p(iz,2))
            if(qv(iz,2).gt.qsatur.or.qc(iz,2).gt.0.) then
              cond(iz)=-min(qc(iz,2)
     :          ,(qsatur-qv(iz,2))/
     :          (1.+4097.93*hlat*qsatur/(cp*(t(iz)-35.86)**2)/(1.-
     :          (1.-rasrv)*esat(t(iz))/p(iz,2))))/dtl
            else
              cond(iz)=0.
            endif
      enddo
      return
      end
      
      
      
      subroutine convert
c
c Evaluation of the autoconversion term (cloud water->rain water)
c
c auto=k1(qc-qc0)
c
c after Kessler(1969)
c
      use alloc_1d, only:
! prongnostic variables:
     : qc,
! tendencies:
     : auto,
! model parameters:
     : nz,
! constants:
     : qco,xk4
      implicit none
      integer iz
      
      do iz=1,nz
        if (qc(iz,2).gt.qco) then
          auto(iz)=xk4*(qc(iz,2)-qco)
        else
          auto(iz)=0.
        endif
      enddo
      
      return
      end
      
      subroutine collect
c
c Evaluation of the collection of raindrops by cloud (cloud water->rain water)
c
c col=const*rho*(3/8)*qc*qr**(7/8)
c
c after Emanuel (1994)
c
      use alloc_1d, only:
! prongnostic variables:
     : qr,qc,
! diagnostic variables:
     : ro,
! tendencies:
     : col,
! model parameters:
     : nz,
! constants:
     : pi,rhol,xnor,xa,gm38,xb,rho0

      implicit none
      real*8 xlr
      integer iz

      do iz=1,nz
          if(qr(iz,2).gt.0.and.qc(iz,2).gt.0) then
            xlr=(ro(iz)*qr(iz,2)/(pi*rhol*xnor))**0.25
	      col(iz)=pi*xnor*xa*qc(iz,2)*gm38*xlr**(xb+3)
     :        *0.25*sqrt(rho0/ro(iz))
          
          else
            col(iz)=0.
          endif
      enddo
      return
      end
      
      subroutine rain
c
c Evaluation of the liquid precipitation velocity (rain water->soil)
c vrain=const*rho**(-3/8)*qr**(1/8)
c after Emanuel(1994)
c
      use alloc_1d, only:
! prongnostic variables:
     : qr,
! diagnostic variables:
     : ro,vrain,
! tendencies:
     : divrain,
! model parameters:
     : nz,dz,
! constants:
     : pi,rhol,xnor,xa,gm48,xb,rho0

      implicit none
      real*8 xlr
      integer iz

      do iz=1,nz
            xlr=(ro(iz)*qr(iz,2)/(pi*rhol*xnor))**0.25
            vrain(iz)=xa*gm48/6.*xlr**xb*sqrt(rho0/ro(iz))
      enddo
      do iz=2,nz-1
            divrain(iz)=(ro(iz+1)*vrain(iz+1)*qr(iz+1,2)
     :       -ro(iz-1)*vrain(iz-1)*qr(iz-1,2))
     :       /(2*dz(iz))/ro(iz)
      enddo
          divrain(1)=divrain(2)
      return
      end
      
      
      subroutine iceini
c--------------------------------------------------------------------------------------------------------------c      
c
c      initiation of ice crystals   (water vapor => cloud ice)                            DC,10.2009
c      after Dudhia,1989
c---------------------------------------------------------------------------------------------------------------c
	use alloc_1d, only:
! prongnostic variables:
     : qci,qv,
! diagnostic variables:
     : t,p,
! tendencies:
     : vini,
! model parameters:
     : nz,dtl   
      
      implicit none
	integer iz
	real(8) xMo,xNice
	real(8),external::qsati,qsat
	
	xMo=1.e-12
      
      do iz=1,nz
	  if (t(iz).le.273.15) then
	    xNice=min(10.**5,0.01*exp(0.6*(273.15-t(iz))))                 !ice crystals concentration after Fletcher
	  else
	    xNice=0.
	  endif
	  if(t(iz).le.273.15) then
	    vini(iz)=max((xMo*xNice-qci(iz,2))/dtl,0.)
	  else
	    vini(iz)=0.
	  endif
		if (qv(iz,2)-qsati(t(iz),p(iz,2)).gt.0) then
          vini(iz)=min(vini(iz),
     :	  (qv(iz,3)-qsati(t(iz),p(iz,2)))/dtl)
	  else
	    vini(iz)=0.
	  endif
	enddo
      return
	end
	
	subroutine icemelt
c-----------------------------------------------------------------------------------------------------------------c
c     cloud ice melting rate if T>0 Celcium  (cloud ice => cloud water)
c
c-----------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qci,
! diagnostic variables:
     : t,
! tendencies:
     : imlt,
! model parameters:
     : nz,dtl   
	implicit none
	integer iz
      do iz=1,nz
	  if (t(iz).gt.273.15) then
	    imlt(iz)=max(0.,qci(iz,2))/dtl
	  else
	    imlt(iz)=0.
	  endif
	enddo
	return
	end
	
	subroutine homofreeze
c-------------------------------------------------------------------------------------------------------------------c
c      homogeneous freezing of cloud water at T < Too = -40. Celcium
c      cloud water => cloud ice
c-------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qc,
! diagnostic variables:
     : t,
! tendencies:
     : hmfrz,
! model parameters:
     : nz,dtl   
	implicit none
	integer iz
	real*8 Too
	
	Too=233.15    
      do iz=1,nz
	  if(t(iz).lt.Too) then
	    hmfrz(iz)=max(0.,qc(iz,2))/dtl
	  else
          hmfrz(iz)=0.
	  endif
	enddo
	return
	end
	
	subroutine accretcwsn
c------------------------------------------------------------------------------------------------------------------c
c       calculation of accretion rate of cloud water by snow  (after Lin et al,1983)
c       (cloud water => snow)
c       if t>0 then it contributes to rain water growth 
c------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qc,qsn,
! diagnostic variables:
     : ro,t,
! tendencies:
     : sacrw, sacrwr,
! model parameters:
     : nz,dtl,
! constants:
     : dsnow,rho0,nso,pi,csnow 
	implicit none
	integer iz
	real*8 gamma
	real(8),external::lamsnow,gamf
	gamma=gamf(3.+dsnow)
	do iz=1,nz
	  if (qsn(iz,2).gt.0) then
	    sacrw(iz)=pi*nso*csnow*qc(iz,2)*gamma
     :	  *sqrt(rho0/ro(iz))/(4.*lamsnow(ro(iz),qsn(iz,2))**(3.+dsnow))
          sacrw(iz)=min(sacrw(iz),qc(iz,2)/dtl)
	  else
		  sacrw(iz)=0. 
	  endif
	  if(t(iz).ge.273.15.and.qc(iz,2).gt.0.) then
          sacrwr(iz)=min(sacrw(iz),qc(iz,2)/dtl)
	  else
          sacrwr(iz)=0.
	  endif
	enddo
	return
	end
	
	subroutine bergeron
c--------------------------------------------------------------------------------------------------------------------c
c        bergeron process (deposition and riming)  after Lin et al, 1983, Hsie et al, 1980 
c        transfer of cloud water to form snow     (sbercw)
c        and transfer from cloud ice to snow      (sberci)
c        and transfer from cloud water to cloud ice (ibercw)
c        only occurs when qsati<qv<qsat
c---------------------------------------------------------------------------------------------------------------------c

      use alloc_1d, only:
! prongnostic variables:
     : qv,qc,qsn,qci,
! diagnostic variables:
     : t,p,ro,
! tendencies:
     : sbercw,sberci,ibercw,
! model parameters:
     : nz,dtl,
! constants:
     : pi

      implicit none
	real*8 mi50,ri50,ui50,ni50,qi50,dt_1,Ein,mn,
     : xNice,mi40,ri60,mi60,ni60,qi60,ui60,mn2
      real(8), external :: a1koenig,a2koenig,qsat,qsati
	integer iz

	mi50=4.8e-10
	mi40=2.46e-10
	ri50=50.e-6
	ri60=60.e-6
	ui50=1.
	ui60=1.
      Ein=0.5
	mn=1.05e-18
	mi60=4.52e-10
	mi40=1.34e-10
      
	do iz=1,nz
	  if(qv(iz,2).ge.qsati(t(iz),p(iz,2))
     :  .and.qv(iz,2).lt.qsat(t(iz),p(iz,2)).and.
     :   t(iz).lt.273.15) then  
	     dt_1=1/(a1koenig(t(iz))*(1.-a2koenig(t(iz))))*
     :	   (mi60**(1.-a2koenig(t(iz)))-mi40**(1.-a2koenig(t(iz))))
	     qi60=qci(iz,2)*dtl/dt_1
	     ni60=qi60/mi60
!	     sbercw(iz)=ni60*(a1koenig(t(iz))*mi60**a2koenig(t(iz))+
!     :     pi*Ein*ro(iz)*qc(iz,2)*ri60**2.*ui60)
           !sbercw(iz)=min(sbercw(iz),qc(iz,2)/dtl)		      
	     sberci(iz)=qci(iz,2)/dt_1
	     sberci(iz)=min(qci(iz,2)/dtl,sberci(iz))
		   xNice=min(10.**5,0.01*exp(0.6*(273.15-t(iz))))
	!      xNice=1000.*exp(12.96*Sice-0.639)
          if(qci(iz,2).eq.0) then
		    mn2=mn
		  else  
	      mn2=qci(iz,2)*ro(iz)/xNice
	    endif
	    if(qc(iz,2).gt.0.and.qci(iz,2).gt.0) then
	      ibercw(iz)=xNice/(1000.*ro(iz))*a1koenig(t(iz))
     :      *mn**a2koenig(t(iz))
	      ibercw(iz)=min(ibercw(iz),(qc(iz,2)/dtl-sbercw(iz)),0.)
	    else
	      ibercw(iz)=0.
	    endif
        else
         sbercw(iz)=0.
	   sberci(iz)=0.
	   ibercw(iz)=0.
        endif
	enddo
	return
	end
	
	subroutine saggregation
c------------------------------------------------------------------------------------------------------------------c
c         ice crystal aggregation rate to form snow 
c         after Lin,1983 and Kessler,1969 
c         cloud ice => snow
c---------------------------------------- --------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qci,
! diagnostic variables:
     : t,
! tendencies:
     : sagg,
! model parameters:
     : nz


	implicit none
	real(8) alfa,qci_0,t0
	integer iz
      t0=273.15
	!qc_i0=1.e-3          !Lin,1983
	qci_0=1.2e-4        !BEM
      do iz=1,nz
	  alfa=1.e-3*exp(0.025*(t(iz)-t0))
	  sagg(iz)=alfa*(qci(iz,2)-qci_0)
	  if(t(iz).lt.t0.and.qci(iz,2).gt.0.) then
	    sagg(iz)=max(0.,sagg(iz))
	  else
	    sagg(iz)=0.
	  endif
	enddo
	return
	end
	
		subroutine accretcisn
c------------------------------------------------------------------------------------------------------------------c
c         accretion rate of cloud ice by snow
c         after Lin,1983
c         cloud ice => snow
c------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qci,qsn,
! diagnostic variables:
     : t,ro,
! tendencies:
     : saci,
! model parameters:
     : nz,dtl,
! constants:
     : rho0,nso,pi,dsnow,csnow

	implicit none
	real(8) t0,esi,gamma
	integer iz
	real(8), external :: gamf,lamsnow
	gamma=gamf(dsnow+3.)
	t0=273.15
	do iz=1,nz
	  esi=exp(0.025*(t(iz)-t0))
	  if (qsn(iz,2).gt.0) then
	    saci(iz)=pi*esi*nso*csnow*qci(iz,2)*gamma
     :	  *sqrt(rho0/ro(iz))*0.25/lamsnow(ro(iz),qsn(iz,2))**(dsnow+3.)
	  else
	    saci(iz)=0.
	  endif     
	  if(t(iz).lt.t0) then
	    saci(iz)=min(saci(iz),qci(iz,2)/dtl)
	  else
	    saci(iz)=0.
	  endif
	enddo
	return
	end
	
	subroutine accret3comp
c------------------------------------------------------------------------------------------------------------------c
c       Three-component process: cloud ice accreted by rain water to form snow
c       raci - sink term for cloud ice
c       iacr - sink term for rain
c       raci+iacr - source term for snow
c       after Lin, 1983
c
c------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qr,qci,qsn,
! diagnostic variables:
     : t,ro,
! tendencies:
     : raci,iacr,
! model parameters:
     : nz,dtl,
! constants:
     : rho0,xnor,xa,xb,pi,rhol

	implicit none

	real(8) eri,Mi,gamma1,gamma2
	integer iz
	real(8), external :: gamf,lamsnow,lamrain
	gamma1=gamf(3.+xb)
	gamma2=gamf(6.+xb)
	eri=1.
	mi=4.19e-13
      do iz=1,nz
	  if(qr(iz,2).gt.0.and.t(iz).lt.273.15) then
	    raci(iz)=pi*eri*xnor*xa*qci(iz,2)*gamma1
     :	  *sqrt(rho0/ro(iz))*0.25/lamrain(ro(iz),qr(iz,2))**(3.+xb)
	    raci(iz)=min(raci(iz),qci(iz,2)/dtl)
	    iacr(iz)=pi**2*eri*xnor*xa*qci(iz,2)*rhol*
     :	  gamma2*sqrt(rho0/ro(iz))/
     :    (24.*Mi*lamrain(ro(iz),qr(iz,2))**(6.+xb))
          iacr(iz)=min(iacr(iz),qr(iz,2)/dtl)
	  else
	    raci(iz)=0.
	    iacr(iz)=0.
	  endif
	enddo
	return
	end
	
	subroutine accretsr
c------------------------------------------------------------------------------------------------------------------c
c      accretion rate of snow at the expense of rain water 
c      after Lin et al, 1983
c------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qr,qsn,
! diagnostic variables:
     : t,ro,
! tendencies:
     : sacrr,
! model parameters:
     : nz,dtl,
! constants:
     : dsnow,brain,pi,nso,xnor,rhol

	implicit none
	integer iz
	real(8),external :: lamsnow,lamrain,Usn,Urn,gamf
	real(8) ls,lr,gamma1,gamma2
	gamma1=gamf(4.+dsnow)
	gamma2=gamf(4.+brain)
      do iz=1,nz
	  if(qr(iz,2).gt.0.and.qsn(iz,2).gt.0.and.t(iz).lt.273.15) then
	    lr=lamrain(ro(iz),qr(iz,2))
	    ls=lamsnow(ro(iz),qsn(iz,2))
	    sacrr(iz)=pi**2*nso*xnor*abs(Usn(ro(iz),qsn(iz,2),
     :	  gamma1)-Urn(ro(iz),qr(iz,2),gamma2))*(rhol/ro(iz))*
     :    (5./(lr**6.*ls)+2./(lr**5.*ls**2.)+0.5/(lr**4.*ls**3.))
          sacrr(iz)=min(sacrr(iz),qr(iz,2)/dtl)
	  else
	    sacrr(iz)=0.
	  endif
	enddo
	return
	end
	
	
	subroutine snowmelt
c------------------------------------------------------------------------------------------------------------------c
c          melting rate of snow to form rain 
c          after Lin et al, 1983
c------------------------------------------------------------------------------------------------------------------c
      use alloc_1d, only:
! prongnostic variables:
     : qsn,qv,
! diagnostic variables:
     : t,ro,p,
! tendencies:
     : smlt,sacrw,sacrr,
! model parameters:
     : nz,dtl,
! constants:
     : nso,pi,dsnow,hlat,csnow,rho0

	implicit none
	real(8) lfus,xmyu,xKa,diffus,cw,Sc,gamma
	integer iz
	real(8), external :: gamf,lamsnow,qsati
	gamma=gamf((dsnow+5.)/2.)
	lfus=3.336e5
	cw=4.187e3
	do iz=1,nz
	  diffus=8.794/(10.**5)*(t(iz)**1.81)/p(iz,2)
	  xmyu=1.496/(10.**6)*(t(iz)**1.5)/(t(iz)+120.)
	  xKa=1.414*10.**3*xmyu
	  Sc=xmyu/(ro(iz)*diffus)
	  if(qsn(iz,2).gt.0.and.t(iz).gt.273.15) then
	    smlt(iz)=-((-2.*pi/ro(iz)/lfus)*(xKa*(t(iz)-273.15)-hlat*diffus
     :	  *ro(iz)*(qsati(t(iz),p(iz,2))-qv(iz,2)))*nso*(0.78/
     :    lamsnow(ro(iz),qsn(iz,2))**2+0.31*Sc**(1./3.)
     :    *gamma*sqrt(csnow)*(rho0/ro(iz))**0.25/
     :    sqrt(xmyu/ro(iz))/lamsnow(ro(iz),qsn(iz,2))**((dsnow+5.)/2.))
     :    -cw*(t(iz)-273.15)/lfus*(sacrw(iz)+sacrr(iz)))
          smlt(iz)=min(smlt(iz),qsn(iz,2)/dtl)
	  else
	    smlt(iz)=0.
	  endif
	enddo
      return
	end
	
	subroutine snowprec
c------------------------------------------------------------------------------c
c           preciptating snow
c           after Lin et al, 1983
c------------------------------------------------------------------------------c
	use alloc_1d, only:
! prongnostic variables:
     : qsn,
! diagnostic variables:
     : ro,vsnow,
! tendencies:
     : divsnow,
! model parameters:
     : nz,dz,
! constants:
     : dsnow
	implicit none
	real(8) gamma
	real(8), external:: Usn,gamf
	integer iz
	gamma=gamf(4+dsnow)
      do iz=1,nz
	    vsnow(iz)=Usn(ro(iz),qsn(iz,2),gamma)
      enddo
      do iz=2,nz-1
            divsnow(iz)=(ro(iz+1)*vsnow(iz+1)*qsn(iz+1,2)
     :       -ro(iz-1)*vsnow(iz-1)*qsn(iz-1,2))
     :       /(2*dz(iz))/ro(iz)
      enddo
          divsnow(1)=divsnow(2)
	return
	end
	
	subroutine saturation
!--------------------------------------------------------------!
!     specific humidity of saturation qs is calculated
!     depending on the fraction of cloud water and cloud ice
!     or on the air temperature	
!--------------------------------------------------------------!
	use alloc_1d, only:
! prongnostic variables:
     : qv,qc,qci,qsn,
! diagnostic variables:
     : t,qs,p,
! model parameters:
     : nz,dz,ifmf,
! constants:
     :Tfrz,Thom
	implicit none
	integer iz
	real, external:: qsat,qsati
	real*8 xNice,xMo,CND,DEP
	if (ifmf.ne.0) then
	xMo=1.e-12
        do iz=1,nz
	    if(qc(iz,3).ne.0.or.qci(iz,3).ne.0.) then
	      qs(iz)=(qc(iz,3)*qsat(t(iz),p(iz,2))+
     :	    qci(iz,3)*qsati(t(iz),p(iz,2)))
     :      /(qc(iz,3)+qci(iz,3))
	    else
	     if(t(iz).le.273.15.and.t(iz).gt.233.15) then	
            CND=(t(iz)-Thom)/(Tfrz-Thom)
	      DEP=(Tfrz-t(iz))/(Tfrz-Thom)
	      qs(iz)=(CND*qsat(t(iz),p(iz,2))+
     :	      DEP*qsati(t(iz),p(iz,2)))/
     :        (CND+DEP)
           endif
	    endif
	    if(t(iz).gt.273.15)	qs(iz)=qsat(t(iz),p(iz,2))
	    if(t(iz).le.233.15)   qs(iz)=qsati(t(iz),p(iz,2))
        enddo
      else
        do iz=1,nz
          qs(iz)=qsat(t(iz),p(iz,2))
        enddo
      endif
	return
	end
      
      
      
      
      function qsat (t,p)
!          humidade especifica de saturacao

      implicit double precision (a-h,o-z)
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      espf=esat(t)/p
      qsat=rasrv*espf/(1.-etvq*espf)
	!qsat=380./p*exp(17.27*(t-273.15)/(t-35.86))
      end
  
	function qsati (t,p)                                                  

      implicit double precision (a-h,o-z)
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      espf=esati(t)/p
      qsati=rasrv*espf/(1.-etvq*espf)
	!qsati=380./p*exp(21.88*(t-273.15)/(t-7.66))
      end
      
      function dqsat (t,p)
      implicit double precision (a-h,o-z)
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      esat=618.78*exp(17.269*(t-273.16)/(t-35.86))
      desat=4097.93*esat/((t-35.86)*(t-35.86))
      dqsat=rasrv*desat*p/((p-etvq*esat)*(p-etvq*esat))
      end



      function esat (t)
      implicit double precision (a-h,o-z)
      esat=618.78*exp(17.269*(t-273.16)/(t-35.86))
	end

      function esati (t)
      implicit double precision (a-h,o-z)
      esati=611.2*10.**(0.43*22.46*(t-273.15)/t)
      end
      
      function a1koenig(t)
	implicit double precision (a-h,o-z)
	t0=273.15
	if(t.gt.t0) a1koenig=0.0
	if(t.le.t0.and.t.gt.t0-1.5) a1koenig=0.7939e-7
	if(t.le.t0-1.5.and.t.gt.t0-2.5) a1koenig=0.7841e-6
	if(t.le.t0-2.5.and.t.gt.t0-3.5) a1koenig=0.3369e-5
	if(t.le.t0-3.5.and.t.gt.t0-4.5) a1koenig=0.4336e-5
	if(t.le.t0-4.5.and.t.gt.t0-5.5) a1koenig=0.5285e-5
	if(t.le.t0-5.5.and.t.gt.t0-6.5) a1koenig=0.3728e-5
	if(t.le.t0-6.5.and.t.gt.t0-7.5) a1koenig=0.1852e-5
	if(t.le.t0-7.5.and.t.gt.t0-8.5) a1koenig=0.2991e-6
	if(t.le.t0-8.5.and.t.gt.t0-9.5) a1koenig=0.4248e-6
	if(t.le.t0-9.5.and.t.gt.t0-10.5) a1koenig=0.7434e-6
	if(t.le.t0-10.5.and.t.gt.t0-11.5) a1koenig=0.1812e-5
	if(t.le.t0-11.5.and.t.gt.t0-12.5) a1koenig=0.4394e-5
	if(t.le.t0-12.5.and.t.gt.t0-13.5) a1koenig=0.9145e-5
	if(t.le.t0-13.5.and.t.gt.t0-14.5) a1koenig=0.1725e-6
	if(t.le.t0-14.5.and.t.gt.t0-15.5) a1koenig=0.3348e-4
	if(t.le.t0-15.5.and.t.gt.t0-16.5) a1koenig=0.1725e-4
	if(t.le.t0-16.5.and.t.gt.t0-17.5) a1koenig=0.9175e-5
	if(t.le.t0-17.5.and.t.gt.t0-18.5) a1koenig=0.4412e-5
	if(t.le.t0-18.5.and.t.gt.t0-19.5) a1koenig=0.2252e-5
	if(t.le.t0-19.5.and.t.gt.t0-20.5) a1koenig=0.9115e-6
	if(t.le.t0-20.5.and.t.gt.t0-21.5) a1koenig=0.4876e-6
	if(t.le.t0-21.5.and.t.gt.t0-22.5) a1koenig=0.3473e-6
	if(t.le.t0-22.5.and.t.gt.t0-23.5) a1koenig=0.4758e-6
	if(t.le.t0-23.5.and.t.gt.t0-24.5) a1koenig=0.6306e-6
	if(t.le.t0-24.5.and.t.gt.t0-25.5) a1koenig=0.8573e-6
	if(t.le.t0-25.5.and.t.gt.t0-26.5) a1koenig=0.7868e-6
	if(t.le.t0-26.5.and.t.gt.t0-27.5) a1koenig=0.7192e-6
	if(t.le.t0-27.5.and.t.gt.t0-28.5) a1koenig=0.6513e-6
	if(t.le.t0-28.5.and.t.gt.t0-29.5) a1koenig=0.5956e-6
	if(t.le.t0-29.5.and.t.gt.t0-30.5) a1koenig=0.5333e-6
	if(t.le.t0-30.5) a1koenig=0.4834e-6	
	end

	function a2koenig(t)
	implicit double precision (a-h,o-z)
	t0=273.15
	if(t.gt.t0) a2koenig=0.0
	if(t.le.t0.and.t.gt.t0-1.5) a2koenig=0.4006
	if(t.le.t0-1.5.and.t.gt.t0-2.5) a2koenig=0.4831
	if(t.le.t0-2.5.and.t.gt.t0-3.5) a2koenig=0.5320
	if(t.le.t0-3.5.and.t.gt.t0-4.5) a2koenig=0.5307
	if(t.le.t0-4.5.and.t.gt.t0-5.5) a2koenig=0.5319
	if(t.le.t0-5.5.and.t.gt.t0-6.5) a2koenig=0.5249
	if(t.le.t0-6.5.and.t.gt.t0-7.5) a2koenig=0.4888
	if(t.le.t0-7.5.and.t.gt.t0-8.5) a2koenig=0.3894
	if(t.le.t0-8.5.and.t.gt.t0-9.5) a2koenig=0.4047
	if(t.le.t0-9.5.and.t.gt.t0-10.5) a2koenig=0.4318
	if(t.le.t0-10.5.and.t.gt.t0-11.5) a2koenig=0.4771
	if(t.le.t0-11.5.and.t.gt.t0-12.5) a2koenig=0.5183
	if(t.le.t0-12.5.and.t.gt.t0-13.5) a2koenig=0.5463
	if(t.le.t0-13.5.and.t.gt.t0-14.5) a2koenig=0.5651
	if(t.le.t0-14.5.and.t.gt.t0-15.5) a2koenig=0.5813
	if(t.le.t0-15.5.and.t.gt.t0-16.5) a2koenig=0.5655
	if(t.le.t0-16.5.and.t.gt.t0-17.5) a2koenig=0.5478
	if(t.le.t0-17.5.and.t.gt.t0-18.5) a2koenig=0.5203
	if(t.le.t0-18.5.and.t.gt.t0-19.5) a2koenig=0.4906
	if(t.le.t0-19.5.and.t.gt.t0-20.5) a2koenig=0.4447
	if(t.le.t0-20.5.and.t.gt.t0-21.5) a2koenig=0.4126
	if(t.le.t0-21.5.and.t.gt.t0-22.5) a2koenig=0.3960
	if(t.le.t0-22.5.and.t.gt.t0-23.5) a2koenig=0.4149
	if(t.le.t0-23.5.and.t.gt.t0-24.5) a2koenig=0.4320
	if(t.le.t0-24.5.and.t.gt.t0-25.5) a2koenig=0.4506
	if(t.le.t0-25.5.and.t.gt.t0-26.5) a2koenig=0.4483
	if(t.le.t0-26.5.and.t.gt.t0-27.5) a2koenig=0.4460
	if(t.le.t0-27.5.and.t.gt.t0-28.5) a2koenig=0.4433
	if(t.le.t0-28.5.and.t.gt.t0-29.5) a2koenig=0.4413
	if(t.le.t0-29.5.and.t.gt.t0-30.5) a2koenig=0.4382
	if(t.le.t0-30.5) a2koenig=0.4361
	end
	
	function gamf(arg)
c----------------------------------------------------------c
c	! Gamma function calculation                                        
c----------------------------------------------------------c
	implicit none
	real*8 xn,slagaemoe,summa,gamf12,xmnojitel,gamf,arg2
	integer i
	real(8), intent(in) :: arg
	real*8 ag(8)
	 data ag/-0.577191652,0.988205891,-0.897056937,0.918206857,
     :    -0.756704078,0.482199394,-0.193527818,0.035868343/
     
       summa=0.
	 slagaemoe=0.
	 xmnojitel=0.
	 gamf12=0.
	 xn=0.
       arg2=arg-1
	 do while(arg2.gt.2)
	 arg2=arg2-1.
	 xn=xn+1
	 enddo
	 do i=1,8
	 
	 slagaemoe=ag(i)*(arg2-1)**i
	 summa=summa+slagaemoe
	 enddo
	 gamf12=1.+summa
	 xmnojitel=arg2+xn
	 do i=1,xn
       xmnojitel=xmnojitel*(arg2+xn-i)
	 enddo
	 gamf=xmnojitel*gamf12
	 summa=0.
	 slagaemoe=0.
	 xmnojitel=0.
	 gamf12=0.
	 xn=0.
	end
	
	function lamsnow(rho,qsnow)
	implicit none
	REAL(8), INTENT(IN) :: rho,qsnow
	real(8) pi,rosn,lamsnow,xnso
c------------------------------------------------------------------------------------------------------------------c
c        slope parameter for snow distribution, after Lin,1983
c------------------------------------------------------------------------------------------------------------------c
	
	pi=3.1421
	rosn=100. !1000.
	xnso=3.*10**6
	lamsnow=(pi*rosn*xnso/(rho*qsnow))**0.25
	end

	function lamrain(rho,qrain)
c------------------------------------------------------------------------------------------------------------------c
c        slope parameter for rain distribution, after Lin,1983
c------------------------------------------------------------------------------------------------------------------c
	implicit none
	real(8) pi,xnro,rhol,lamrain,rho,qrain
	pi=3.1421
	xnro=8.*10**6
	rhol=1.e3
	lamrain=(pi*rhol*xnro/(rho*qrain))**0.25
	end
	
	function Usn(rho,qsnow,gamma)
c------------------------------------------------------------------------------------------------------------------c
c     mass-weight mean terminate velocity of precipitating snow, after Lin, 1983
c------------------------------------------------------------------------------------------------------------------c
      implicit none
	real(8), intent(in):: qsnow,rho,gamma
	real(8) csnow,dsnow,rho0,Usn
	real(8),external:: lamsnow
	csnow=4.83
	dsnow=0.25
      rho0=1.23
	Usn=csnow*gamma*sqrt(rho0/rho)
     :/(6.*lamsnow(rho,qsnow)**dsnow)
	end
	
	function Urn(rho,qrain,gamma)
c------------------------------------------------------------------------------------------------------------------c
c     mass-weight mean terminate velocity of precipitating rain, after Lin, 1983
c------------------------------------------------------------------------------------------------------------------c
      implicit none
	real(8) arain,brain,rho0,Urn,qrain,rho,gamma
	real(8),external:: gamf,lamrain
	arain=842.
	brain=0.8
      rho0=1.23
	Urn=arain*gamma*sqrt(rho0/rho)
     :/(6.*lamrain(rho,qrain)**brain)
	end

