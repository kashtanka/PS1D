      subroutine surf_layer_t
      use alloc_1d
      use ice_mod, only:
     : Tsi
      implicit none
      real*8 uvs,ta,qa,ps,z0,cdu,hflux,Elatent,
     : Ribulk,dzita,rhoa,z1,tau,water,tems,ustar,Tstar,
     : qstar,temsfc
      real,external:: qsat
      integer i
      
      water=0.
      ta=th(1,2)
      if (seaice.eq.1) then
         ts = Tsi
	else
           do i = 1,36
              if (nstep*dt.ge.gabls_tim(i).and.
     :           nstep*dt.lt.gabls_tim(i+1)) then
                 ts = gabls_ts(i) + (nstep*dt-gabls_tim(i))/3600.*
     :                (gabls_ts(i+1) - gabls_ts(i))
              endif
           enddo
!          if(dt*nstep.lt.3600*icetime) then
!             ts = 250.
!             water = 0.
!          else
!             ts = 273.15
!             water = 1.
!          endif
           
      endif
      write(0,*) 'ts = ',ts
 !     ts=235.
      z0=0.01
 !     ts = !241.15-dt*nstep/3600.*0.25
      
!      do i = 1,100
!      if(dt*nstep.gt.3600*12.and.dt*nstep.gt.3600.*icetime*i.and.
!     :  dt*nstep.le.3600*(icetime+0.01)*i)
!     : then
!         ts=268.15
!         water=1.
!      goto 111
!      endif
!      enddo
! 111  continue
      th(0,3)=ts*(p00/p(1,2))**0.286
      th(0,2)=ts*(p00/p(1,2))**0.286
      th(0,1)=ts*(p00/p(1,2))**0.286
           
      uvs=sqrt(u(1,2)**2.+v(1,2)**2.)
      ps=p(1,2)
      if (qif.ne.0) then
        qa=qv(1,2)
      else
        qa=0.
      endif
      rhoa=ro(1)
       
      call surf_scheme3(water,
     :       ta,ts,qa,
     :       uvs,ps,z_sl,z0,
     :       cdu,hflux,Elatent,ustar,qif,Ribulk,dzita,Tstar,rhoa,z1
     :       ,tau,qstar)
     
      cdm=cdu
      ust_s=ustar
      tst_s=Tstar
      if (qif.ne.0) then
         qst_s=qstar
      else
         qst_s=0.
      endif
      dzits=dzita
      h=hflux
      le=Elatent
      Fv=-ust_s*tst_s-0.61*t(1)*ust_s*qst_s
      
      if (seaice.eq.1.and.frac.lt.1) then
         water = 1
         ts = 271.35
         call surf_scheme3(water,
     :       ta,ts,qa,
     :       uvs,ps,z_sl,z0,
     :       cdu,hflux,Elatent,ustar,qif,Ribulk,dzita,Tstar,rhoa,z1
     :       ,tau,qstar)

         cdm2 = cdu
         ust_s2 = ustar
         tst_s2 = Tstar
         if (qif.ne.0) then
            qst_s2 = qstar
         else
            qst_s2 = 0.
         endif
         dzits2 = dzita
         h2 = hflux
         le2 = Elatent
         Fv2 = -ust_s2*tst_s2-0.61*t(1)*ust_s2*qst_s2
       else
          cdm2 = 0.
          ust_s2 = 0.
          tst_s2 = 0.
          qst_s2 = 0.
          dzits2 = 0.
          h2 = 0.
          le2 = 0.
          Fv2 = 0.
         
       endif
         
      end
      
      subroutine surf_scheme3(xsea,temp2,temp1,q2,uv2,p,z2,z0,cdu,hflux,
     :	                   Elatent,ustar,qif,Ribulk,dzita,Tstar,ro,z1,
     :                     tau,qstar)
	implicit none
	real(8) t2,t1,q2,u2,p,z2,z1,xsea
	real(8) k,g,pi,cp,nu,Le,Lsub,z0min,aMagw,bMagw,
     :aMagi,bMagi,am,ah,bm,bh,B_m,B_h,ch,q1,aMag,bMag
	real(8) dq,dT,Ribulk,p00,l1,ro,Tsr, dU,esat1,u1,uv2,temp1,temp2,z0
	real(8) L,dzita1,dzita2,xst1,xst2,unifU1,unifU2,unifT1,unifT2,
     :	unifq1,unifq2,ustar,Tstar,qstar,R,b0,b1,b2,xx,zt,dzita1T,Lit,
     :    dL,x1,x2,kanzasU1,kanzasU2,yU1,yU2,convU1,convU2,hflux,Elatent
     :    ,tau,yT1,yT2,convT1,convT2,kanzasT1,kanzasT2,Uneust,w,cdu,qif
     :    ,dzita
	integer paramzt,param_uf,conv, paramz1,strat,iter,metras
        !***************************************************************************!
	  ! options for z1 - roughness length for momentum 
	  ! paramz1
	  ! 1 - z1 = z0   - read from file
	  ! 2 - METRAS Charnock
	  ! 3 - Charnock's formulae for sea surface
	  !***************************************************************************!
	  ! paramzt - parameterization of roughness parameter for heat flux
	  ! 0 - z0t = z0m
	  ! 1 - don't remember who
	  ! 2 - Zilitinkevich
	  ! 3 - Andreas - for sea ice and snow
	  ! 4 - Beljaars
	  ! 5 - metras
	  !***************************************************************************!
	  !param_uf - universal functions for stable stratification
	  ! 1 - Grachev (SHEBA) (2007)
	  ! 2 - Beljaars and Holtslag (1991)
	  ! 3 - Cheng and Brutsaert (2005)
	  ! 4 - Log-linear
	  ! 5 - Zilitinkevich and Esau (2006)
	  ! 6 - Holtslag and de Bruin (1988)
	  ! 7 - log-linear as in GABLS
	  !****************************************************************************!
	  ! metras - universal functions for unstable stratification from METRAS (Businger)
	  ! 1 - on
	  ! 0 - off
	  !****************************************************************************!
	  ! conv - adding to wind speed in unstable stratification due to convection
	  ! 1 - on
	  ! 0 - off
	  !*****************************************************************************! 
!	  write(0,*)xsea,temp1
	  if (xsea.eq.1) then
	  paramzt=5
	  paramz1=2
	  else
	  paramzt=0
	  paramz1=1
	  endif
	  conv=0
	  param_uf=2.
	  metras=1
	!-----------------------------------------------------------------------------!
	  p00=1.e5
	  k=0.4
	  g=9.81
	  pi=3.14159
	  cp=1005.
	  nu=1.5/10**5
	  Le=2.501e6
	  Lsub=2.837e6
	  z0min=1.5e-5
	  aMagw = 7.6326   ! coefficient for Magnus formula for water surface
      bMagw = 241.9      ! coefficient for Magnus formula for water surface
      aMagi = 9.5        ! coefficient for Magnus formula for ice   surface
      bMagi = 265.5      ! coefficient for Magnus formula for ice   surface
	  !----------------SHEBA const---------------!
	  am=5.
	  ah=5.
	  bm=am/6.5
	  bh=5.
	  B_m=((1-bm)/bm)**(1./3.)
	  B_h=sqrt(5.)
	  ch=3.
	  !------------------------------------------!
	
	  iter=0.
	  u1=0.
	  t1=temp1
	  t2=temp2
	  u2=uv2
	  if(u2.eq.0) u2=0.1
	  z1=max(z0,0.00001)
!	  if(t1.ge.273.15) then
        if(t1.ge.273.15.or.xsea.gt.0) then
	  l1=0.
	  else
	  l1=1.
	  endif	  
	  if(l1.eq.0.) then
	  aMag=aMagw
	  bMag=bMagw
	  else
        aMag=aMagi
	  bMag=bMagi
	  endif
	  
	  t1=t1*(p00/p)**0.286      ! t1 - fixed real surface temperature (T) being
	                            ! converted into potential temperature (theta)
!	  t2=t2**(p00/p)**0.286
        esat1= 610.7*10.**(aMag*(t1-273.15)/(bMag+(t1-273.15)))
        q1  = 0.622/p*esat1

	  Tsr=(t1+t2)/2.
	  dT=t2-t1
	  if(qif.gt.0) then
	  dq=q2-q1
	else
	  dq=0.
	q1=0.
	q2=0.
	  endif
	  dU=u2-u1
	  !-------Bulk Richardson Number------------!
        Ribulk=z2*g*((t2*(1.+0.61*q2))-(t1*(1.+0.61*q1)))
     :  /(0.5*(t1*(1.+0.61*q1)+t2*(1.+0.61*q2)))/u2**2
	  !
	  !-------DEFINE STRATIFICATION-------------!
	  ! 1 - stable
	  ! 2 - unstable
	  ! neutral is not necessary
	  !-----------------------------------------!
	  if(dT+0.61*dq*Tsr.gt.0.) strat=1           
	  if(dT+0.61*dq*Tsr.le.0.) strat=2     
        !-----------------------------------------!
	if (strat.eq.1) then   
	  L=10000.    
20	  iter=iter+1
        dzita2=z2/L
	  dzita1=z1/L

	                     ! STABLE STRATIFICATION
      if(param_uf==1) then
	   xst1=(1.+dzita1)**(1./3.)
	   xst2=(1.+dzita2)**(1./3.)
	   unifU1=-3.*am/bm*(xst1-1)+am*B_m/(2.*bm)*(2*dlog((xst1+B_m)
     :   /(1+B_m))-dlog((xst1**2-xst1*B_m+B_m**2)/(1-B_m+B_m**2))
     :   +2*sqrt(3.)*(atan((2.*xst1-B_m)/(sqrt(3.)*B_m))
     :   -atan((2.-B_m)/(sqrt(3.)*B_m))))
	   unifU2=-3.*am/bm*(xst2-1)+am*B_m/(2.*bm)*(2*dlog((xst2+B_m)
     :   /(1+B_m))-dlog((xst2**2-xst1*B_m+B_m**2)/(1-B_m+B_m**2))
     :   +2*sqrt(3.)*(atan((2.*xst2-B_m)/(sqrt(3.)*B_m))
     :   -atan((2.-B_m)/(sqrt(3.)*B_m))))
	 else if(param_uf==2) then
	   unifU1=-(1.*dzita1+0.667*(dzita1-5./0.35)
     :   *dexp(-0.35*dzita1)+0.667*5./0.35)
	   unifU2=-(1.*dzita2+0.667*(dzita2-5./0.35)
     :   *dexp(-0.35*dzita2)+0.667*5./0.35)
	 else if(param_uf==3) then
	   unifU1=-6.1*dlog(dzita1+(1+dzita1**2.5)**(1./2.5))
	   unifU2=-6.1*dlog(dzita2+(1+dzita2**2.5)**(1./2.5))
	 else if(param_uf==4) then
	   unifU1=-5.*dzita1
	   unifU2=-5.*dzita2
!           write(0,*) unifU2
      
	 else if(param_uf==5) then
	   unifU1=-(3.*dzita1**(5./6.))
	   unifU2=-(3.*dzita2**(5./6.))
	 else if(param_uf==6) then
         unifU1=-(0.7*dzita1+0.75*(dzita1-5./0.35)
     :   *dexp(-0.35*dzita1)+0.75*5./0.35)
	   unifU2=-(0.7*dzita2+0.75*(dzita2-5./0.35)
     :   *dexp(-0.35*dzita2)+0.75*5./0.35)
       else if(param_uf==7) then
	   unifU1=-4.8*dzita1
	   unifU2=-4.8*dzita2
	   endif
       ustar=k*dU/(dlog(z2/z1)-unifU2+unifU1)
	 if(paramz1.eq.1) then
	   z1=z1
	 elseif (paramz1.eq.2) then
	  z1=0.0185*ustar**2./g
	 elseif (paramz1.eq.3) then
	   z1 = dmin1(dmax1(0.111*nu/ustar + 
     :   0.0144*ustar**2./g, 1.d-5),1.1d-1)
	 endif
	  R=ustar*z1/nu
	  if (R.lt.0.135) then
	  b0=1.25
	  b1=0.
	  b2=0.
	  endif
	  if(R.gt.0.135.and.R.lt.2.5) then
	  b0=0.149
	  b1=-0.55
	  b2=0.
	  endif
	  if(R.gt.2.5) then
	  b0=0.317
	  b1=-0.565
	  b2=-0.183
	  endif
	  if (paramzt.eq.0) zt=z1
	  if (paramzt.eq.1) zt=z1/exp(6.5)
	  if (paramzt.eq.2) zt=z1/exp(0.13*R**0.45)
	  if (paramzt.eq.3) zt=z1*exp(b0+b1*dlog(R)+b2*(dlog(R))**2.)
	  if (paramzt.eq.4) then
	  if(R .le. 0.111) then
			xx = - 2.43
		else
			if(0.111 .lt. R .and. R .le. 16.3) then
				xx = 0.83*log(R) - 0.6
			else
				xx = 0.49 * R**0.45
			end if
		end if
		zt = z1*exp(-xx)
		endif
	  if (paramzt.eq.5) zt=0.1*z1
	  dzita1T=zt/L
	  if (param_uf==1) then
	  unifT1=-bh/2.*dlog(1+ch*dzita1T+dzita1T**2)+
     :  (-ah/B_h+bh*ch/(2.*B_h))*(dlog((2.*dzita1T+ch-B_h)
     :  /(2*dzita1T+ch+B_h))-dlog((ch-B_h)/(ch+B_h)))
        unifT2=-bh/2.*dlog(1+ch*dzita2+dzita2**2)+
     :  (-ah/B_h+bh*ch/(2.*B_h))*(dlog((2.*dzita2+ch-B_h)
     :  /(2*dzita2+ch+B_h))-dlog((ch-B_h)/(ch+B_h)))
	  else if(param_uf==2) then
	  unifT1=-((1.+2./3.*1.*dzita1T)**(3./2.)+
     :  0.667*(dzita1T-5./0.35)*exp(-0.35*dzita1T)+0.667*5./0.35-1.)
        unifT2=-((1+2./3.*1.*dzita2)**(3./2.)+
     :  0.667*(dzita2-5./0.35)*exp(-0.35*dzita2)+0.667*5./0.35-1.)
	  else if (param_uf==3) then
	  unifT1=-5.3*dlog(dzita1T+(1.+dzita1T**1.1)**(1./1.1))
        unifT2=-5.3*dlog(dzita2+(1.+dzita2**1.1)**(1./1.1))
	  else if (param_uf==4) then
	  unifT1=-5.*dzita1T
        unifT2=-5.*dzita2
	  else if (param_uf==5) then
	  unifT1=-(2.5*dzita1T**(4./5.))
        unifT2=-(2.5*dzita2**(4./5.))
        else if (param_uf==6) then
	  unifT1=-(0.7*dzita1T+0.75*(dzita1T-5./0.35)
     :  *dexp(-0.35*dzita1T)+0.75*5./0.35)
        unifT2=-(0.7*dzita2+0.75*(dzita2-5./0.35)
     :  *dexp(-0.35*dzita2)+0.75*5./0.35)
        else if (param_uf==7) then
	  unifT1=-7.8*dzita1T
        unifT2=-7.8*dzita2
	  endif
	  unifq1=unifT1
	  unifq2=unifT2
	  Tstar=(k)*dT/(dlog(z2/zt)-unifT2+unifT1)
	  qstar=k*dq/(dlog(z2/zt)-unifq2+unifq1)
	  Lit=ustar**2/(k*g*(Tstar+0.61*Tsr*qstar))*Tsr
	  Lit = sign(1.,Lit)*max(abs(Lit),1.e-7) 
	  dL=dabs(L-Lit)

	  if (iter.gt.25) goto 33
	  
	  if (dL.lt.5) then
	     if(L.gt.z2) then
                if(dL.lt.1) goto 33
                L= (Lit+L)/2.
             else
                if(dL.lt.0.00001) goto 33
                L=Lit
             endif
          endif
	  if(iter.gt.1) then
	     L= (Lit+L)/2.
	  else
	     L=Lit
	  endif
          goto 20
	  endif

      
	  if (strat.eq.2) then                   !UNSTABLE STRATIFICATION
	  L=-10000.
	  !-----------------------------------!
	  if(conv.eq.1) then
	  ustar=k*dU/dlog(z2/z1)
	  tstar=k*dT/dlog(z2/zt)
	  w=(-g/Tsr*tstar*ustar*1000.)**(1/3)
	  Uneust=sqrt(u2**2+(w*1.)**2)
	  dU=Uneust-u1
        endif

	  !-----------------------------------!
	  
 40     iter=iter+1
        dzita2=z2/L
	  dzita1=z1/L 	  
	  x2=(1.-16.*dzita2)**(1./4.)
        x1=(1.-16.*dzita1)**(1./4.)
	  kanzasU2=2.*dlog(0.5*(1.+x2))+dlog(0.5*(1+x2**2))-2.*atan(x2)+pi/2.
	  kanzasU1=2.*dlog(0.5*(1.+x1))+dlog(0.5*(1+x1**2))-2.*atan(x1)+pi/2.	  
        yU2=(1.-10.15*dzita2)**(1./3.)
	  yU1=(1.-10.15*dzita1)**(1./3.)	  
	  convU2=(3./2.)*dlog((1./3.)*(yU2**2+yU2+1.))-sqrt(3.)
     :  *dlog((2.*yU2+1.)/sqrt(3.))+pi/sqrt(3.)
        convU1=(3./2.)*dlog((1./3.)*(yU1**2+yU1+1.))-sqrt(3.)
     :  *dlog((2.*yU1+1.)/sqrt(3.))+pi/sqrt(3.)     
        if (metras.eq.0.) then
	  unifU2=(kanzasU2+dzita2**2*convU2)/(1+dzita2**2)
        unifU1=(kanzasU1+dzita1**2*convU1)/(1+dzita1**2)  
	else 
	  unifU2=kanzasU2
	  unifU1=kanzasU1
	  endif    
	  ustar=k*dU/(dlog(z2/z1)-unifU2+unifU1)
       if(paramz1.eq.1) then
	   z1=z1
	 elseif (paramz1.eq.2) then
	  z1=0.0185*ustar**2./g
	 elseif (paramz1.eq.3) then
	   z1 = dmin1(dmax1(0.111*nu/ustar + 
     :   0.0144*ustar**2./g, 1.d-5),1.1d-1)
	 endif
	  R=ustar*z1/nu
	  if (R.lt.0.135) then
	  b0=1.25
	  b1=0.
	  b2=0.
	  endif
	  if(R.gt.0.135.and.R.lt.2.5) then
	  b0=0.149
	  b1=-0.55
	  b2=0.
	  endif
	  if(R.gt.2.5) then
	  b0=0.317
	  b1=-0.565
	  b2=-0.183
	  endif
	  if (paramzt.eq.0) zt=z1
	  if (paramzt.eq.1) zt=z1/exp(6.5)
	  if (paramzt.eq.2) zt=z1/exp(0.13*R**0.45)
	  if (paramzt.eq.3) zt=z1*exp(b0+b1*dlog(R)+b2*(dlog(R))**2)
	  if (paramzt.eq.4) then
	  if(R .le. 0.111) then
			xx = - 2.43
		else
			if(0.111 .lt. R .and. R .le. 16.3) then
				xx = 0.83*log(R) - 0.6
			else
				xx = 0.49 * R**0.45
			end if
		end if
		zt = z1*exp(-xx)
		endif
	 if (paramzt.eq.5) zt=0.1*z1
      dzita1T=zt/L
	  kanzasT2=2.*dlog(0.5*(1.+sqrt(1.-16.*dzita2)))
      kanzasT1=2.*dlog(0.5*(1.+sqrt(1.-16.*dzita1T)))
	  yT2=(1.-34.15*dzita2)**(1./3.)
	  yT1=(1.-34.15*dzita1T)**(1./3.)
	  convT2=(3./2.)*dlog((1./3.)*(yT2**2+yT2+1.))
     :	  -sqrt(3.)*dlog((2.*yT2+1.)/sqrt(3.))+pi/sqrt(3.)
	  convT1=(3./2.)*dlog((1./3.)*(yT1**2+yT1+1.))
     :	  -sqrt(3.)*dlog((2.*yT1+1.)/sqrt(3.))+pi/sqrt(3.)
	  if (metras.eq.0) then
	  unifT2=(kanzasT2+dzita2**2.*convT2)/(1+dzita2**2.)
	  unifT1=(kanzasT1+dzita1T**2.*convT1)/(1+dzita1T**2.)
	  else
	  unifT2=kanzasT2
	  unifT1=kanzasT1
	  endif
	  unifq2=unifT2
	  unifq1=unifT1
	  Tstar=(k)*dT/(dlog(z2/zt)-unifT2+unifT1)
	  qstar=k*dq/(dlog(z2/zt)-unifq2+unifq1)
	  Lit=ustar**2./(k*g*(Tstar+0.61*Tsr*qstar))*Tsr
	  Lit = sign(1.,Lit)*max(abs(Lit),1.e-3) 
      dL=dabs(L-Lit)
	  
	   if (iter.gt.10) goto 35
	  
	  if (dL.lt.5) then
	     if(L.gt.z2) then
		 if(dL.lt.1) goto 35
		 L=(Lit+L)/2.
             else
		 if(dL.lt.0.2) goto 35
		 L=Lit
             endif
          endif
	  if(iter.gt.1) then
	     L=(Lit+L)/2.
	  else
	     L=Lit
	  endif
      goto 40
	  
	  endif
      
	  if (strat.eq.3) then
      ustar=k*dU/dlog(z2/z1)
	  tstar=k*dT/dlog(z2/zt)
	  goto 30
	  endif

33    L=Lit
       dzita2=z2/L
	   dzita1=z1/L 	  
       if(param_uf==1) then
	   xst1=(1.+dzita1)**(1./3.)
	   xst2=(1.+dzita2)**(1./3.)
	   unifU1=-3.*am/bm*(xst1-1.)+am*B_m/(2.*bm)*(2*dlog((xst1+B_m)
     :   /(1+B_m))-dlog((xst1**2-xst1*B_m+B_m**2.)/(1-B_m+B_m**2))+
     :   2.*sqrt(3.)*(atan((2.*xst1-B_m)/(sqrt(3.)*B_m))
     :   -atan((2.-B_m)/(sqrt(3.)*B_m))))
	   unifU2=-3.*am/bm*(xst2-1.)+am*B_m/(2.*bm)*(2.*dlog((xst2+B_m)
     :   /(1+B_m))-dlog((xst2**2.-xst1*B_m+B_m**2.)/(1-B_m+B_m**2))+
     :   2.*sqrt(3.)*(atan((2.*xst2-B_m)/(sqrt(3.)*B_m))
     :   -atan((2.-B_m)/(sqrt(3.)*B_m))))
	   else if(param_uf==2) then
	   unifU1=-(1.*dzita1+0.667*(dzita1-5./0.35)*
     :   dexp(-0.35*dzita1)+0.667*5./0.35)
	   unifU2=-(1.*dzita2+0.667*(dzita2-5./0.35)*
     :   dexp(-0.35*dzita2)+0.667*5./0.35)
	   else if(param_uf==3) then
	   unifU1=-6.1*dlog(dzita1+(1+dzita1**2.5)**(1./2.5))
	   unifU2=-6.1*dlog(dzita2+(1+dzita2**2.5)**(1./2.5))
	   else if(param_uf==4) then
	   unifU1=-5.*dzita1
	   unifU2=-5.*dzita2
	   else if(param_uf==5) then
	   unifU1=-(3.*dzita1**(5./6.))
	   unifU2=-(3.*dzita2**(5./6.))
	   else if(param_uf==6) then
         unifU1=-(0.7*dzita1+0.75*(dzita1-5./0.35)*
     :   dexp(-0.35*dzita1)+0.75*5./0.35)
	   unifU2=-(0.7*dzita2+0.75*(dzita2-5./0.35)*
     :   dexp(-0.35*dzita2)+0.75*5./0.35)
         else if(param_uf==7) then
	   unifU1=-4.8*dzita1
	   unifU2=-4.8*dzita2
	   endif
	  ustar=k*dU/(dlog(z2/z1)-unifU2+unifU1)
	  R=ustar*z1/nu
	  if (R.lt.0.135) then
	  b0=1.25
	  b1=0.
	  b2=0.
	  endif
	  if(R.gt.0.135.and.R.lt.2.5) then
	  b0=0.149
	  b1=-0.55
	  b2=0.
	  endif
	  if(R.gt.2.5) then
	  b0=0.317
	  b1=-0.565
	  b2=-0.183
	  endif
	  if (paramzt.eq.0) zt=z1
	  if (paramzt.eq.1) zt=z1/exp(6.5)
	  if (paramzt.eq.2) zt=z1/exp(0.13*R**0.45)
	  if (paramzt.eq.3) zt=z1*exp(b0+b1*dlog(R)+b2*(dlog(R))**2)
	  if (paramzt.eq.4) then
	  if(R .le. 0.111) then
			xx = - 2.43
		else
			if(0.111 .lt. R .and. R .le. 16.3) then
				xx = 0.83*log(R) - 0.6
			else
				xx = 0.49 * R**0.45
			end if
		end if
		zt = z1*exp(-xx)
		endif
	 if (paramzt.eq.5) zt=0.1*z1
		dzita1T=zt/L
	  if (param_uf==1) then
	    unifT1=-bh/2.*dlog(1+ch*dzita1T+dzita1T**2)+
     :    (-ah/B_h+bh*ch/(2.*B_h))*(dlog((2*dzita1T+ch-B_h)
     :    /(2*dzita1T+ch+B_h))-dlog((ch-B_h)/(ch+B_h)))
          unifT2=-bh/2.*dlog(1+ch*dzita2+dzita2**2)
     :	+(-ah/B_h+bh*ch/(2.*B_h))*(dlog((2*dzita2+ch-B_h)
     :    /(2*dzita2+ch+B_h))-dlog((ch-B_h)/(ch+B_h)))
	  else if(param_uf==2) then
	    unifT1=-((1+2./3.*1.*dzita1T)**(3./2.)+
     :    0.667*(dzita1T-5./0.35)*exp(-0.35*dzita1T)+0.667*5./0.35-1.)
          unifT2=-((1+2./3.*1.*dzita2)**(3./2.)+
     :    0.667*(dzita2-5./0.35)*exp(-0.35*dzita2)+0.667*5./0.35-1.)
	  else if (param_uf==3) then
	    unifT1=-5.3*dlog(dzita1T+(1+dzita1T**1.1)**(1./1.1))
          unifT2=-5.3*dlog(dzita2+(1+dzita2**1.1)**(1./1.1))
	  else if (param_uf==4) then
          unifT1=-5.*dzita1T
          unifT2=-5.*dzita2
	  else if (param_uf==5) then
	    unifT1=-(2.5*dzita1T**(4./5.))
          unifT2=-(2.5*dzita2**(4./5.))
	  else if (param_uf==6) then
	    unifT1=-(0.7*dzita1T+0.75*(dzita1T-5./0.35)
     :    *dexp(-0.35*dzita1T)+0.75*5./0.35)
          unifT2=-(0.7*dzita2+0.75*(dzita2-5./0.35)
     :    *dexp(-0.35*dzita2)+0.75*5./0.35)
        else if (param_uf==7) then
	    unifT1=-7.8*dzita1T
          unifT2=-7.8*dzita2
	  endif
	  unifq1=unifT1
	  unifq2=unifT2
	  Tstar=(k)*dT/(dlog(z2/zt)-unifT2+unifT1)
	  qstar=k*dq/(dlog(z2/zt)-unifq2+unifq1)
	  goto 30

35      L=Lit
        dzita2=z2/L
	  dzita1=z1/L 	  
	  x2=(1.-16.*dzita2)**(1./4.)
        x1=(1.-16.*dzita1)**(1./4.)
	  kanzasU2=2*dlog(0.5*(1+x2))+dlog(0.5*(1+x2**2))-2*atan(x2)+pi/2.
	  kanzasU1=2*dlog(0.5*(1+x1))+dlog(0.5*(1+x1**2))-2*atan(x1)+pi/2.	  
        yU2=(1.-10.15*dzita2)**(1./3.)
	  yU1=(1.-10.15*dzita1)**(1./3.)	  
	  convU2=(3./2.)*dlog((1./3.)*(yU2**2+yU2+1))-sqrt(3.)
     :  *dlog((2*yU2+1)/sqrt(3.))+pi/sqrt(3.)
        convU1=(3./2.)*dlog((1./3.)*(yU1**2+yU1+1))-sqrt(3.)
     :  *dlog((2*yU1+1)/sqrt(3.))+pi/sqrt(3.)     
        if (metras.eq.0.) then
	  unifU2=(kanzasU2+dzita2**2*convU2)/(1+dzita2**2)
        unifU1=(kanzasU1+dzita1**2*convU1)/(1+dzita1**2)  
	  else 
	  unifU2=kanzasU2
	  unifU1=kanzasU1
	  endif      
	  ustar=k*dU/(dlog(z2/z1)-unifU2+unifU1)
	  R=ustar*z1/nu
	  if (R.lt.0.135) then
	  b0=1.25
	  b1=0.
	  b2=0.
	  endif
	  if(R.gt.0.135.and.R.lt.2.5) then
	  b0=0.149
	  b1=-0.55
	  b2=0.
	  endif
	  if(R.gt.2.5) then
	  b0=0.317
	  b1=-0.565
	  b2=-0.183
	  endif
	  if (paramzt.eq.0) zt=z1
	  if (paramzt.eq.1) zt=z1/exp(6.5)
	  if (paramzt.eq.2) zt=z1/exp(0.13*R**0.45)
	  if (paramzt.eq.3) zt=z1*exp(b0+b1*dlog(R)+b2*(dlog(R))**2)
	  if (paramzt.eq.4) then
	  if(R .le. 0.111) then
			xx = - 2.43
		else
			if(0.111 .lt. R .and. R .le. 16.3) then
				xx = 0.83*log(R) - 0.6
			else
				xx = 0.49 * R**0.45
			end if
		end if
		zt = z1*exp(-xx)
		endif
	  if (paramzt.eq.5) zt=0.1*z1
        dzita1T=zt/L
	  kanzasT2=2.*dlog(0.5*(1.+sqrt(1.-16.*dzita2)))
        kanzasT1=2.*dlog(0.5*(1.+sqrt(1.-16.*dzita1T)))
	  yT2=(1.-34.15*dzita2)**(1./3.)
	  yT1=(1.-34.15*dzita1T)**(1./3.)
	  convT2=(3./2.)*dlog((1./3.)*(yT2**2+yT2+1))
     :  -sqrt(3.)*dlog((2*yT2+1.)/sqrt(3.))+pi/sqrt(3.)
	  convT1=(3./2.)*dlog((1./3.)*(yT1**2+yT1+1))
     :  -sqrt(3.)*dlog((2*yT1+1.)/sqrt(3.))+pi/sqrt(3.)
	  if (metras.eq.0) then
	  unifT2=(kanzasT2+dzita2**2*convT2)/(1+dzita2**2)
	  unifT1=(kanzasT1+dzita1T**2*convT1)/(1+dzita1T**2)
	  else
	  unifT2=kanzasT2
	  unifT1=kanzasT1
	  endif
	  unifq1=unifT1
	  unifq2=unifT2
	  Tstar=(k)*dT/(dlog(z2/zt)-unifT2+unifT1)
	  qstar=k*dq/(dlog(z2/zt)-unifq2+unifq1)
	  goto 30

30    continue
        hflux=-1.*ro*cp*ustar*Tstar
	  if (l1==0) then
	  Elatent=-1.*Le*ro*ustar*qstar
	  else
	  Elatent=-1.*Lsub*ro*ustar*qstar
	  endif
	  tau=ro*ustar**2
	  cdu=ustar**2/u2
	  dzita=dzita2

!	write(*,*) 'surf_layer',dzita,dT

	  
	  iter=0
	  
31    continue
      
	
	
	
	return
	end
