      subroutine radiation
      use alloc_1d
      implicit none
      integer iz
      real*8,parameter :: xkp = 85.
      real*8,parameter :: DIV = 3.75e-6
      real*8,parameter :: F0 = 70.
      real*8,parameter :: F1 = 22.
      real*8,parameter:: gq = 1.7e-6
      real*8 zi,Q1,Q0,Q2,R3,mdR,dR2,LWP,LWP2,qmax,dRmax
      real*8 Fplus
!-----------------------------------------!
!-!    net longwave radiative flux for  !-!
!-!    Sc experiment at half-levels     !-!
!-----------------------------------------!
      Q0=0.
      Q2=0.
      LWP=0.
      qmax=0.
      do iz=1,nz-1
      Q0=Q0+xkp*ro(iz)*qc(iz,2)*0.5*(dz(iz)+dz(iz+1))
      enddo

      do iz=2,nz-1
!--------determine qmax for LWP calc--------------!
        if (z(iz).le.hbl) then
           LWP = LWP+ro(iz)*qc(iz,2)*1000.*dz(iz)
           if (qc(iz,2).gt.qmax) qmax = qc(iz,2)
        endif       
         if(qv(iz,2)+qc(iz,2).ge.8.e-3.and.
     :      qv(iz+1,2)+qc(iz+1,2).lt.8.e-3) then
            zi=0.5*(z(iz)+z(iz+1))
            !zi = hbl
            endif
            zi=max(100.,zi)
      enddo
      LWP2=ro(1)*qmax**2./(2.*gq)*1000.
!      zi=800.
      Q1=Q0
      rfl(0)=F0*exp(-Q1)+F1
      do iz=2,nz-1
         Q1=Q1-xkp*ro(iz)*qc(iz-1,2)*0.5*(dz(iz)+dz(iz+1))
         Q2=Q2+xkp*ro(iz)*qc(iz-1,2)*0.5*(dz(iz)+dz(iz+1))
         if(0.5*(z(iz)+z(iz-1)).ge.zi) then
         R3=0.5*(ro(iz)+ro(iz-1))
     :           *cp*DIV*((0.5*(z(iz-2)+z(iz-1))-zi)**(4./3.)/4.
     :           +zi*(0.5*(z(iz-1)+z(iz-2))-zi)**(1./3.))
         else
         R3=0.
         endif
         
         rfl(iz)=F0*exp(-Q1)+F1*exp(-Q2) +R3
      enddo
      dR=0.
      mdR=0.
      do iz=1,nz-1
        rad(iz)=-(rfl(iz+1)-rfl(iz))/(0.5*(dz(iz+1)+dz(iz)))/cp/ro(iz)
!         dR=(rfl(iz)-rfl(iz-1))/cp/ro(iz)
!         if (dR.gt.mdR) mdR=dR
         if (0.5*(z(iz)+z(iz-1)).lt.hbl.and.
     :            0.5*(z(iz)+z(iz+1)).gt.hbl) then
            dR=(rfl(iz)-rfl(iz-1))/cp/ro(iz)
            Fplus = rfl(iz)
         endif
      enddo
      dRmax = 0.
      do iz=1,nz-1
         if(qc(iz,2).gt.0) then
            if ((Fplus - rfl(iz))/cp/ro(iz).gt.dRmax) 
     :          dRmax = (Fplus - rfl(iz))/cp/ro(iz)
         endif
      enddo
      dRmax = 55/cp/ro(1)
       dR2 = Fplus*(1. - exp(-0.03*sqrt(LWP)))
       dR=dRmax
!      rad=0.
      
      end
