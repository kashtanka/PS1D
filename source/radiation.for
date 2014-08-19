      subroutine radiation
      use alloc_1d
      implicit none
      integer iz
      real*8,parameter :: xkp = 85.
      real*8,parameter :: DIV = 3.75e-6
      real*8,parameter :: F0 = 70.
      real*8,parameter :: F1 = 22.
      real*8 zi,Q1,Q0,Q2,R3,mdR,dR2
!-----------------------------------------!
!-!    net longwave radiative flux for  !-!
!-!    Sc experiment at half-levels     !-!
!-----------------------------------------!
      Q0=0.
      Q2=0.
      do iz=1,nz-1
      Q0=Q0+xkp*ro(iz)*qc(iz,2)*0.5*(dz(iz)+dz(iz+1))
      enddo

      do iz=1,nz-1
         if(qv(iz,2)+qc(iz,2).ge.8.e-3.and.
     :      qv(iz+1,2)+qc(iz+1,2).lt.8.e-3) then
            zi=z(iz)
            endif
            zi=max(100.,zi)
      enddo
!      zi=800.
      Q1=Q0
      rfl(0)=F0*exp(-Q1)+F1
      do iz=1,nz
         Q1=Q1-xkp*ro(iz)*qc(iz,2)*0.5*(dz(iz)+dz(iz+1))
         Q2=Q2+xkp*ro(iz)*qc(iz,2)*0.5*(dz(iz)+dz(iz+1))
         if(z(iz).ge.zi) then
         R3=0.5*(ro(iz)+ro(iz-1))
     :           *cp*DIV*((z(iz)-zi)**(4./3.)/4.+zi*(z(iz)-zi)**(1./3.))
         else
         R3=0.
         endif
         
         rfl(iz)=F0*exp(-Q1)+F1*exp(-Q2)+R3
!         write(0,*) z(iz),rfl(iz),Q1
      enddo
      dR=0.
      mdR=0.
      do iz=1,nz-1
         rad(iz)=-(rfl(iz)-rfl(iz-1))/(0.5*(dz(iz+1)+dz(iz)))/cp/ro(iz)
!         dR=(rfl(iz)-rfl(iz-1))/cp/ro(iz)
!         if (dR.gt.mdR) mdR=dR
!         write(0,*) z(iz),rad(iz),qc(iz,2)
         if (z(iz).lt.hbl.and.z(iz+1).ge.hbl) then
            dR=(rfl(iz)-rfl(iz-3))/cp/ro(iz)
!            write(0,*) dR,(rfl(iz)-rfl(iz-3))/cp/ro(iz)
         endif
      enddo
!      write(0,*) dR2,mdR
!      rad=0.
      
      end
