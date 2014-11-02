      subroutine implicit_dif(F,F1,difk,sfx,difunt)
      use alloc_1d, only:
     : nz,dz,ust_s,tst_s,dtl
      implicit none
      real*8 F(1:nz),difk(1:nz),difunt(1:nz),F1(1:nz)
      real*8 D1,kappa1,nu1,Dn,kappa2,nu2,sfx,ed
      real*8 A(1:nz),B(1:nz),C(1:nz)
      integer iz
      ed=1.

! boundary conditions for forward phase
         
         D1=dtl*difk(2)/(0.5*(dz(1)+dz(2))*dz(2))
         kappa1=D1/(D1+ed)
         nu1=-D1*dz(2)/difk(2)/(D1+ed)*sfx+1./(D1+1.)*F(1)
! boundary conditions for backward phase
         Dn=dtl*difk(nz)/(dz(nz)**2.)
         kappa2=Dn/(Dn+ed)
         nu2=F(nz)/(Dn+ed)
! coefficients
         do iz = 2,nz-1
            A(iz)=dtl*difk(iz)/(0.5*(dz(iz)+dz(iz+1))*dz(iz))
            B(iz)=dtl*difk(iz+1)/(0.5*(dz(iz)+dz(iz+1))*dz(iz+1))
            C(iz)=A(iz)+B(iz)+ed
         enddo
         call progonka(A,B,C,F,F1,kappa1,nu1,kappa2,nu2,difunt)
      end
