      subroutine progonka(A,B,C,F,kappa1,nu1,kappa2,nu2,difunt)
      use alloc_1d, only:
     :  nz,dtl
      implicit none
      integer i
!     input arrays and params
      real*8 A(1:nz),B(1:nz),C(1:nz),F(1:nz),difunt(1:nz)
      real*8 kappa1,nu1,kappa2,nu2
!     output
      real*8 y(1:nz)
      real*8 alpha(1:nz), beta(1:nz)
      y=0.      

! PRYAMOY HOD
      alpha(2)=kappa1
      beta(2)=nu1
      do i = 2,nz-1
         alpha(i+1)=B(i)/(C(i)-alpha(i)*A(i))
         beta(i+1)=(A(i)*beta(i)+F(i))/(C(i)-alpha(i)*A(i))
      enddo
! OBRATNIY HOD
      y(nz)=(nu2+kappa2*beta(nz))/(1.-kappa2*alpha(nz))
      do i = nz-1,1,-1
         y(i)=alpha(i+1)*y(i+1)+beta(i+1)
      enddo
      do i = 1,nz
         difunt(i)=(y(i)-F(i))/dtl
         write(0,*) y(i),F(i)
      enddo
      end
