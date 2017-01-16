      subroutine richardson
      use alloc_1d
      implicit none
      real*8 def,xn,zero,rich
      real,external:: qsat
      integer iz

      ! def, shear - wind shear
      ! rich,ri - richardson number
      zero = 1.e-9
      
      do iz=2,nz-2
         def=dsqrt(((u(iz+1,2)-u(iz-1,2))
     :       /(z(iz+1)-z(iz-1)))**2.+
     :       ((v(iz+1,2)-v(iz-1,2))/(z(iz+1)-z(iz-1)))**2.)
 !       def=dsqrt(((u(iz,2)-u(iz-1,2))/(z(iz)-z(iz-1)))**2.+
 !    :         ((v(iz,2)-v(iz-1,2))/(z(iz)-z(iz-1)))**2.)
         if(qif.ne.0) then 
             if(ifwr.eq.0.) then              ! virtual theta
                xn=g*(th(iz+1,2)*(1.+0.61*qv(iz,2))
     :             -th(iz-1,2)*(1.+0.61*qv(iz,2)))
     :             /(z(iz+1)-z(iz-1))/
     :             th(iz,2)/(1.+0.61*qv(iz,2))
             else   
!                if(qv(iz,2)
!     :             .ge.0.9*qsat(t(iz),p(iz,2))) then
               if(qc(iz,2).gt.0.or.qc(iz+1,2).gt.0.or.qc(iz-1,2).gt.0)
     :          then
            xn=g*(1.+hlat*qsat(t(iz),p(iz,2))/r/t(iz))
     :    *(1.+0.622*hlat**2.*qsat(t(iz),p(iz,2))/cp/r/t(iz)**2.)**(-1.)
     :      *((log(th(iz+1,2))
     :      -log(th(iz-1,2)))/(z(iz+1)-z(iz-1))+
     :       hlat/cp/t(iz)*(qsat(t(iz+1),p(iz+1,2))-
     :       qsat(t(iz-1),p(iz-1,2))))/(z(iz+1)-z(iz-1))-
     :       g*(qv(iz+1,2)+qc(iz+1,2)-qv(iz-1,2)-qc(iz-1,2))/
     :       (z(iz+1)-z(iz-1))
!                    xn=g*((th(iz+1,2)-hlat/cp*qc(iz+1,2))
!     :                 -(th(iz-1,2)-hlat/cp*qc(iz-1,2)))
!     :                 /(z(iz+1)-z(iz-1))/th(iz,2)
 !              xn=g*(th(iz+1,2)
 !    :      -th(iz-1,2))
 !    :      /(z(iz+1)-z(iz-1))/
 !    :     th(iz,2)
                else
                   xn=g*(th(iz+1,2)*(1.+0.61*qv(iz+1,2))
     :                -th(iz-1,2)*(1.+0.61*qv(iz-1,2)))
     :                /(z(iz+1)-z(iz-1))/
     :                th(iz,2)/(1.+0.61*qv(iz,2))
                endif
            endif
         else
            xn=g*(th(iz+1,2)-th(iz-1,2))/(z(iz+1)-z(iz-1))/th(iz,2)
         endif
         rich=xn/(def*def+zero)
         ri(iz) = rich    ! richardson number
         shear(iz) = def  ! wind shear
      enddo
         
      end subroutine richardson
