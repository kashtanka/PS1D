      subroutine vert_adv
      use alloc_1d
      implicit none
      integer iz
      real gambl, gamfa,thplus,thminus,zist

      do iz=2,nz-1                        ! vertical advection
         vat(iz)= w(iz)*(th(iz+1,2)-th(iz,2))/
     :            (0.5*(dz(iz+1)+dz(iz)))
!0.5*(w(iz+1)+w(iz))*(0.5*(th(iz+1,1)+th(iz,1))
!     :           -0.5*(th(iz,1)+th(iz-1,1)))/
!     :           (0.5*(dz(iz+1)+dz(iz)))
!w(iz)*(th(iz+1,1)-
!     :   th(iz,1))/
!     :   (0.5*(dz(iz)+dz(iz+1)))
         vau(iz)=w(iz)*(u(iz+1,2)-
     :   u(iz,2))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vav(iz)=w(iz)*(v(iz+1,2)-
     :   v(iz,2))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vaq(iz)=w(iz)*(qv(iz+1,2)-
     :   qv(iz,2))/
     :   (0.5*(dz(iz)+dz(iz+1)))
         vaqc(iz)=w(iz)*(qc(iz+1,2)-
     :   qc(iz,2))/
     :   (0.5*(dz(iz)+dz(iz+1)))
      enddo
      vat(nz)=0.

      vat2=vat
      vaq2=vaq
      vaqc2=vaqc
      do iz = 2, nz-1
 !        if(z(iz).gt.hbl.and.z(iz-1).lt.hbl) then
 !         if (0.5*(z(iz)+z(iz+1)).gt.zi_rec
 !    :       .and.0.5*(z(iz)+z(iz-1)).le.zi_rec) then
         if (0.5*(z(iz)+z(iz+1))
     :      .gt.hbl.and.0.5*(z(iz-1)+z(iz)).le.hbl) then
       
            gambl=(th(iz-1,1)-th(iz-2,1))/dz(iz-1)
            gamfa=(th(iz+2,1)-th(iz+1,1))/dz(iz+2)
            vat2(iz)= w(iz)*delta_th/dz(iz)   !(th(iz+1,1) -0.5*dz(iz+1)*gamfa
 !    :              -(th(iz-1,1) +0.5*gambl*dz(iz)
 !    :                ))/(2*dz(iz))*w(iz)
            vat2(iz-1) = w(iz-1)*gambl
!            if(zi_rec+(0.5*(w(iz)+w(iz+1)))*dtl/2
!     :                           .lt.0.5*(z(iz)+z(iz-1))) then
!               zist = zi_rec+0.5*(w(iz)+w(iz+1))*dtl/2
!               thplus = th(iz+1,2)-gamfa*(z(iz+1)-zi_rec)
!               thminus = th(iz-1,2)+gambl*(zi_rec - z(iz-1))
!               vat(iz)= -(thplus + gamfa*(z(iz)-zist) - th(iz,1))/dtl
!               vat(iz-1)=-(((thminus + gambl*0.5*(z(iz-1)+z(iz-2)))*
!     :                   (zist -0.5*(z(iz-1)+z(iz-2))) +
!     :                     (thplus + gamfa*0.5*(z(iz)+z(iz-1)))*
!     :                     (0.5*(z(iz)+z(iz-1))-zist))
!     :                  /(0.5*(z(iz)-z(iz-2))) - th(iz-1,1))/dtl
!               else
!
!            vat(iz) = 0.5*(w(iz)+w(iz+1))
!     :                 *(delta_th + gambl*(zi_rec-0.5*(z(iz)+z(iz+1)))
!     :                   +gamfa*(0.5*(z(iz+1)+z(iz+2))-zi_rec))/dz(iz)
!     :                -w_e*(delta_th)/dz(iz)
!!            write(0,*) 'wzi=',0.5*(w(iz)+w(iz+1))
!            write(0,*) 'we=', wth_h/(delta_th+0.61*th(iz,2)*delta_qv)
! 0.5*(w(iz)+w(iz+1))*(th(iz+1,1)-(thzi+delta_th))/
 !    :            (z(iz+1) - zi_rec)
!            vat(iz-1) =0.5*(w(iz-1)+w(iz))*gambl+              !*(thzi-th(iz-1,1))/
!     :            (zi_rec - z(iz-1))+
!     :             w_e*(delta_th)/dz(iz-1)
!            endif
            gambl = (qv(iz-1,1)-qv(iz-2,1))/dz(iz-1)
            gamfa = (qv(iz+2,1)-qv(iz+1,1))/dz(iz+2)
!            vaq(iz) =   0.5*(w(iz)+w(iz+1))
!     :                   *(delta_qv+ gambl*(zi_rec-0.5*(z(iz)+z(iz+1)))
!     :                   +gamfa*(0.5*(z(iz+1)+z(iz+2))-zi_rec))/dz(iz)
!     :                  -w_e
!     :                  *(delta_qv)/dz(iz)
            xint1=vat(iz)+vat(iz-1)
            xint2=vat2(iz)+vat2(iz-1)
!0.5*(w(iz)+w(iz+1))*(qv(iz+1,1)-(qvzi+delta_qv))/
!     :            (z(iz+1) - zi_rec)
!            vaq(iz-1) =0.5*(w(iz-1)+w(iz))*gambl+    !*(qvzi-qv(iz-1,1))/
!     :            (zi_rec - z(iz-1))
!     :                  +w_e
!     :                   *(delta_qv)/dz(iz-1)
            vaq2(iz)= (qv(iz+1,1) !-0.5*dz(iz+1)*gamfa
     :              -(qv(iz-1,1) !+0.5*gambl*dz(iz)
     :               ))/(2*dz(iz))*w(iz)
            vaq2(iz-1) = w(iz-1)*gambl
            gambl = (qc(iz-1,2)-qc(iz-2,2))/dz(iz-1)
            gamfa = (qc(iz+2,2)-qc(iz+1,2))/dz(iz+2)
            vaqc2(iz) = 0.5*(w(iz)+w(iz+1))
     :                  *(qc(iz+1,1)-(qczi+delta_qc)+
     :                gambl*(zi_rec-0.5*(z(iz)+z(iz+1)))
     :                   +gamfa*(0.5*(z(iz+1)+z(iz+2))-zi_rec))/
     :            (z(iz+1) - zi_rec)
            vaqc2(iz-1) =0.5*(w(iz-1)+w(iz))*(qczi-qc(iz-1,1))/
     :            (zi_rec - z(iz-1))
            
         endif
      enddo
      vat=vat2
!      vaq=0.
!      vaqc=0.
      end

