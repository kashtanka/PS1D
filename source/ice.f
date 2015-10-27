      subroutine ice_model
      use ice_mod
      use alloc_1d, only:
     : cp,ro,dt,ust_s,tst_s,th
      implicit none
      character*100 outfile
      real turb
      
      integer iz,it
      integer hour 
      

c --------- conductive flux of heat at half levels in ice --------c
         do iz = 2,nl
            fc(iz) = Ki*(Ti(iz)-Ti(iz-1))/dzi
         enddo
c--------- conductive flux of heat at half levels in snow --------c
         if (snow_h.gt.0) then
            do iz = 2,nls
               fcs(iz) = Ks*(Tsn(iz) - Tsn(iz-1))/dzs
            enddo
         endif
c--------- fluxes at the boundaries ------------------------------c
         turb = -ro(1)*cp*ust_s*tst_s
c--------- upper ice and upper/lower snow boundaries -------------c
         if (snow_h.gt.0) then
            T_int = (Ki*Ti(1)/dzi+Ks*Tsn(nls)/dzs)/(Ki/dzi+Ks/dzs)
            fc(1) = -Ki*(T_int - Ti(1))/(0.5*dzi)
            fcs(nls+1) = -Ks*(Tsn(nls) - T_int)/(0.5*dzs)
            fcs(1) = Ks*(Tsn(1)-Tsi)/(0.5*dzs)             
         else
            fc(1) = Ki*(Ti(1)-Tsi)/(0.5*dzi)
         endif
c------------lower ice boundary-----------------------------------c
            fc(nl+1) = Ki*(Tb - Ti(nl))/(0.5*dzi) ! think about it
c------------flux divergency and time integration-------------------c
         do iz = 1,nl
            dfc(iz) = (fc(iz+1)-fc(iz))/dzi
            Ti(iz) = Ti(iz) + dt*dfc(iz)/(roi*ci)
            !write(0,*) iz,Ti(iz)
         enddo
         if (snow_h.gt.0) then
            do iz = 1,nls
               dfcs(iz) = (fcs(iz+1)-fcs(iz))/dzs
               Tsn(iz) = Tsn(iz) + dt*dfcs(iz)/(rosn*ci)
               !write(0,*) iz, Tsn(iz)
            enddo
         endif
c----------- time integration ------------------------------------c
!	   if (nstep*dt.lt.36000) then
         
         
 !        else
!	   do iz = 1,nl
!            Ti(iz) = Ti(iz) + dt*dfc(iz)/(roi*ci)
!         enddo
!	   endif

      
      
      end
