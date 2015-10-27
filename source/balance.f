      subroutine balance
      use ice_mod,only:
     : Ti,Tsi,Ki,dzi,
     : eps,sig,dzs,Tsn,Ks,
     : snow_h
      use alloc_1d, only:
     : th,ro,cp,ust_s,tst_s,
     : u,v
      implicit none
      real A,B,C
      integer i
      real f,fp
      real Ch,uvs
      real dz, T1
      real K_th

      uvs = sqrt(u(1,2)**2.+v(1,2)**2.)
!--------ice or snow on top-------------!
      if (snow_h.ne.0.) then
         T1 = Tsn(1)
         K_th = Ks
         dz = dzs
      else
         T1 = Ti(1)
         K_th = Ki
         dz = dzi
      endif
!-------first guess---------------------!      
      Tsi = T1
!--------begin iteration----------------!
      do i = 1,5
         call surf_layer_t
         Ch = - tst_s*ust_s/(uvs*(Tsi-th(1,2)))
         A = eps*sig
         B = ro(1)*cp*Ch*uvs + K_th/(0.5*dz)
         C = -ro(1)*cp*Ch*th(1,2)*uvs - K_th/(0.5*dz)*T1
         f = A*Tsi**4. + B*Tsi + C
         fp = 4.*A*Tsi**3. + B
         Tsi = Tsi - f/fp    
      enddo      

      end
