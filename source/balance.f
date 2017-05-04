      subroutine balance
      use ice_mod,only:
     : Ti,Tsi,Ki,dzi,
     : eps,sig,dzs,Tsn,Ks,
     : snow_h
      use alloc_1d, only:
     : th,ro,cp,ust_s,tst_s,
     : u,v,ust_s2,tst_s2,nstep,dt,LW,
     :  rad_par
      implicit none
      real A,B,C
      integer i
      real f,fp
      real Ch,uvs,Ch2
      real dz, T1
      real K_th
      real nn, aKL, bKL ! nn-total cloud amount 0-1; aKL,bKL Konig-Langlo koefficients 
      real hr
      hr = nstep*dt/3600
!      if (hr.le.48) then
!         nn=0.
!      else
!         nn = 0.
!      endif
      nn = 0.9
      aKL=0.765
      bKL=0.22 

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
      if(nstep.eq.1) Tsi = T1
!--------begin iteration----------------! 
!      do i = 1,5
!         call surf_layer_t
!         Ch = - tst_s*ust_s/(uvs*(Tsi-th(1,2)))
!         A = eps*sig
!         B = ro(1)*cp*Ch*uvs + K_th/(0.5*dz)
!         C = -ro(1)*cp*Ch*th(1,2)*uvs - K_th/(0.5*dz)*T1 - (aKL+bKL*nn
!     :       **3)*eps*sig*(th(1,2))**4
!         f = A*Tsi**4. + B*Tsi + C
!         fp = 4.*A*Tsi**3. + B
!         Tsi = Tsi - f/fp    
!      enddo      

cccc-  VERSION2 WHERE LW BALANCE IS PRESCRIBED -cccc
      if(rad_par.ne.2) then
      LW = -eps*sig*Tsi**4. + (aKL+bKL*nn**3)
     :       *eps*sig*(th(1,2))**4
      endif
      do i = 1,5
         call surf_layer_t
         Ch = - tst_s*ust_s/(uvs*(Tsi-th(1,2)))
         Tsi = (LW + K_th/(0.5*dz)*T1 + 
     :         ro(1)*cp*Ch*th(1,2)*uvs)/
     :          (ro(1)*cp*Ch*uvs + K_th/(0.5*dz))
      enddo
   ! write(0,*) Tsi
      end
