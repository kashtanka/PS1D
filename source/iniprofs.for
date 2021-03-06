      subroutine iniprofs
      use alloc_1d
      implicit none
      integer iz,id,it
      real dtd,dzd,dvd,dud,dqd,height,dtfc
      real*8 F, Fwt,pm,C1,C2,DIV
      real dzz,dtt,u1,u2
      real,external:: qsati
      do id=0,ndat-1
      dtd=thdat(id+1)-thdat(id)
      dzd=zthdat(id+1)-zthdat(id)
      dud=usdat(id+1)-usdat(id)
      dvd=vsdat(id+1)-vsdat(id)
      dqd=qvsdat(id+1)-qvsdat(id)
      do iz=0,nz
      height=(z(iz)-z_sl)
      if (height.ge.zthdat(id).and.
     : height.le.zthdat(id+1)) then
      th0(iz)=thdat(id)+dtd*(height-zthdat(id))/dzd
      u0(iz)=usdat(id)+dud*(height-zthdat(id))/dzd
      v0(iz)=vsdat(id)+dvd*(height-zthdat(id))/dzd
      qv0(iz)=qvsdat(id)+dqd*(height-zthdat(id))/dzd
      endif
      
      
      enddo
      enddo

      if (iftf) then
         do id = 0,ntfc
            dtfc = tfcdat(id+1) - tfcdat(id)
            dzd=zfc(id+1)-zfc(id)
            do iz=0,nz
               height=(z(iz)-z_sl)
               if (height.ge.zfc(id).and.
     :             height.le.zfc(id+1)) then
                  tfc(iz)=tfcdat(id)+dtfc*(height - zfc(id))/dzd
                  tfc(iz)=tfc(iz)*dt/3600.
                  write(0,*) height,tfc(iz)
               endif
            enddo
         enddo
      endif
      th0(0)=th0(1)
      
      do iz=0,nz
      write(0,*)'t0',z(iz)-z_sl,th0(iz)
      th(iz,1)=th0(iz)
      th(iz,2)=th0(iz)
      th(iz,3)=th0(iz)
      u(iz,1)=u0(iz)
      u(iz,2)=u0(iz)
      u(iz,3)=u0(iz)
      v(iz,1)=v0(iz)
      v(iz,2)=v0(iz)
      v(iz,3)=v0(iz)
      if(qif.gt.0) then
      qv(iz,1)=qv0(iz)
      qv(iz,2)=qv0(iz)
      qv(iz,3)=qv0(iz)
      else
      qv(iz,1)=0.
      qv(iz,2)=0.
      qv(iz,3)=0.
      endif
      enddo
      if(inifile.eq.0) then
      !!-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-!!
      !! initialization for marine Sc case   !!
      !!-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-!!
      do iz = 0,nz
         if (z(iz)-z_sl.le.840) then
            qv(iz,2)= 9.e-3
            th(iz,2)=289.
         else
            qv(iz,2)= 1.5e-3
            th(iz,2)=297.5+(z(iz)-z_sl-840.)**(1./3.)
         endif
         qv(:,1)=qv(:,2)
         qv(:,3)=qv(:,2)
         th(:,1)=th(:,2)
         th(:,3)=th(:,2)
         u(:,:)=ug
         v(:,:)=vg
      enddo
      endif
!     vertical wind speed at half-levels
      DIV=3.75e-6      ! large-scale divergence
      w(1)=0.
      do iz=2,nz
         w(iz)=-DIV*0.5*(z(iz-1)+z(iz)-z_sl)
      enddo
      
      !---calculate initial pressure profile to define ptop given pa----!
      !    pressure is calculated using hydrostatic equation and
      !    definition of potential temperature. 
      
      p(1,2)=pa   !known value
      write(0,*) p(1,2)
      do iz = 2,nz
         dzz = z(iz)-z(iz-1)
         dtt = th(iz,2) - th(iz-1,2)
         u2 = th(iz,2)
         u1 = th(iz-1,2)
         if( dtt.ne.0) then
         p(iz,2) = p(iz-1,2)**akapa - g*akapa*p00**akapa/r*dzz/dtt*
     :        (log(u2)-log(u1))
         else
          p(iz,2) = p(iz-1,2)**akapa - g*akapa*p00**akapa/r*dzz
     :        /th(iz-1,2)  
         endif
         p(iz,2) = p(iz,2)**(1./akapa)
!         write(0,*) 'dzz=','dtt=',dzz,dtt,th(iz,2)

!      do iz=1,nz-1
!      C1=g*(z(iz+1)-z(iz))
!     :    0.5*p00**akapa/r/(0.5*(th(iz+1,2)+th(iz,2))) !constant during iteration at this z-level
!      C2=p(iz,2)                                                  !constant during iteration at this z-level
      
!      pm=p(iz,2)-1000.                                          !first guess
!      F=pm+C1*pm**(1.-akapa)-C2                                 ! function value with first guess
!      Fwt=1.+C1*(1.-akapa)*pm**(-akapa)                         ! derivative value with first guess
      !-----iteration------------!
!      do it=1,10
!      pm=pm-F/Fwt
!      F=pm+C1*pm**(1.-akapa)-C2
!      Fwt=1.+C1*(1.-akapa)*pm**(-akapa)
!      end do
!      p(iz+1,2)=2.*pm-p(iz,2)
      enddo
      p(:,1)=p(:,2)
      ptop=p(nz,2)
      
      do iz =1,nz
      t(iz)=th(iz,2)*(p(iz,2)/p00)**akapa
      write(0,*) z(iz),p(iz,2),t(iz)
      !if (qif.ne.0) then
!      qv(iz,1)=0.8*qsati(t(iz),p(iz,1))
!      qv(iz,2)=0.8*qsati(t(iz),p(iz,1))
!      qv(iz,3)=0.8*qsati(t(iz),p(iz,1))
      !endif
      !vgeos(iz)=z(iz)/z(nz)*vg
      !v(iz,1)=vgeos(iz)
      !v(iz,2)=vgeos(iz)
      !v(iz,3)=vgeos(iz)
      enddo
      end
