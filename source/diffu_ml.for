      subroutine diffu_ml
      use alloc_1d
      implicit none
      real,external:: qsat
      real*8 wstar,thv,w_m3,wg1,wg2,zf,t_hl
      real*8 pot_hl,alpha,Pr,tstar_f,ws,Pr0,eps2
      real*8 ws0,wsm,eps
      integer iz
      real*8 dift_hl(1:nz),difk_hl(1:nz)
      real, parameter:: B2 = 5.
      real, parameter:: karm =0.4
      real, parameter:: al = 3.
      eps=z_sl/hbl
      if(qif.gt.0) then
        thv=th(1,2)+0.61*th(1,2)*qv(1,2)
      else
        thv=th(1,2)
      endif
      !--------scaling parameters----------------------------------
      wstar=(g/thv*Fv*hbl)**(1./3.) 
      ws0 = (ust_s**3. + 7.*eps*karm*wstar**3.)**(1./3.)
      wsm=(ust_s**3. + 7.*karm*wstar**3.*0.5)**(1./3.)
      tstar_f=-ust_s*tst_s/wstar
!---------ENTRAINMENT----------------------------------!      
      w_m3=(wstar**3.+B2*ust_s**3.) +hbl*dR*g/thv
      wth_h= -6.*w_m3/hbl*2.5
       do iz = 2, nz-1
         if (0.5*(z(iz)+z(iz+1))
     :      .ge.hbl.and.0.5*(z(iz-1)+z(iz)).lt.hbl) then
            w_e = wth_h/(delta_th+0.61*th(iz,2)*delta_qv)
            wthl_h = -w_e*(delta_th - hlatcp*delta_qc)
            wqt_h = -w_e*(delta_qv) ! + delta_qc)
        endif
      enddo
!-------------------------------------------------------!
!--------------FLUXES OF CONSERVATIVE VARS--------------!
      fthl(1) = ust_s*tst_s
      fqt(1) = ust_s*qst_s
      do iz = 2,nz
         zf = 0.5*(z(iz)+z(iz-1))
         wg1 = max(0.,(hbl - zf)/(hbl - 0.5*z(iz)))
         wg2 = min(1.,(zf - 0.5*z(iz))/(hbl - 0.5*z(iz)))
         if (0.5*(z(iz)+z(iz-1)).le.hbl) then
            fthl(iz) = fthl(1)*wg1 + wthl_h*wg2
            fqt(iz) = fqt(1)*wg1 + wqt_h*wg2
         else
            fthl(iz) = 0.
            fqt(iz) = 0.
         endif
         if (qc(iz,1).gt.0.and.qc(iz-1,1).gt.0) then
            t_hl = 0.5*(t(iz)+t(iz-1))
            pot_hl = t_hl/(0.5*(th(iz,2)+th(iz-1,2)))
            alpha = hlat*0.5*(qsat(t(iz),p(iz,2))
     :                 +qsat(t(iz-1),p(iz-1,2)))/rv/t_hl/t_hl
            h3(iz) = (fthl(iz)
     :                 +hlatcp*fqt(iz)/pot_hl)    !1./pot_hl*
     :                    /(1.+hlatcp*alpha)
            wq3(iz) = alpha*h3(iz)*pot_hl
            wq3c(iz) = fqt(iz) - wq3(iz)
         else
            h3(iz) = fthl(iz)
            wq3(iz) = fqt(iz)
            wq3c(iz) = 0.
         endif
      enddo

      h3(1)=ust_s*tst_s ! surface flux of heat
         wq3(1)=ust_s*qst_s !surface flux of moisture
      
      do iz=1,nz-1
         difunt(iz)=(h3(iz+1)
     :       -h3(iz))/(0.5*(dz(iz+1)+dz(iz)))
         difunqv(iz)=(wq3(iz+1)
     :       -wq3(iz))/(0.5*(dz(iz+1)+dz(iz)))
      enddo
       
      if (ifwr.ne.0) then
         do iz=1,nz-1
            difunqc(iz)=(wq3c(iz+1)-wq3c(iz))/(0.5*(dz(iz+1)+dz(iz)))
         enddo
      endif

      do iz=2,nz-2
      !-------------K-profile----------------------------------------!
          if(Fv.gt.0.and.0.5*(z(iz+1)+z(iz+2)).le.hbl) then
	    ws=(ust_s**3.+7.*karm*wstar**3.*z(iz)/hbl)**(1./3.) 
	    Pr=1.+(Pr0-1.)*exp(-al*(z(iz)-eps*hbl)**2./hbl**2.)
            difk(iz)=karm*ws*z(iz)*(1.-eps2*z(iz)/hbl)**2.
            dift(iz)= difk(iz)/Pr
          !  difk(iz) = max(difk(iz),difk2(iz))
          !  dift(iz) = max(dift(iz),dift2(iz))
          !  write(0,*) iz,difk(iz)
         else
           dift(iz)=0.
           difk(iz)=0.
         endif
          if(0.5*(z(iz+2)+z(iz+1)).le.hbl) then
           dift3(iz) = 0.85*karm*(dR*hbl*g/thv)**(1./3.)*
     :               (z(iz))**2.
     :              /(hbl)*(1.-(z(iz))
     :              /(hbl))**0.5
           difk3(iz) = 0.75*dift3(iz)
         else 
           dift3(iz)=0.
           difk3(iz)=0.
         endif
         difk(iz)=difk(iz)+difk3(iz)
         dift(iz)=dift(iz)+dift3(iz)
         !write(0,*) dift3(iz),dift(iz)
      enddo

      do iz = 2,nz
         dift_hl(iz)=0.5*(dift(iz-1)+dift(iz))
         difk_hl(iz)=0.5*(difk(iz-1)+difk(iz))
      enddo

      do iz=2,nz
         def13(iz)=(u(iz,2)-u(iz-1,2))/dz(iz)
         def23(iz)=(v(iz,2)-v(iz-1,2))/dz(iz)
         def13(iz)=def13(iz)*difk_hl(iz)   ! momentum flux
         def23(iz)=def23(iz)*difk_hl(iz)
      enddo
      def13(1)=cdm*u(1,2) !surface flux of momentum
      def23(1)=cdm*v(1,2)      
      do iz=1,nz-1
         difunu(iz)=(def13(iz+1)-def13(iz))/(0.5*(dz(iz+1)+dz(iz)))
         difunv(iz)=(def23(iz+1)-def23(iz))/(0.5*(dz(iz+1)+dz(iz)))
      enddo

      end
