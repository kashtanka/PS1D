      subroutine allocvar
      use alloc_1d
      implicit none
      if(allocated(u)) deallocate(u)
          allocate(u(0:nz,1:3))
          u=0.
      if(allocated(v)) deallocate(v)
          allocate(v(0:nz,1:3))
          v=0.
      if(allocated(th)) deallocate(th)
          allocate(th(0:nz,1:3))
          th=0.
      if(allocated(qv)) deallocate(qv)
          allocate(qv(0:nz,1:3))
          qv=0.
      if(allocated(qc)) deallocate(qc)
          allocate(qc(0:nz,1:3))
          qc=0.
      if(allocated(qr)) deallocate(qr)
          allocate(qr(0:nz,1:3))
          qr=0.
      if(allocated(qci)) deallocate(qci)
          allocate(qci(0:nz,1:3))
          qci=0.
      if(allocated(qsn)) deallocate(qsn)
          allocate(qsn(0:nz,1:3))
          qsn=0.
      if(allocated(qs)) deallocate(qs)
          allocate(qs(0:nz))
          qs=0.
      if(allocated(w)) deallocate(w)
          allocate(w(0:nz))
          w=0.
      if(allocated(rad)) deallocate(rad)
          allocate(rad(0:nz))
          rad=0.
      if(allocated(rfl)) deallocate(rfl)
          allocate(rfl(0:nz))
          rfl=0.
      if(allocated(flx_lw2)) deallocate(flx_lw2)
          allocate(flx_lw2(0:nz))
          flx_lw2=0.
      if(allocated(vat)) deallocate(vat)
          allocate(vat(0:nz))
          vat=0.
      if(allocated(vaq)) deallocate(vaq)
          allocate(vaq(0:nz))
          vaq=0.
      if(allocated(vaqc)) deallocate(vaqc)
          allocate(vaqc(0:nz))
          vaqc=0.
      if(allocated(vat2)) deallocate(vat2)
          allocate(vat2(0:nz))
          vat2=0.
      if(allocated(vaq2)) deallocate(vaq2)
          allocate(vaq2(0:nz))
          vaq2=0.
      if(allocated(vaqc2)) deallocate(vaqc2)
          allocate(vaqc2(0:nz))
          vaqc2=0.
      if(allocated(tfc)) deallocate(tfc)
          allocate(tfc(0:nz))
          tfc=0.
      if(allocated(vau)) deallocate(vau)
          allocate(vau(0:nz))
          vau=0.
      if(allocated(vav)) deallocate(vav)
          allocate(vav(0:nz))
          vav=0.
      if(allocated(p)) deallocate(p)
          allocate(p(0:nz,1:2))
          p=0.
      if(allocated(t)) deallocate(t)
          allocate(t(0:nz))
          t=0.
      if(allocated(ro)) deallocate(ro)
          allocate(ro(0:nz))
          ro=0.
      if(allocated(dpdy)) deallocate(dpdy)
          allocate(dpdy(0:nz))
          dpdy=0.
      if(allocated(vgeos)) deallocate(vgeos)
          allocate(vgeos(0:nz))
          vgeos=0.
      if(allocated(difunu)) deallocate(difunu)
          allocate(difunu(0:nz))
          difunu=0.
      if(allocated(condensat)) deallocate(condensat)
          allocate(condensat(0:nz))
          condensat=0.
      if(allocated(sublim)) deallocate(sublim)
          allocate(sublim(0:nz))
          sublim=0.
      if(allocated(difunv)) deallocate(difunv)
          allocate(difunv(0:nz))
          difunv=0.
      if(allocated(difunt)) deallocate(difunt)
          allocate(difunt(0:nz))
          difunt=0.
      if(allocated(difunqv)) deallocate(difunqv)
          allocate(difunqv(0:nz))
          difunqv=0.
      if(allocated(difunqr)) deallocate(difunqr)
          allocate(difunqr(0:nz))
          difunqr=0.
      if(allocated(difunqc)) deallocate(difunqc)
          allocate(difunqc(0:nz))
          difunqc=0.
      if(allocated(difunqci)) deallocate(difunqci)
          allocate(difunqci(0:nz))
          difunqci=0.
      if(allocated(difunqsn)) deallocate(difunqsn)
          allocate(difunqsn(0:nz))
          difunqsn=0.
!--------MICROPHYSICS VARIABLES-----------------!          
      if(allocated(cond)) deallocate(cond)
          allocate(cond(0:nz))
          cond=0.
      if(allocated(evap)) deallocate(evap)
          allocate(evap(0:nz))
          evap=0.
      if(allocated(auto)) deallocate(auto)
          allocate(auto(0:nz))
          auto=0.
      if(allocated(col)) deallocate(col)
          allocate(col(0:nz))
          col=0.
      if(allocated(sberci)) deallocate(sberci)
          allocate(sberci(0:nz))
          sberci=0.
      if(allocated(ibercw)) deallocate(ibercw)
          allocate(ibercw(0:nz))
          ibercw=0.
      if(allocated(sbercw)) deallocate(sbercw)
          allocate(sbercw(0:nz))
          sbercw=0.
      if(allocated(iacr)) deallocate(iacr)
          allocate(iacr(0:nz))
          iacr=0.
      if(allocated(raci)) deallocate(raci)
          allocate(raci(0:nz))
          raci=0.
      if(allocated(saci)) deallocate(saci)
          allocate(saci(0:nz))
          saci=0.
      if(allocated(sagg)) deallocate(sagg)
          allocate(sagg(0:nz))
          sagg=0.
      if(allocated(sacrw)) deallocate(sacrw)
          allocate(sacrw(0:nz))
          sacrw=0.
      if(allocated(sacrwr)) deallocate(sacrwr)
          allocate(sacrwr(0:nz))
          sacrwr=0.
      if(allocated(hmfrz)) deallocate(hmfrz)
          allocate(hmfrz(0:nz))
          hmfrz=0.
      if(allocated(smlt)) deallocate(smlt)
          allocate(smlt(0:nz))
          smlt=0.
      if(allocated(imlt)) deallocate(imlt)
          allocate(imlt(0:nz))
          imlt=0.
      if(allocated(sacrr)) deallocate(sacrr)
          allocate(sacrr(0:nz))
          sacrr=0.
      if(allocated(vdepi)) deallocate(vdepi)
          allocate(vdepi(0:nz))
          vdepi=0.
      if(allocated(vini)) deallocate(vini)
          allocate(vini(0:nz))
          vini=0.
      if(allocated(vdeps)) deallocate(vdeps)
          allocate(vdeps(0:nz))
          vdeps=0.
      if(allocated(divrain)) deallocate(divrain)
          allocate(divrain(0:nz))
          divrain=0.
      if(allocated(divsnow)) deallocate(divsnow)
          allocate(divsnow(0:nz))
          divsnow=0.
      if(allocated(vrain)) deallocate(vrain)
          allocate(vrain(0:nz))
          vrain=0.
      if(allocated(vsnow)) deallocate(vsnow)
          allocate(vsnow(0:nz))
          vsnow=0.
!------------------------------------------------!
      
      if(allocated(def13)) deallocate(def13)
          allocate(def13(0:nz))
          def13=0.
      if(allocated(def23)) deallocate(def23)
          allocate(def23(0:nz))
          def23=0.
      if(allocated(def13c)) deallocate(def13c)
          allocate(def13c(0:nz))
          def13c=0.
      if(allocated(def23c)) deallocate(def23c)
          allocate(def23c(0:nz))
          def23c=0.
      if(allocated(dz)) deallocate(dz)
          allocate(dz(0:nz))
          dz=0.
      if(allocated(ri)) deallocate(ri)
          allocate(ri(0:nz))
          ri=0.
      if(allocated(shear)) deallocate(shear)
          allocate(shear(0:nz))
          shear=0.
      if(allocated(difk)) deallocate(difk)
          allocate(difk(0:nz))
          difk=0.
      if(allocated(dift)) deallocate(dift)
          allocate(dift(0:nz))
          dift=0.
      if(allocated(dif_qc)) deallocate(dif_qc)
          allocate(dif_qc(0:nz))
          dif_qc=0.
      if(allocated(difk3)) deallocate(difk3)
          allocate(difk3(0:nz))
          difk3=0.
      if(allocated(dift3)) deallocate(dift3)
          allocate(dift3(0:nz))
          dift3=0.
      if(allocated(difk2)) deallocate(difk2)
          allocate(difk2(0:nz))
          difk2=0.
      if(allocated(dift2)) deallocate(dift2)
          allocate(dift2(0:nz))
          dift2=0.
      if(allocated(th0)) deallocate(th0)
          allocate(th0(0:nz))
          th0=0.
       if(allocated(u0)) deallocate(u0)
          allocate(u0(0:nz))
          u0=0.
      if(allocated(v0)) deallocate(v0)
          allocate(v0(0:nz))
          v0=0.
      if(allocated(qv0)) deallocate(qv0)
          allocate(qv0(0:nz))
          qv0=0.
      if(allocated(qc0)) deallocate(qc0)
          allocate(qc0(0:nz))
          qc0=0.
      if(allocated(qci0)) deallocate(qci0)
          allocate(qci0(0:nz))
          qci0=0.
      if(allocated(qsn0)) deallocate(qsn0)
          allocate(qsn0(0:nz))
          qsn0=0.
      if(allocated(qr0)) deallocate(qr0)
          allocate(qr0(0:nz))
          qr0=0.
      if(allocated(h3)) deallocate(h3)
          allocate(h3(0:nz))
          h3=0.
      if(allocated(fthl)) deallocate(fthl)
          allocate(fthl(0:nz))
          fthl=0.
      if(allocated(fqt)) deallocate(fqt)
          allocate(fqt(0:nz))
          fqt=0.
      if(allocated(HF)) deallocate(HF)
          allocate(HF(0:nz))
          HF=0.
       if(allocated(HF2)) deallocate(HF2)
          allocate(HF2(0:nz))
          HF2=0.
      if(allocated(h3c)) deallocate(h3c)
          allocate(h3c(0:nz))
          h3c=0.
      if(allocated(h3e)) deallocate(h3e)
          allocate(h3e(0:nz))
          h3e=0.
      if(allocated(wq3)) deallocate(wq3)
          allocate(wq3(0:nz))
          wq3=0.
      if(allocated(wq3e)) deallocate(wq3e)
          allocate(wq3e(0:nz))
          wq3e=0.
      if(allocated(wq3c)) deallocate(wq3c)
          allocate(wq3c(0:nz))
          wq3c=0.
      if(allocated(wq3_c)) deallocate(wq3_c)
          allocate(wq3_c(0:nz))
          wq3_c=0.
      if(allocated(wq3ci)) deallocate(wq3ci)
          allocate(wq3ci(0:nz))
          wq3ci=0.
      if(allocated(wq3r)) deallocate(wq3r)
          allocate(wq3r(0:nz))
          wq3r=0.
      if(allocated(wq3sn)) deallocate(wq3sn)
          allocate(wq3sn(0:nz))
          wq3sn=0.
      if(allocated(sh3)) deallocate(sh3)
          allocate(sh3(0:nz))
          sh3=0.
      if(allocated(ht)) deallocate(ht)
          allocate(ht(0:nz))
          ht=0.
      if(allocated(mom)) deallocate(mom)
          allocate(mom(0:nz))
          mom=0.
      if(allocated(ug_bar)) deallocate(ug_bar)
          allocate(ug_bar(0:nz))
          ug_bar=0.   
      cond2=0.
      if(allocated(z)) goto 77
          allocate(z(0:nz))
          z=0.
77    continue
      h=0.
      le=0.
      end
