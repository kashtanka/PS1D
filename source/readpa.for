      subroutine readpa(fnmap,fnfor,fnout
     :   ,itheta,iwind,scrout,fngrd,difloc,nonloc_tm86,
     : nonloc_ls96,nonloc_noh03,loc_inm,nonloc_lock)
      use alloc_1d
      implicit none
      logical scrout,dif3d,diftx
	logical difloc,nonloc_tm86,nonloc_ls96,nonloc_noh03,
     :          loc_inm,nonloc_lock  
      character*80 fnmap,fnfor,fext,fnout,fngrd
      character*80 line
      character*20 keyword,uppercase
      integer i,icha,ichz,itext,inperr,ichlef,ichrig,nchark,ii,
     : itheta,iwind,iz
      open(29,file='ps1d.dat',status='old')
      open(27,status='scratch')
      icha=ichar('0')
      ichz=ichar('z')

      itext=0
      inperr=0
      do while(inperr.eq.0)
        read(29,'(a)') line
c        write(*,*) line
c
c isolate keyword
c
        ichlef=1
        do i=1,80
          if(line(i:i).ne.' ') exit
          ichlef=i
        enddo
        if(ichlef.eq.80) cycle
        ichrig=ichlef
        do i=ichlef+1,80
          if(ichar(line(i:i)).lt.icha.or.ichar(line(i:i)).gt.ichz) exit
          ichrig=i
        enddo
        keyword=line(ichlef:ichrig)
        nchark=ichrig-ichlef+1
        !write(*,*) keyword,nchark,ichlef,ichrig
        write(27,'(a)') line(ichrig+1:80)
        !write(*,'(a)') line(ichrig+1:80)
        backspace(27)
        if(keyword(1:3).eq.'end') then
          exit
        elseif(keyword(1:1).eq.'#') then
          cycle
        elseif(keyword(1:nchark).eq.'ns') then
          read(27,*,iostat=inperr) ns
          write(*,*) 'ns=',ns
        elseif(keyword(1:nchark).eq.'pa') then
          read(27,*,iostat=inperr) pa
          write(*,*) 'pa=',pa
          p0=pa
        elseif(keyword(1:nchark).eq.'nz') then
          read(27,*,iostat=inperr) nz
          write(*,*) 'nz=',nz
        elseif(keyword(1:nchark).eq.'dt') then
          read(27,*,iostat=inperr) dt
          write(*,*) 'dt=',dt
        elseif(keyword(1:nchark).eq.'ntime') then
          read(27,*,iostat=inperr) ntime
          write(*,*) 'ntime=',ntime
        elseif(keyword(1:nchark).eq.'zgrid') then
          if(allocated(z)) deallocate(z)
          allocate(z(0:nz))
          z=0.
          read (29,*) (z(iz),iz=1,nz)
          z_sl=0.5*z(1)
          write(0,*) z          
        elseif(keyword(1:nchark).eq.'ztop') then
          read(27,*,iostat=inperr) ztop
          write(*,*) 'ztop',ztop
        elseif(keyword(1:nchark).eq.'iodif') then
          read(27,*,iostat=inperr) iodif
          write(*,*) 'iodif=',iodif
          if (iodif.eq.3) difloc=.true.
          if (iodif.eq.4) nonloc_tm86=.true.
          if (iodif.eq.5) nonloc_ls96=.true.
          if (iodif.eq.6) nonloc_noh03=.true.
          if (iodif.eq.7) nonloc_lock=.true.
          if (iodif.eq.8) loc_inm=.true.
        elseif(keyword(1:nchark).eq.'ifhle') then
          read(27,*,iostat=inperr) ifhle
          write(*,*) 'ifhle=',ifhle
        elseif(keyword(1:nchark).eq.'icetime') then
          read(27,*,iostat=inperr) icetime
          write(*,*) 'icetime=',icetime
        elseif(keyword(1:nchark).eq.'ug') then
          read(27,*,iostat=inperr) ug
          write(*,*) 'ug=',ug
        elseif(keyword(1:nchark).eq.'vg') then
          read(27,*,iostat=inperr) vg
          write(*,*) 'vg=',vg
        elseif(keyword(1:nchark).eq.'ts') then
          read(27,*,iostat=inperr) ts
          write(*,*) 'ts=',ts
        elseif(keyword(1:nchark).eq.'phi') then
          read(27,*,iostat=inperr) phi
          write(*,*) 'phi=',phi
        elseif(keyword(1:nchark).eq.'tfix') then
          read(27,*,iostat=inperr) tfix
          write(*,*) 'tfix=',tfix
        elseif(keyword(1:nchark).eq.'qif') then
          read(27,*,iostat=inperr) qif
          write(*,*) 'qif=',qif
        elseif(keyword(1:nchark).eq.'ifwr') then
          read(27,*,iostat=inperr) ifwr
          write(*,*) 'ifwr=',ifwr
        elseif(keyword(1:nchark).eq.'ifmf') then
          read(27,*,iostat=inperr) ifmf
          write(*,*) 'ifmf=',ifmf
        elseif(keyword(1:nchark).eq.'profile') then
          read(29,*,iostat=inperr) pts0,ioreft,iorefq,iorefu,iorefv
          read(29,*,iostat=inperr) ndth,ndus,ndvs,ndqvs
          ndat=ndth
          call allocref
          do ii=0,ndat
            read (29,*,iostat=inperr) zthdat(ii),thdat(ii)
     :        ,zusdat(ii),usdat(ii)
     :        ,zvsdat(ii),vsdat(ii)
     :        ,zvsdat(ii),qvsdat(ii)
            write(0,*)'thdat',ii,thdat(ii)
          enddo
          dzpro=zthdat(ndat)/npro
          pressdat(0)=pa
          write(*,*) 'Profile given'
        endif
        enddo
      end subroutine readpa
      
      subroutine allocref
      use alloc_1d
      
      if(allocated(zthdat)) deallocate(zthdat)
          allocate(zthdat(0:ndat))
          if(allocated(thdat)) deallocate(thdat)
          allocate(thdat(0:ndat))
          if(allocated(zusdat)) deallocate(zusdat)
          allocate(zusdat(0:ndat))
          if(allocated(usdat)) deallocate(usdat)
          allocate(usdat(0:ndat))
          if(allocated(zvsdat)) deallocate(zvsdat)
          allocate(zvsdat(0:ndat))
          if(allocated(vsdat)) deallocate(vsdat)
          allocate(vsdat(0:ndat))
          if(allocated(psdat)) deallocate(psdat)
          allocate(psdat(0:ndat))
          if(allocated(ptdat)) deallocate(ptdat)
          allocate(ptdat(0:ndat))
          if(allocated(zqvsdat)) deallocate(zqvsdat)
          allocate(zqvsdat(0:ndat))
          if(allocated(qvsdat)) deallocate(qvsdat)
          allocate(qvsdat(0:ndat))
          if(allocated(pressdat)) deallocate(pressdat)
          allocate(pressdat(0:ndat))
      
      end
