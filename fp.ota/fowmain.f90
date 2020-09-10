program fow

  use fowcomm
  use fowprep
  use foworbit
  use foworbitclassify
  use fowdistribution

  use fpwrite

  use plinit
  use plparm
  use equnit_mod
  use fpinit
  use fpmenu
  use fpcomm
  use fpprep
  use fpsub

  use libmtx

  implicit none

  integer :: IERR=0,nr,nth,np,nsa,nstp
  integer :: mode(3)
  real(rkind),allocatable :: check_orbit(:,:),lorentz(:),lorentzg(:),beta(:),betag(:),fu(:,:,:,:),fI(:,:,:,:),J(:,:,:,:)

  call mtx_initialize

  CALL pl_init
  CALL eq_init
  CALL fp_init
  CALL pl_parm(1,'plparm',IERR)
  CALL eqparm(1,'eqparm',IERR)
  CALL fp_parm(1,'fpparm',IERR)

  call fp_broadcast
  call fp_prep(IERR)

  call fow_read_namelist
  call fow_allocate
  call fow_prep

  call fow_orbit_construct(orbit_m)

  allocate(fu(npmax,nthmax,nrmax,nsamax),fI(npmax,nthmax,nrmax,nsamax),J(npmax,nthmax,nrmax,nsamax))
  do nsa = 1, nsamax
    do nr = 1, nrmax
      do nth = 1, nthmax
        do np = 1, npmax
          fu(np,nth,nr,nsa)=FPMXWL(PM(NP,NSA),NR,NSA)
        end do
      end do
    end do
  end do
  call fpcsv2D(fu(:,:,1,1),"./csv/fu.csv")

  write(*,*)"1"
  call fow_orbit_jacobian(J,orbit_m)
  call fow_distribution_u2I(fI, fu, orbit_m)
  write(*,*)"2"
  call fpcsv2D(J(:,:,1,1),"./csv/J.csv")
  call fpcsv2D(fI(:,:,1,1),"./csv/fI.csv")

  nsa=2
  mode=[0,1,0]
  allocate(check_orbit(nthmax+mode(2),nrmax+mode(3)))
  do nr=1,nrmax+mode(3)
    do nth=1,nthmax+mode(2)
      do np=npmax+mode(1),1,-1
        if(trapped(np,nth,nr,2,[0,0,0]))then
          if(mode(1)==0)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)**2)
          if(mode(1)==1)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)**2)
          exit
        elseif(np==1.and.(.not.trapped(np,nth,nr,2,mode)))then
          check_orbit(nth,nr)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_trapped_ion.csv")
  deallocate(check_orbit)

  mode=[0,1,0]
  allocate(check_orbit(nthmax+mode(2),nrmax+mode(3)))
  do nr=1,nrmax+mode(3)
    do nth=1,nthmax+mode(2)
      do np=npmax+mode(1),1,-1
        if(.not.forbitten(np,nth,nr,2,mode))then
          if(mode(1)==0)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)**2)
          if(mode(1)==1)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)**2)
          exit
        elseif(np==1.and.(forbitten(np,nth,nr,2,mode)))then
          check_orbit(nth,nr)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_forbitten_ion.csv")
  deallocate(check_orbit)

  nsa=1
  mode=[0,1,0]
  allocate(check_orbit(nthmax+mode(2),nrmax+mode(3)))
  do nr=1,nrmax+mode(3)
    do nth=1,nthmax+mode(2)
      do np=npmax+mode(1),1,-1
        if(trapped(np,nth,nr,1,[0,0,0]))then
          if(mode(1)==0)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)**2)
          if(mode(1)==1)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)**2)
          exit
        elseif(np==1.and.(.not.trapped(np,nth,nr,1,mode)))then
          check_orbit(nth,nr)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_trapped_ele.csv")
  deallocate(check_orbit)

  mode=[0,1,0]
  allocate(check_orbit(nthmax+mode(2),nrmax+mode(3)))
  do nr=1,nrmax+mode(3)
    do nth=1,nthmax+mode(2)
      do np=npmax+mode(1),1,-1
        if(.not.forbitten(np,nth,nr,1,mode))then
          if(mode(1)==0)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)**2)
          if(mode(1)==1)check_orbit(nth,nr)=sqrt(1.d0-1.d0/sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)**2)
          exit
        elseif(np==1.and.(forbitten(np,nth,nr,1,mode)))then
          check_orbit(nth,nr)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_forbitten_ele.csv")

  call fpcsv1D(psimg,"./csv/psimg.csv")
  call fpcsv1D(Fpsig,"./csv/Fpsig.csv")
  call fpcsv1D(Boutg,"./csv/Boutg.csv")
  call fpcsv1D(Bing,"./csv/Bing.csv")
  call fpcsv1D(psim,"./csv/psim.csv")
  call fpcsv1D(Fpsi,"./csv/Fpsi.csv")
  call fpcsv1D(Bout,"./csv/Bout.csv")
  call fpcsv1D(Bin,"./csv/Bin.csv")

  deallocate(check_orbit)

  write(*,*)"end"
  call fow_deallocate

end program fow
