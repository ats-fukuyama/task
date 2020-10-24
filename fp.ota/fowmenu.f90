! fowmenu.f90

MODULE fowmenu

  PRIVATE
  PUBLIC fow_menu

CONTAINS

  SUBROUTINE fow_menu
    use fowcomm
    use fowprep
    use foworbit
    use foworbitclassify
    use fowdistribution

    USE fpcomm,ONLY: PM,FNS0
    USE fpinit,ONLY: fp_parm,fp_broadcast
    USE fpprep,ONLY: fp_prep
    use fpwrite
    USE fpsub,ONLY: FPMXWL
    use equnit_mod
    use libmpi

    implicit none

    CHARACTER(LEN=1)::  KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: MODE,IERR
    integer :: nr,nth,np,nsa
    real(rkind),allocatable :: fu(:,:,:,:),fI(:,:,:,:),J(:,:,:,:)

1   CONTINUE
    IF(nrank.EQ.0) THEN
       ierr=0
       WRITE(6,601)
601    FORMAT('## FOW MENU: R:RUN Q:QUIT')
       CALL TASK_KLIN(LINE,KID,MODE,fp_parm)
    ENDIF
    CALL mtx_barrier
    CALL mtx_broadcast_character(KID,1)
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.NE.1) GOTO 1


    SELECT CASE(kid)
    CASE('R')
       call fp_broadcast
       call fp_prep(IERR)

       call fow_read_namelist
       call fow_allocate
       call fow_prep

       call fow_orbit_construct(orbit_m)

       allocate(fu(npmax,nthmax,nrmax,nsamax), &
            fI(npmax,nthmax,nrmax,nsamax), &
            J(npmax,nthmax,nrmax,nsamax))
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
     
       do nr = 1, nrmax
          write(*,*)"FNS0",FNS0(1,1,NR,2)
       end do

       call fow_distribution_maxwellian_inCOM(fI)
       call fpcsv2D(fI(:,:,2,2),"./csv/fI_center.csv")
       call fpcsv2D(fI(:,:,nrmax,2),"./csv/fI_edge.csv")
       call fpcsv2D(fI(:,:,nrmax/2,2),"./csv/fI_quarter.csv")
       call fpcsv1D(fI(:,nthmax/4,nrmax/2,2),"./csv/fInth4.csv")
       call fpcsv1D(fI(:,nthmax/3,nrmax/2,2),"./csv/fInth3.csv")
       call fpcsv1D(fI(:,1,nrmax/2,2),"./csv/fInth0.csv")

    ! write(*,*)"1"
    ! call fow_orbit_jacobian(J,orbit_m)
    ! write(*,*)"2"
    ! call fpcsv2D(J(:,:,1,1),"./csv/J.csv")

       write(*,*)"debug"
       call fow_debug

       write(*,*)"end"
       call fow_deallocate
    CASE('Q')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE fow_menu

  subroutine fow_debug

  use fowcomm
  use foworbitclassify
  use fowdistribution

  use fpcomm
  use fpwrite

  implicit none

  integer :: IERR=0,nr,nth,np,nsa,nstp
  integer :: mode(3)
  real(rkind),allocatable :: check_orbit(:,:),lorentz(:),lorentzg(:),beta(:),betag(:),mean_r(:,:,:,:)&
                            ,trapped_boundary(:,:,:), forbitten_boundary(:,:,:), X_boundary(:,:,:)
  real(rkind) :: mean_psip

  write(*,*)"mean"
  allocate(mean_r(npmax,nthmax,nrmax,nsamax))
  do nsa = 1, nsamax
    do nr = 1, nrmax
      do nth = 1, nthmax
        do np = 1, npmax
          if ( forbitten(np,nth,nr,nsa,[0,0,0]) ) then
            mean_r(np,nth,nr,nsa) = 0.d0
          else
            write(*,*)1
            mean_psip = sum(orbit_m(np,nth,nr,nsa)%psip)/size(orbit_m(np,nth,nr,nsa)%psip)
            call fow_get_ra_from_psip(mean_r(np,nth,nr,nsa), mean_psip)
            mean_r(np,nth,nr,nsa) = rm(nr)-mean_r(np,nth,nr,nsa)
          end if
        end do
      end do
    end do
  end do

  call fpcsv2D(mean_r(:,:,2,2),"./csv/mean_r_center.csv")
  call fpcsv2D(mean_r(:,:,nrmax/2,2),"./csv/mean_r_quarter.csv")
  call fpcsv2D(mean_r(:,:,nrmax,2),"./csv/mean_r_edge.csv")

  write(*,*)"trapped_ion"
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

  write(*,*)"forbitten_ion"
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
  write(*,*)"trapped_ele"
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

  write(*,*)"forbitten_ele"
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

  write(*,*)"equiv_variable"
  call fpcsv1D(psimg,"./csv/psimg.csv")
  call fpcsv1D(Fpsig,"./csv/Fpsig.csv")
  call fpcsv1D(Boutg,"./csv/Boutg.csv")
  call fpcsv1D(Bing,"./csv/Bing.csv")
  call fpcsv1D(psim,"./csv/psim.csv")
  call fpcsv1D(Fpsi,"./csv/Fpsi.csv")
  call fpcsv1D(Bout,"./csv/Bout.csv")
  call fpcsv1D(Bin,"./csv/Bin.csv")
  call fpcsv1D(pm(:,1),"./csv/pm_ele.csv")
  call fpcsv1D(pg(:,1),"./csv/pg_ele.csv")
  call fpcsv1D(pm(:,2),"./csv/pm_ion.csv")
  call fpcsv1D(pg(:,2),"./csv/pg_ion.csv")
  call fpcsv1D(xi,"./csv/xi.csv")
  call fpcsv1D(xig,"./csv/xig.csv")

  allocate(trapped_boundary(nthmax,nrmax,nsamax),forbitten_boundary(nthmax,nrmax,nsamax))
  call fow_trapped_boundary(trapped_boundary)
  call fow_forbitten_boundary(forbitten_boundary)

  do nsa = 1, nsamax
    do nr = 1, nrmax
      do nth = 1, nthmax
        if ( trapped_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
          trapped_boundary(nth,nr,nsa) = pm(npmax,nsa)
        end if
        if ( forbitten_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
          forbitten_boundary(nth,nr,nsa) = pm(npmax,nsa)
        end if
      end do
    end do  
  end do

  call fpcsv2D(trapped_boundary(:,:,2),"./csv/trapped_boundary.csv")
  call fpcsv1D(trapped_boundary(:,2,2),"./csv/trapped_boundary_center.csv")
  call fpcsv1D(trapped_boundary(:,nrmax/2,2),"./csv/trapped_boundary_quarter.csv")
  call fpcsv1D(trapped_boundary(:,nrmax,2),"./csv/trapped_boundary_edge.csv")

  call fpcsv2D(forbitten_boundary(:,:,2),"./csv/forbitten_boundary.csv")
  call fpcsv1D(forbitten_boundary(:,2,2),"./csv/forbitten_boundary_center.csv")
  call fpcsv1D(forbitten_boundary(:,nrmax/2,2),"./csv/forbitten_boundary_quarter.csv")
  call fpcsv1D(forbitten_boundary(:,nrmax,2),"./csv/forbitten_boundary_edge.csv")

  ! allocate(X_boundary(nrmax,2,2))
  ! call fow_stagnation_type(X_boundary)
  
  ! do nr = 1, nrmax
  !   write(*,*)"-nr",nr,X_boundary(nr,2,1),X_boundary(nr,2,2)
  ! end do
  ! do nr = 1, nrmax
  !   write(*,*)"+nr",nr,X_boundary(nr,1,1),X_boundary(nr,1,2)
  ! end do
  ! write(*,*)"test"
  ! call fow_test

end subroutine fow_debug
END MODULE fowmenu
