module eqread_mod
  use equ_params, only : eqfile
  implicit none
  real(8) :: sum_dl
  real(8), dimension(:), allocatable :: AJphVRL

  namelist /txequ/ eqfile

contains

!=======================================================================
!  Initialize equilibrium parameters
!=======================================================================
  subroutine eqinit

    call txequparf

  end subroutine eqinit


  SUBROUTINE TXEQUPARF
    use tx_commons, only : ierr
    integer(4) :: IST, KL
    logical :: LEX
    character(len=6) :: KPNAME = 'txparm'

    INQUIRE(FILE=KPNAME,EXIST=LEX)
    IF(.NOT.LEX) RETURN

    OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD')
    IF(IST > 0) THEN
       WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
       RETURN
    END IF
    READ(25,txequ,IOSTAT=IST)
    IF(IST > 0) THEN
       WRITE(6,*) 'XX PARM FILE READ ERROR'
       RETURN
    ELSE IF(IST < 0) THEN
       WRITE(6,*) 'XX PARM FILE EOF ERROR'
       RETURN
    END IF
    CALL KTRIM(KPNAME,KL)
    WRITE(6,*) &
         &     '## FILE (',KPNAME(1:KL),') IS ASSIGNED FOR TXEQU PARM INPUT'
    IERR=0

  END SUBROUTINE TXEQUPARF

!=======================================================================
!  Interface subroutine
!=======================================================================
  subroutine intequ
    use mod_spln
    use tx_interface, only : dfdx
    use tx_commons, only : ieqread, nrmax, nra, Pi, Pisq, rMU0, rr, ra, bb, rbvt, Rax, Zax, &
         & surflcfs, rho, rhov, epst, aat, rrt, ckt, suft, sst, vro, vlt, art, ait, bit, bbrt, &
         & elip, trig, PsitV, PsiV, hdt, fipol, sdt, bbt, rIPs, ft
    use equ_params, only : hiv, aav, rrv, ckv, shv, ssv, vlv, arv, aiv, biv, brv, epsv, elipv, &
         & trigv, qqv, siv, nv, rmaj, rpla, raxis, zaxis, pds, fds, sdw, ftv
!    use tx_core_module, only : intg_area
    integer(4) :: n, nr, ierr, nrmaxx, nrax, nrs
    real(8), parameter :: fourPisq = 4.d0 * Pi * Pi
    real(8) :: rhonrs, asdt, bsdt, sdtvac
    real(8), allocatable :: U(:,:), U0(:), deriv(:), rho_v(:), fipolv(:), hdv(:), zzv(:), &
         &                  pdst(:), fdst(:), qqt(:), zz(:), zzfunc(:)
    real(8) :: aitken2p

    nrmaxx = nrmax + 1
    nrax   = nra   + 1

    rhonrs = 0.1d0
    do nr = 0, nrmax
       if( rho(nr) > rhonrs ) then
          nrs = nr - 1
          exit
       end if
    end do

    ! Read equilibrium data in MEUDAS/SELENE format.

    call eqdsk

    allocate(hdv(nv),fipolv(nv))

    ! hdv: dpsit/dV = q * dpsi/dV = qqv * sdw
    !    qqv and sdw from MEUDAS/SELENE often diverge near the magnetic axis
    !    and anomaly behavior of qqv stems from that of sdw.
    !    Therefore, hdv = qqv * sdw may not have odd behavior near the magnetic axis.

    hdv(1:nv) = qqv(1:nv) * sdw(1:nv)

    ! Exclude anomaly at axis and compensate it by linear extrapolation

    hdv(1) = aitken2p(vlv(1),hdv(2),hdv(3),hdv(4),vlv(2),vlv(3),vlv(4))

    ! hiv = toroidal magnetic flux on the equilibrium mesh
    ! hiv is calculated by integration with spline
    allocate(U(4,nv),U0(nv),deriv(nv))

    call SPL1D  (vlv,hdv,deriv,U,nv,0,ierr)
    call SPL1DI0(vlv,U,U0,nv,ierr)
    do n = 1, nv
       call SPL1DI(vlv(n),hiv(n),vlv,U,U0,nv,ierr)
    end do

    deallocate(U,U0,deriv)

    ! ******************************************************************
    !     fipol: R B_t,  Poloidal current function
    !        I = 4 Pi^2 * dpsit/dV / <R^{-2}>
    !          = fourPisq * hdv / aav
    !
    !     fipol is at first calculated so that PsitV is calculated
    !        based on fipol from the aspect of consistency, because
    !        in txcalv fipol is always calculated by differentiating PsitV.
    ! ******************************************************************

    ! fipol and rho on the equilibrium grid

    fipolv(1:nv) = fourPisq * hdv(1:nv) / aav(1:nv) ! equ. grid

    allocate(rho_v(nv),zzv(nv),zz(0:nrmax))
    rho_v(1:nv) = sqrt( hiv(1:nv) / hiv(nv) ) ! equ. grid

    ! fipol: poloidal current function, assumed to be constant outside the plasma initially
    call spln( fipol, rho,   nrax,   fipolv,rho_v,    nv,   1 )
    fipol(nrax:nrmax) = fipol(nra)

    ! epst: inverse aspect ratio, proportional to rho: epst = (Rmax-Rmin)/(Rmax+Rmin)
    call spln( epst,  rho,   nrmaxx, epsv,  rho_v,    nv, 151 )
    ! suft: (interim definition) <|nabla psi|>, 
    !      proportional to rho near the axis
    !      proportional to rho^{-1} outside the plasma surface
    call spln( suft(0:nra), rho(0:nra), nrax, shv, rho_v, nv, 1 ) ! shv = <|nabla psi|>

!    call spln( vlt,   rho,  nrmaxx,  vlv,   rho_v, nv, 122 )
!    call spln( aat,   rho,  nrmaxx,  aav,   rho_v, nv, 121 )
    ! vlt: volume, proportional to rho**2
    zzv(:)    = rho_v(:)**2
    zz(:)     = rho(:)**2
    call spln( vlt,   zz,   nrmaxx,  vlv,   zzv,   nv, 151 )
    ! aat: <R^-2>, proportional to (sqrt(1-epst**2))^{-1}
    zzv(1:nv) = 1.d0 / sqrt( 1.d0 - epsv(1:nv)**2 )
    zz(:)     = 1.d0 / sqrt( 1.d0 - epst(:)**2    )
    call spln( aat,   zz,   nrmaxx,  aav,   zzv , nv, 151 )

    ! hdt: dpsit/dV
    hdt(:) = fipol(:) * aat(:) / fourPisq

    allocate(U(4,0:nrmax),U0(0:nrmax),deriv(0:nrmax),qqt(0:nrmax),zzfunc(0:nrmax))

    ! === Volume integration of I<R^{-2}> to obtain psit ===

    call SPL1D  (vlt,hdt,deriv,U,nrmaxx,0,ierr)
    call SPL1DI0(vlt,U,U0,nrmaxx,ierr)
    do nr = 0, nrmax
       call SPL1DI(vlt(nr),PsitV(nr),vlt,U,U0,nrmaxx,ierr)
    end do

    !  === ckt(NRMAX) is related to sdt via plasma current ===

    call spln( ckt,   vlt,   nrmaxx, ckv,   vlv, nv, 151 )

    ! === Volume integration of I<R^{-2}>/q to obtain psi ===

    ! qqt: safety factor from equilibrium data (not transferred to TX)
    call spln( qqt(0:nra),PsitV(0:nra), nrax, qqv,   hiv, nv, 1 )
    ! sdt: dpsi/dV = (dpsit/dV) / q
    sdt(0:nra) = hdt(0:nra) / qqt(0:nra)

    sdtvac = 2.d0 * Pi * rMU0 * rIPs * 1.D6 / ckt(nrmax)

!!$    ! interpolate sdt between sdt(nra) at rho=rhoa and sdtvac at rho=rhob
!!$    ! sdt should be proportional to epst^{-2} because outside the plasma surface
!!$    !    sdt = dpsi/dV = <|nabla psi|>/<|nabla V|> = epst^{-1}/epst = epst^{-2}
!!$    asdt = (epst(nra)*epst(nrmax))**2*(sdt(nra)-sdtvac)  / (epst(nrmax)**2-epst(nra)**2)
!!$    bsdt = (sdtvac*epst(nrmax)**2-sdt(nra)*epst(nra)**2) / (epst(nrmax)**2-epst(nra)**2)
!!$    do nr = nrax, nrmax
!!$       sdt(nr) = asdt / epst(nr)**2 + bsdt
!!$    end do

    zzfunc(0:nra) = ckt(0:nra) * sdt(0:nra)**2
    zzfunc(nrax)  = ckt(nrmax) * sdtvac**2
    zz    (0:nra) = rho(0:nra)
    zz    (nrax)  = rho(nrmax)
    call spln( zzfunc, rho, nrmaxx, zzfunc(0:nrax), zz(0:nrax), nrax+1, 151 )
    sdt(nrax:nrmax) = sqrt( zzfunc(nrax:nrmax) / ckt(nrax:nrmax) )

    call SPL1D  (vlt,sdt,deriv,U,nrmaxx,0,ierr)
    call SPL1DI0(vlt,U,U0,nrmaxx,ierr)

    do nr = 0, nrmax
       call SPL1DI(vlt(nr),PsiV(nr),vlt,U,U0,nrmaxx,ierr)
    end do

    ! suft: <|nabla V|> = <|nabla psi|> / sdt
    !       <|nabla V|> is proportional to epst

    suft(0:nra) = suft(0:nra) / sdt(0:nra)
    call spln( suft, epst, nrmaxx, suft(0:nra), epst(0:nra), nrax, 151 )

!!$    do nr=0,nrmax
!!$       write(201,'(F8.5,1P8E15.7)') rho(nr),PsitV(nr)/PsitV(nra),vlt(nr)/vlt(nra),epst(nr),suft(nr),sdt(nr),ckt(nr),ckt(nr)*sdt(nr)**2
!!$    end do

    deallocate(U,U0,deriv,hdv,fipolv,qqt,zzfunc)

    ! ******************************************************************
    !    Mapping & extrapolate
    ! ******************************************************************

    allocate(pdst(0:nrmax),fdst(0:nrmax),zzfunc(1:nv))

    call spln( rrt,  PsitV, nrmaxx, rrv,   hiv, nv, 151 )
    call spln( sst,  PsitV, nrmaxx, ssv,   hiv, nv, 151 ) ! ssv = <|nabla V|^2>
    call spln( art,  PsitV, nrmaxx, arv,   hiv, nv, 152 )
    call spln( ait,  PsitV, nrmaxx, aiv,   hiv, nv, 151 )
    call spln( bit,  PsitV, nrmaxx, biv,   hiv, nv, 151 )
    call spln( bbrt, PsitV, nrmaxx, brv,   hiv, nv, 151 )
    call spln( elip, PsitV, nrmaxx, elipv, hiv, nv, 151 )
    call spln( trig, PsitV, nrmaxx, trigv, hiv, nv, 151 )
    call spln( pdst, PsitV, nrmaxx, pds,   hiv, nv,   0 ) ! Not extrapolated
    call spln( fdst, PsitV, nrmaxx, fds,   hiv, nv,   0 ) ! Not extrapolated

    ! ft: Trapped particle fraction (ft^4 almost linearly scales with V or PsitV.)

    zzfunc(1:nv) = ftv(1:nv)**4
    call spln( ft,   PsitV, nrmaxx, zzfunc,hiv, nv, 151 )
    ft(:) = sqrt(sqrt(ft(:)))
    where( ft > 1.d0 ) ft = 1.d0

!!$    ! ******************************************************************
!!$    !     Replace profiles near the axis by Least Square Method
!!$    ! 
!!$    !   Due to coarse grids near the axis on the equilibrium grid, 
!!$    !     mapping profiles on the equilibrium grid onto the transport grid
!!$    !     may gives rise to the shape of profiles different from what is
!!$    !     supposed to be.
!!$    !   To avoid this, LSM is applied to the variables whose ideal profile
!!$    !     has already been known near the axis.
!!$    !   For example, the volume (vlt) is known to behave a function of rho^2.
!!$    ! ******************************************************************
!!$
!!$    !                mode,exp,coord_in,      val_in,   coord_out,     val_out
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),   rrt(0:nrs), rho(0:nrs),   rrt(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),   sst(0:nrs), rho(0:nrs),   sst(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),   art(0:nrs), rho(0:nrs),   art(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),   ait(0:nrs), rho(0:nrs),   ait(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),   bit(0:nrs), rho(0:nrs),   bit(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),  bbrt(0:nrs), rho(0:nrs),  bbrt(0:nrs))
!!$!    call lesq6_wrapper(2, 2, rho(0:nrs),  PsiV(0:nrs), rho(0:nrs),  PsiV(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),  pdst(0:nrs), rho(0:nrs),  pdst(0:nrs))
!!$    call lesq6_wrapper(2, 2, rho(0:nrs),  fdst(0:nrs), rho(0:nrs),  fdst(0:nrs))

    ! Rescale

    PsiV(:)  = PsiV(:) - PsiV(0)

    ! vro: dV/drho

    vro(0) = 0.d0
    vro(1:nrmax) = 2.d0 * PsitV(nrmax) * rho(1:nrmax) / hdt(1:nrmax)

    ! <B^2>
    
    bbt(:) = fourPisq*fourPisq * hdt(:)*hdt(:) / aat(:) + ckt(:) * sdt(:)*sdt(:)

    ! Replace geometric parameters given by namelist or default
    !   with those from an equilibrium

    rr   = rmaj
    ra   = rpla
    bb   = rbvt / rr
    Rax  = raxis
    Zax  = zaxis
    surflcfs = sum_dl * 2.d0 * Pi * rr

    ! Mesh

    call txmesh

    ! ******************************************************************

    ! Accuracy of elongation and triangularity near the axis significantly
    !   depends upon the number (or the size) of the radial mesh.
    ! Avoiding their non-physical values requires the replacement of the values
    !   inside rho~0.075 (when ivdm=201) 
    !   inside rho~0.15  (when ivdm= 51)

    do nr = nrmax, 0, -1
       if( rho(nr) < 0.1d0 ) then
          nrs = nr + 1
          elip(0:nrs-1) = elip(nrs)
!          trig(0:nrs-1) = trig(nrs)
          exit
       end if
    end do

    ! <j_phi/R>: toroidal current density
    ! electron current: AJphVRL = - e ne <ueph/R>

    if(.not. allocated(AJphVRL)) allocate(AJphVRL(0:NRMAX))
    AJphVRL(0:NRA) = - ( aat(0:NRA) * fdst(0:NRA) + pdst(0:NRA) ) / rMU0
    AJphVRL(NRA+1:NRMAX) = 0.d0

!!$    do nr=0,nrmax
!!$       write(200,'(F8.5,1P21E15.7)') rho(nr), &
!!$            & vlt(nr)/vlt(nra),fipol(nr),epst(nr),suft(nr),aat(nr),hdt(nr),psitv(nr),ckt(nr),sdt(nr),psiv(nr), &
!!$            & rrt(nr),sst(nr),art(nr),ait(nr),bit(nr),bbrt(nr),elip(nr),trig(nr),vro(nr),bbt(nr), &
!!$            & AJphVRL(nr)
!!$    end do

    ! Confirm total plasma current [A]
    !write(6,*) intg_area(ajphvrl/ait)

    deallocate(rho_v,pdst,fdst,zzv,zz,zzfunc)

  end subroutine intequ

!=======================================================================
!  Mesh 
!=======================================================================

  subroutine txmesh
    use tx_commons, only : r, ra, rb, rr, ravl, rbvl, vlt, Pisq, rho, rhov, rhosq, vv, hv, &
         &                 NEMAX, NRMAX, NRA, rhob
    integer(4) nr

    ! rhov = volume rho [m]

    !--- Create various kinds of mesh that depends upon equilibrium ---
  
    rhosq(:) = rho(:)**2
    rhov(:)  = sqrt( vlt(:) / ( 2.d0 * Pisq * rr ) )

    vv(:)  = vlt(:)

    !  Mesh interval

    hv(1:NEMAX) = vlt(1:NRMAX) - vlt(0:NRMAX-1)

    !------------------------------------------------------------------

    ! ravl: volume-averaged minor radius at the plasma surface

    ravl = rhov(nra)

    r(:)     = rho(:) * ravl

    ! rbvl: volume-averaged minor radius at the virtual wall

    rbvl = rhov(nrmax)

    !   Virtual wall radius (m)

    rb   = rhob * ra   ! geometric

  end subroutine txmesh

!=======================================================================
! Allocate and deallocate arrays associated with equilibria
!=======================================================================

  subroutine alloc_equ(mode)
    use equ_params

    integer(4), intent(in) :: mode

    select case(mode)
    case(1)
       ! *** Allocate arrays that are widely used. ***
       allocate(rg(irdm),zg(izdm2),psi(irzdm2),rbp(irzdm2))
       allocate(pds(ivdm),fds(ivdm),vlv(ivdm),qqv(ivdm),prv(ivdm))
       allocate(csu(isrzdm),rsu(isrzdm),zsu(isrzdm))
       allocate(hiv(ivdm),siv(ivdm),siw(ivdm),sdw(ivdm),ckv(ivdm),ssv(ivdm),aav(ivdm),rrv(ivdm), &
            &   rbv(ivdm),arv(ivdm),bbv(ivdm),biv(ivdm),r2b2v(ivdm),shv(ivdm),grbm2v(ivdm), &
            &   rov(ivdm),aiv(ivdm),brv(ivdm),epsv(ivdm),elipv(ivdm),trigv(ivdm),ftv(ivdm))
       
    case(2)
       ! *** Allocate arrays that are used only when reading an equilibium data. ***
       ! *** They are specific to MEUDAS/SELENE.                                 ***
       allocate(ieqout(10),ieqerr(10),icp(10),cp(10))
       allocate(ivac(0:nsfix),ncoil(0:nsfix),cvac(0:nsfix),rvac(0:nsfix),zvac(0:nsfix))
       allocate(rcoil(100,icvdm),zcoil(100,icvdm),ccoil(100,icvdm),rlimt(200),zlimt(200))

    case(-1)
       ! *** Deallocate arrays that are widely used ***
       if(allocated(rg)) then
          deallocate(rg,zg,psi,rbp)
          deallocate(pds,fds,vlv,qqv,prv)
          deallocate(csu,rsu,zsu)
          deallocate(hiv,siv,siw,sdw,ckv,ssv,aav,rrv, &
               &     rbv,arv,bbv,biv,r2b2v,shv,grbm2v,rov, &
               &     aiv,brv,epsv,elipv,trigv,ftv)
       end if

    case(-2)
       ! *** Deallocate arrays that are used only when reading an equilibium data ***
       deallocate(ieqout,ieqerr,icp,cp)
       deallocate(ivac,ncoil,cvac,rvac,zvac)
       deallocate(rcoil,zcoil,ccoil,rlimt,zlimt)

    end select

  end subroutine alloc_equ

!=======================================================================
!             eqdsk               level=32       date=81.12.23
!=======================================================================
  subroutine eqdsk
!=======================================================================
    use tx_commons, only : cnpi => PI, rMU0
    use equ_params
    logical :: lex
    integer(4) :: i, j, n, ir, iz, ist, igeqdsk, ieqtmp = 22
    real(8) :: rsep, zsep, rrmax, rrmin, zzmax, zzmin, rzmax, rzmin, rpmax, rpmin, betp, zzlam
    real(8) :: dpsi, btv2, dsr, dsz, dxx, dl
    character(len=100) :: command, tmpfile, content
!=======================================================================
    igeqdsk = 0

    inquire(file=eqfile,exist=lex)
    if( lex .neqv. .true. ) stop 'No equilibrium file.'
!-----------------------------------------------------------------------

    call alloc_equ(1)
    call alloc_equ(2)

    ! check eqfile with the format of whether binary (MEUDAS) or ascii (G EQDSK)
    tmpfile = trim(eqfile)//".tmp"
    command = "file "//trim(eqfile)//" > "//trim(tmpfile)
    call system(trim(command))
    open(ieqtmp,file=tmpfile,iostat=ist,form='formatted')
    if(ist /= 0) stop 'file open error.'
    read(ieqtmp,'(A)') content
    i = index(content,'ASCII')
    close(ieqtmp,status='delete')

    if( i == 0 ) then
       ! *** Standard MEUDAS/SELENE format (unformatted) ***

       open(ieqrd,file=eqfile,iostat=ist,form='unformatted',action='read')
       if(ist /= 0) stop 'file open error.'

       read(ieqrd,iostat=ist)nsr,nsz  &
            &        ,(rg(i),i=1,nsr),(zg(j),j=1,nsz)  &
            &        ,(psi(i),i=1,nsr*nsz)  &
            &        ,nv,(pds(i),i=1,nv),(fds(i),i=1,nv)  &
            &        ,(vlv(i),i=1,nv),(qqv(i),i=1,nv),(prv(i),i=1,nv)  &
            &        ,btv,ttcu,ttpr,bets,beta,betj  &
            &        ,(icp(i),i=1,10),(cp(i),i=1,10)  &
            &        ,saxis,raxis,zaxis,ell,trg  &
            &        ,nsu,(rsu(i),i=1,nsu),(zsu(i),i=1,nsu),(csu(i),i=1,nsu)  &
            &        ,isep,dsep,rsep,zsep  &
            &        ,(ivac(i),i=0,nsfix)  &
            &        ,(rvac(i),i=0,nsfix),(zvac(i),i=0,nsfix)  &
            &        ,(cvac(i),i=0,nsfix),(ncoil(i),i=1,nsfix)  &
            &        ,((rcoil(i,j),i=1,100),j=1,icvdm)  &
            &        ,((zcoil(i,j),i=1,100),j=1,icvdm)  &
            &        ,((ccoil(i,j),i=1,100),j=1,icvdm)  &
            &        ,ilimt,(rlimt(i),i=1,100),(zlimt(i),i=1,100)
       close(ieqrd)
       if( ist > 0 ) then ! error
          stop '===== dataio/read error ====='
       else if( ist < 0 ) then ! end
          backspace ieqrd
          write(6,"(5x,'dataio/end:backspace')")
          read(*,*)
       end if
    else
       ! *** Standard EFIT format, G EQDSK (formatted) ***
       !   We assume that any ascii file is regarded as G EQDSK format files.

       open(ieqrd,file=eqfile,iostat=ist,form='formatted',action='read')
       if(ist /= 0) stop 'file open error.'
       call read_geqdsk(igeqdsk)
       close(ieqrd)
    end if
!-----------------------------------------------------------------------
    write(6,40)eqfile,ieqrd,jeqrd
40  format(//5x,'diskio:',a19,':','read(',i3,'/',i3,')')
!-----------------------------------------------------------------------
!=======================================================================
    sigcu=1.d0
      iudsym=0
    nr=nsr
      nrm=nr-1
      nz=(nsz-1)/2+1
      nszm=nsz-1
      nvm=nv-1
!-----	  
    dr=rg(2)-rg(1)
    dz=zg(2)-zg(1)
    dr2i=1.d0/(2.d0*dr)
    dz2i=1.d0/(2.d0*dz)
    ddri=1.d0/(dr*dr)
    ddzi=1.d0/(dz*dz)
!-----
    iraxis=(raxis-rg(1))/dr+1.d0
    izaxis=(zaxis-zg(1))/dz+1.d0
    iaxis=nr*(izaxis-1)+iraxis
!-----
!    call eqrbp
!    call eqlin
!-----
    qaxis=qqv(1)
    qsurf=qqv(nv)
    dpsi=-saxis/dfloat(nv-1)
    btv2=btv*btv
    rbv(nv)=btv
    do n=nv-1,1,-1
       btv2=btv2-dpsi*(fds(n+1)+fds(n))
       rbv(n)=sqrt(btv2)
    end do
    do iz=2,nsz-1
       do ir=2,nsr-1
          i=nsr*(iz-1)+ir
          dsr=(psi(i+1)-psi(i-1))/(2.d0*dr)
          dsz=(psi(i+nsr)-psi(i-nsr))/(2.d0*dz)
          rbp(i)=sqrt(dsr*dsr+dsz*dsz)
       end do
    end do
    dxx=dpsi*0.5d0*(2.d0*cnpi)**2
    hiv(1)=0.d0
    do n=2,nv
       hiv(n)=hiv(n-1)+dxx*(qqv(n)+qqv(n-1))
    end do
!-----
    call eqrbp
    call eqlin
    do n = nv, 1, -1
       siv(n) = siw(n)
    end do
!-----
    cpl=0.d0
    sum_dl=0.d0
    rrmax=rsu(1)   ! rrmax = R_max
    rrmin=rsu(1)   ! rrmin = R_min
    zzmax=zsu(1)   ! zzmax = Z_max
    zzmin=zsu(1)   ! zzmin = Z_min
    rzmax=-1000.d0 ! rzmax = R at Z_max
    rzmin= 1000.d0 ! rzmin = R at Z_min
    do i=2,nsu
       dl=sqrt((rsu(i)-rsu(i-1))**2+(zsu(i)-zsu(i-1))**2)
       sum_dl=sum_dl+dl
       cpl=cpl+0.5d0*(csu(i)+csu(i-1))*dl
       if(rsu(i) > rrmax)rrmax=rsu(i)
       if(rsu(i) < rrmin)rrmin=rsu(i)
       if(zsu(i) > zzmax)then
          zzmax=zsu(i)
          rzmax=rsu(i)
       endif
       if(zsu(i) < zzmin)then
          zzmin=zsu(i)
          rzmin=rsu(i)
       endif
    end do
    ccpl=cpl/(rMU0*1.d6)
    rpmax=rrmax
    rpmin=rrmin
    rmaj=0.5d0*(rpmax+rpmin)
    rpla=0.5d0*(rpmax-rpmin)
!--
    elipup=(zzmax-zaxis)/rpla
    trigup=(rmaj -rzmax)/rpla
    elipdw=(zaxis-zzmin)/rpla
    trigdw=(rmaj -rzmin)/rpla
    elip=0.5d0*(elipup+elipdw)
    trig=0.5d0*(trigup+trigdw)
!-----
    qqj=5.d0*rpla**2*btv/rmaj**2/ttcu*(1.d0+elip*2)/2.d0
!-----------------------------------------------------------------------
!!$    psi0=0.05d0*saxis
!!$    call eqtrc(psi0,raxis,-1,irend,izend)
!!$    q95=1.d0/((2.d0*zpi)**2)*btv*zaav/zsdw
!!$!-----
!!$    psi0=0.05d0*saxis
!!$    x=nvm*(saxis-psi0)/saxis+1.d0
!!$    ix=x
!!$    if(ix.ge.nv)ix=nvm
!!$    x=x-dfloat(ix)
!!$    q95=qqv(ix)+(qqv(ix+1)-qqv(ix))*x
!
!      bets=2.d0*sqrt(tps/vlv(nv))/(rbv(nv)/rmaj)**2
!      beta=2.d0*tpr/(vlv(nv)*(rbv(nv)/rmaj)**2)
!      betj=4.d0*tpr/(rmaj*cpl**2)
!      betp=2.d0*tpr/vlv(nv)/(cpl/zlen)**2
!-----To avoid error
    betp=0.d0
    zzlam=0.d0
!-----
    do n = 1, nv
       prv(n) = prv(n) / (rMU0*1.d6)
!-----Minor radius
       rov(n) = sqrt( vlv(n) / ( 2.0_8 * cnpi**2 * rmaj ) )
    enddo

    if( igeqdsk /= 0 ) ccpl = ttcu

!=======================================================================
    write(6,12) betj,beta,bets,ttpr,prfac  &
         &     ,ccpl,ttcu,btv,betp,qaxis,qsurf,q95,qqj  &
         &     ,raxis,zaxis,saxis  &
         &     ,rmaj,rpla,rov(nv),rpmax,rpmin  &
         &     ,vlv(nv),elip,trig,el95  &
         &     ,elipup,trigup,elipdw,trigdw  &
         &     ,zzlp,zzli,zzlam  &
         &     ,(cvac(i),i=0,3)
    write(6,13)(i,cvac(i),i=4,icvdm)
    write(6,*)
12  format(//10x,'===== equilibrium ===== '//  &
           & 10x,'beta-j=',f10.5,2x,'beta-t=',f10.5,2x,'beta-s=',f10.5  &
           &, 2x,'tpress=',f10.5,2x,'p0/pav=',f10.5  &
           &/10x,'cpl   =',f10.5,2x,'tcur  =',f10.5,2x,'btv   =',f10.5  &
           &, 2x,'beta-p=',f10.5  &
           &/10x,'qaxis =',f10.5,2x,'qsurf =',f10.5  &
           & ,2x,'q95   =',f10.5,2x,'qqj   =',f10.5  &
           &/10x,'raxis =',f10.5,2x,'zaxis =',f10.5,2x,'saxis =',f10.5  &
           &/10x,'rmaj  =',f10.5,2x,'rpla  =',f10.5  &
           & ,2x,'apave =',f10.5,2x,'rpmax =',f10.5,2x,'rpmin =',f10.5  &
           &/10x,'volume=',f10.5,2x,'ellip =',f10.5,2x,'trig  =',f10.5  &
           & ,2x,'el95  =',f10.5  &
           &/10x,'elipup=',f10.5,2x,'trigup=',f10.5,2x,'elipdw=',f10.5  &
           & ,2x,'trigdw=',f10.5  &
           &/10x,'lp    =',f10.5,2x,'li    =',f10.5,2x,'lam   =',f10.5  &
           &/10x,'s-surf=',f10.5,2x,'coil 1=',f10.5,2x,'coil 2=',f10.5  &
           & ,2x,'coil 3=',f10.5)
13  format(10x,'coil',i2,'=',f10.5,2x,'coil',i2,'=',f10.5  &
         & ,2x,'coil',i2,'=',f10.5,2x,'coil',i2,'=',f10.5)
!-----------------------------------------------------------------------
    if( isep /= 0 )then
       write(6,14)rsep,zsep,rspmx,rspmn,psep
14     format(/  &
            & 10x,'rsep  =',f10.5,2x,'zsep  =',f10.5,2x,'rspmx ='  ,f10.5  &
            & ,2x,'rspmn =',f10.5,2x,'psep  =',f10.5/)
    endif
!=======================================================================

    call alloc_equ(-2)

  end subroutine eqdsk

!=======================================================================
!             eqrbp               revised on oct.,1986
!=======================================================================
  subroutine eqrbp
!=======================================================================
!     load rbp                                                     jaeri
!=======================================================================
    use equ_params
    integer(4) :: i, ir, iz
    real(8) :: x, y
!-----------------------------------------------------------------------
    i=nrm
    do iz=2,nszm
       i=i+2
       do ir=2,nrm
          i=i+1
          x=dr2i*(psi(i+ 1)-psi(i- 1))
          y=dz2i*(psi(i+nr)-psi(i-nr))
          rbp(i)=sqrt(x*x+y*y)
       end do
    end do
  end subroutine eqrbp

!=======================================================================
!             eqlin               revised        nov.,1986
!=======================================================================
  subroutine eqlin
!=======================================================================
!     line integrals
!=======================================================================
    use tx_commons, only : zpi => PI
    use equ_params
    ! intf: num. of division in the direction of lambda for trapped particle fraction
    integer(4), parameter :: intf = 100
    integer(4) :: i,i1,i2,i3,i4,ir,iz,istep,irst,ist,jr,k,ll,lm,lp,n,j
    real(8) :: dpsi,psi0,x,r0,z0,r1,z1,s1,s2,s3,s4,dl &
         &     ,bp0,bl0,ds0,ck0,ss0,vl0,aa0,rr0,bb0,bi0,sh0 &
         &     ,bp1,bl1,ds1,ck1,ss1,vl1,aa1,rr1,bb1,bi1,sh1 &
         &     ,ai0,ai1,br0,br1,bm0,bm1,zzmax,zzmin,rzmax,rzmin &
         &     ,rrmax,rrmin,rmajl,rplal,sdw2
    real(8) :: fintx, hsq, h
    integer(4), dimension(:), allocatable :: nsul
    real(8), dimension(:), allocatable :: bmax,fint,flam,dll,zbl
!=======================================================================
    ieqerr(1)=0
    allocate(bmax(ivdm),fint(0:intf),flam(0:intf))
    allocate(nsul(isrzdm),dll(isrzdm),zbl(isrzdm))
    do j = 0, intf
       flam(j)  = real( j, 8 ) / real( intf, 8 )
    end do
!-----------------------------------------------------------------------
    istep=0
999 nsu=0
    dpsi=saxis/float(nvm)
    do jr=iraxis,nrm
       irst=jr
       ist= (izaxis-1)*nr + jr
       if(sigcu*psi(ist+1) > 0.d0) exit
    end do
!-----------------------------------------------------------------------
    do n=nv,2,-1
!-----------------------------------------------------------------------
!..for double integral to compute trapped particle fraction
       nsul(n)=0
       fint(:) = 0.d0
!-----------------------------------------------------------------------
       psi0=dpsi*dfloat(nv-n)
5      siw(n)=psi0
!-----------------------------------------------------------------------
!..search starting point
       i=ist+1
       ir=irst+1
       iz=izaxis
10   i=i-1
       ir=ir-1
       if(ir < iraxis) go to 800
       if((psi0-psi(i))*(psi0-psi(i+1)) > 0.d0) go to 10
       ist=i
       irst=ir
!-----------------------------------------------------------------------
!..initialize line integrals
       arv(n)=0.d0
       vlv(n)=0.d0
       sdw(n)=0.d0
       ckv(n)=0.d0
       ssv(n)=0.d0
       aav(n)=0.d0
       rrv(n)=0.d0
       bbv(n)=0.d0
       biv(n)=0.d0
       shv(n)=0.d0
       grbm2v(n)=0.d0
       aiv(n)=0.d0
       brv(n)=0.d0
       bmax(n)=0.d0

       x=(psi0-psi(i))/(psi(i+1)-psi(i))
       r1=rg(ir)+x*dr
       z1=zg(iz)
       bp1=(rbp(i)+x*(rbp(i+1)-rbp(i)))/r1
       vl1=r1*r1
       ds1=1.d0/bp1
       ck1=bp1
       ss1=vl1*bp1
       aa1=ds1/vl1
       rr1=vl1/bp1
       bb1=rbv(n)*rbv(n)/(r1*r1)+bp1*bp1
       bi1=ds1/bb1
       bm1=sqrt(bb1)
       br1=bm1*ds1
       bb1=bb1*ds1
       sh1=r1
       ai1=ds1/r1

       rrmax=r1
       rrmin=r1
       zzmax=z1
       zzmin=z1
       rzmax=-1000.d0
       rzmin= 1000.d0

!..trace contour
       k=1
       nsul(n)=1
       if(n == nv) then
          nsu=1
          rsu(1)=r1
          zsu(1)=z1
          csu(1)=bp1*sigcu
       endif
20     r0=r1
       z0=z1
       bp0=bp1
       vl0=vl1
       ds0=ds1
       ck0=ck1
       ss0=ss1
       aa0=aa1
       rr0=rr1
       bi0=bi1
       bb0=bb1
       sh0=sh1
       ai0=ai1
       br0=br1
       bm0=bm1
       i1=i
       i2=i1+1
       i3=i2-nr
       i4=i1-nr
       s1=psi0-psi(i1)
       s2=psi0-psi(i2)
       s3=psi0-psi(i3)
       s4=psi0-psi(i4)
       if(s1*s2 < 0.d0 .and. k /= 1) go to 30
       if(s2*s3 < 0.d0 .and. k /= 2) go to 40
       if(s3*s4 < 0.d0 .and. k /= 3) go to 50
       if(s4*s1 < 0.d0 .and. k /= 4) go to 60
       psi0=psi0+1.d-08*saxis
       go to 5
30     x=s1/(s1-s2)
       r1=rg(ir)+dr*x
       z1=zg(iz)
       bp1=(rbp(i1)+x*(rbp(i2)-rbp(i1)))/r1
       i=i+nr
       iz=iz+1
       k=3
       go to 70
40     x=s2/(s2-s3)
       r1=rg(ir+1)
       z1=zg(iz)-dz*x
       bp1=(rbp(i2)+x*(rbp(i3)-rbp(i2)))/r1
       i=i+1
       ir=ir+1
       k=4
       go to 70
50     x=s4/(s4-s3)
       r1=rg(ir)+dr*x
       z1=zg(iz-1)
       bp1=(rbp(i4)+x*(rbp(i3)-rbp(i4)))/r1
       i=i-nr
       iz=iz-1
       k=1
       go to 70
60     x=s1/(s1-s4)
       r1=rg(ir)
       z1=zg(iz)-dz*x
       bp1=(rbp(i1)+x*(rbp(i4)-rbp(i1)))/r1
       i=i-1
       ir=ir-1
       k=2
!..check
70     if( ir <=  1   ) go to 810
       if( ir >= nrm  ) go to 810
       if( iz <=  1   ) go to 810
       if( iz >= nszm ) go to 810
!..line integrals
       vl1=r1*r1
       ds1=1.d0/bp1
       ck1=bp1
       ss1=vl1*bp1
       aa1=ds1/vl1
       rr1=vl1/bp1
       bb1=rbv(n)*rbv(n)/(r1*r1)+bp1*bp1
       bi1=ds1/bb1
       bm1=sqrt(bb1)
       br1=bm1*ds1
       bb1=bb1*ds1
       sh1=r1
       ai1=ds1/r1
       dl=sqrt((r1-r0)*(r1-r0)+(z1-z0)*(z1-z0))
       arv(n)=arv(n)+(r0+r1)*(z0-z1)*0.5d0
       vlv(n)=vlv(n)+(vl0+vl1)*(z0-z1)*0.5d0
       sdw(n)=sdw(n)+dl*(ds0+ds1)*0.5d0
       ckv(n)=ckv(n)+dl*(ck0+ck1)*0.5d0
       ssv(n)=ssv(n)+dl*(ss0+ss1)*0.5d0
       aav(n)=aav(n)+dl*(aa0+aa1)*0.5d0
       rrv(n)=rrv(n)+dl*(rr0+rr1)*0.5d0
       bbv(n)=bbv(n)+dl*(bb0+bb1)*0.5d0
       biv(n)=biv(n)+dl*(bi0+bi1)*0.5d0
       shv(n)=shv(n)+dl*(sh0+sh1)*0.5d0
       grbm2v(n)=grbm2v(n)+dl*(r0**2*bp0**2*bi0 &
            &                 +r1**2*bp1**2*bi1)*0.5d0
       aiv(n)=aiv(n)+dl*(ai0+ai1)*0.5d0
       brv(n)=brv(n)+dl*(br0+br1)*0.5d0

       dll(nsul(n)) = dl*(ds0+ds1)*0.5d0 ! dlp/Bp
       zbl(nsul(n)) = 0.5d0*(bm0+bm1)    ! B

       bmax(n) = max(bmax(n),zbl(nsul(n)))

!..geometric factors

       if(r1 > rrmax) rrmax=r1 ! rrmax = R_max
       if(r1 < rrmin) rrmin=r1 ! rrmin = R_min
       if(z1 > zzmax) then
          zzmax=z1 ! zzmax = Z_max
          rzmax=r1 ! rzmax = R at Z_max
       end if
       if(z1 < zzmin) then
          zzmin=z1 ! zzmin = Z_min
          rzmin=r1 ! rzmin = R at Z_min
       end if

       nsul(n)=nsul(n)+1
       if(n == nv) then
          nsu=nsu+1
          rsu(nsu)=r1
          zsu(nsu)=z1
          csu(nsu)=bp1*sigcu
       endif
!..
80     if( iudsym == 1 .and. iz > nz) go to 100
       if( i /= ist ) go to 20
!=======================================================================
100    if( iudsym == 1 ) then
          vlv(n)=2.d0*vlv(n)
          arv(n)=2.d0*arv(n)
          sdw(n)=2.d0*sdw(n)
          ckv(n)=2.d0*ckv(n)
          ssv(n)=2.d0*ssv(n)
          aav(n)=2.d0*aav(n)
          rrv(n)=2.d0*rrv(n)
          bbv(n)=2.d0*bbv(n)
          biv(n)=2.d0*biv(n)
          shv(n)=2.d0*shv(n)
          grbm2v(n)=2.d0*grbm2v(n)
          aiv(n)=2.d0*aiv(n)
          brv(n)=2.d0*brv(n)
          if( n == nv ) then
             lp=nsu
             lm=nsu
             do ll=2,nsu-1
                lp=lp+1
                lm=lm-1
                rsu(lp)= rsu(lm)
                zsu(lp)=-zsu(lm)
                csu(lp)= csu(lm)
             end do
             nsu=2*nsu-1
          endif
       endif
!-----
       vlv(n)=zpi*vlv(n)
       sdw(n)=1.d0/(2.d0*zpi*sdw(n))
       ckv(n)=2.d0*zpi*ckv(n)/sdw(n)
       r2b2v(n)=2.d0*zpi*ssv(n)*sdw(n)
       ssv(n)=2.d0*zpi*ssv(n)/sdw(n)
       aav(n)=2.d0*zpi*aav(n)*sdw(n)
       rrv(n)=2.d0*zpi*rrv(n)*sdw(n)
       bbv(n)=2.d0*zpi*bbv(n)*sdw(n)
       biv(n)=2.d0*zpi*biv(n)*sdw(n)
       shv(n)=2.d0*zpi*shv(n)*sdw(n)
       grbm2v(n)=2.d0*zpi*grbm2v(n)*sdw(n)
       aiv(n)=2.d0*zpi*aiv(n)*sdw(n)
       brv(n)=2.d0*zpi*brv(n)*sdw(n)
       sdw(n)=sdw(n)*sigcu

       rmajl=0.5d0*(rrmax+rrmin)
       rplal=0.5d0*(rrmax-rrmin)
       epsv(n) =(rrmax-rrmin)/(rrmax+rrmin)
       elipv(n)=(zzmax-zzmin)/(rrmax-rrmin)
       trigv(n)=(rmajl-0.5d0*(rzmax+rzmin))/rplal

       ! *** trapped particle fraction **********************************
       !
       !   ft = 1 - 0.75 <h^2> int_0^1 flam dflam / <sqrt(1 - flam * h)>
       !
       ! ****************************************************************

       nsul(n) = nsul(n) - 1
       do i = 1, nsul(n)
          h = zbl(i) / bmax(n) ! h
          do j = 0, intf
             if( 1.d0 - flam(j) * h < 0.d0 ) then
                write(6,'(a,1pe26.18,a)') &
                     &        '** WARNING  IN SQRT(DX),DX < 0.0(DX=',1.d0-flam(j)*h,').'
             endif
             fint(j) = fint(j) + sqrt( abs( 1.d0 - flam(j) * h ) ) * dll(i)
          enddo
       enddo

       fintx = 0.d0
       do j = 1, intf
          fintx = fintx + ( flam(j)**2 - flam(j-1)**2 ) / ( fint(j) + fint(j-1) )
       enddo
       hsq     = bbv(n) / bmax(n)**2 ! <h^2>
       fintx   = 0.75d0 * hsq * fintx / (2.d0 * zpi * sdw(n))
       ftv(n) = 1.d0 - fintx

    end do
!..evaluate on the axis
    siw(1)=saxis
    arv(1)=0.d0
    vlv(1)=0.d0
!    sdw(1)=-dpsi*(vlv(3)**2-2.d0*vlv(2)**2)/(vlv(3)*vlv(2)*(vlv(3)-vlv(2)))
    sdw2   = 0.5_8 * ( ( siw(2) - siw(1) ) / vlv(2) &
         &           + ( siw(3) - siw(2) ) / ( vlv(3) - vlv(2) ) )
    sdw(1) = ( vlv(3) * sdw2 - 2.0_8 * vlv(2) / vlv(3) &
         &           * ( siw(3) - siw(1) ) ) / ( vlv(3) - 2.0_8 * vlv(2) )
    ckv(1)=0.d0
    r2b2v(1)=0.d0
    ssv(1)=0.d0
    aav(1)=1.d0/(raxis*raxis)
    rrv(1)=raxis*raxis
    bbv(1)=(rbv(1)/raxis)**2
    biv(1)=1.d0/bbv(1)
    shv(1)=0.d0
    grbm2v(1)=0.d0
    aiv(1)=1.d0/raxis
    brv(1)=rbv(1)/raxis
    nsu = nsu - 1

    epsv(1)=0.d0
    elipv(1)=1.d0
    trigv(1)=0.d0

    ftv(1) =0.d0

    deallocate(bmax,fint,flam,nsul,dll,zbl)
    return
!*
800 stop 'eqlin : no starting point                '
810 stop 'eqlin : out of range                     '
  end subroutine eqlin

!=======================================================================
! Read G EQDSK files
!=======================================================================

  subroutine read_geqdsk(igeqdsk)
    use tx_commons, only : rMU0
    use equ_params
    integer(4), intent(inout) :: igeqdsk
    integer(4) :: i, ist, idum, nw, nh, nbbbs, limitr, kvtor, nmass
    real(8) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, &
         &     current, xdum, rvtor
    real(8), dimension(:), allocatable :: fpol, pres, ffprim, pprime, psirz, qpsi
    real(8), dimension(:), allocatable :: rbbbs, zbbbs, rlim, zlim
    real(8), dimension(:), allocatable :: pressw, pwprim, dmion, rhovn
    character(len=8) :: case(6)
    real(8) :: drg, dzg

    igeqdsk = 1

    !---- read EFIT-format
    read(ieqrd,2000,iostat=ist) (case(i),i=1,6), idum, nw, nh
    read(ieqrd,2020,iostat=ist)  rdim, zdim, rcentr, rleft, zmid
    read(ieqrd,2020,iostat=ist)  rmaxis, zmaxis, simag, sibry, bcentr
    read(ieqrd,2020,iostat=ist)  current, simag, xdum, rmaxis, xdum
    read(ieqrd,2020,iostat=ist)  zmaxis, xdum, sibry, xdum, xdum

    allocate(fpol(nw), pres(nw), ffprim(nw), pprime(nw), qpsi(nw))
    read(ieqrd,2020,iostat=ist) (fpol(i),i=1,nw) 
    read(ieqrd,2020,iostat=ist) (pres(i),i=1,nw)
    read(ieqrd,2020,iostat=ist) (ffprim(i),i=1,nw)
    read(ieqrd,2020,iostat=ist) (pprime(i),i=1,nw)
    allocate(psirz(nw*nh))
    read(ieqrd,2020,iostat=ist) (psirz(i),i=1,nw*nh)
    read(ieqrd,2020,iostat=ist) (qpsi(i),i=1,nw)

    read(ieqrd,2022,iostat=ist)  nbbbs,limitr
    allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))
    read(ieqrd,2020,iostat=ist) (rbbbs(i),zbbbs(i),i=1,nbbbs)
    read(ieqrd,2020,iostat=ist) (rlim(i),zlim(i),i=1,limitr)

    if( ist /= 0 ) stop 'G-EQDSK file read error! Processing end ...'

    !---- additional data
    allocate(pressw(nw),pwprim(nw),dmion(nw),rhovn(nw))
    read(ieqrd,2024,iostat=ist)  kvtor,rvtor,nmass
    if ( kvtor > 0 ) then
       read(ieqrd,2020,iostat=ist) (pressw(i),i=1,nw)
       read(ieqrd,2020,iostat=ist) (pwprim(i),i=1,nw)
    endif
    if ( nmass > 0 ) then
       read(ieqrd,2020,iostat=ist) (dmion(i),i=1,nw)
    endif
    if( ist < 0 ) stop 'G-EQDSK file read error! Processing end ...'
    read(ieqrd,2020,iostat=ist) (rhovn(i),i=1,nw)

    nsr = nw
    nsz = nh

    !---- make rg
    rg(1)  = rleft
    rg(nw) = rdim + rleft
    drg    = rdim / real(nw-1,8)
    do i = 2,nw-1
       rg(i) = rg(1) + drg*(i-1)
    enddo

    !---- make zg
    zg(1)  = zmid - 0.5d0 * zdim
    zg(nh) = zmid + 0.5d0 * zdim
    dzg    = zdim / real(nh-1,8)
    do i = 2,nh-1
       zg(i) = zg(1) + dzg*(i-1)
    enddo

    psi(1:nsr*nsz) = psirz(1:nsr*nsz) - sibry

    nv  = nw

    pds(1:nv) = pprime(1:nv) * rMU0
    fds(1:nv) = ffprim(1:nv)
    ! vlv will be computed in eqlin.
    qqv(1:nv) = qpsi(1:nv)
    prv(1:nv) = pres(1:nv) * rMU0

    btv   = rcentr * bcentr
    ttcu  = current * 1.d-6
    ! no values associated with ttpr, bets, beta, betj
    ttpr  = 0.d0
    bets  = 0.d0
    beta  = 0.d0
    betj  = 0.d0
    saxis = simag - sibry
    raxis = rmaxis
    zaxis = zmaxis

    nsu        = nbbbs
    rsu(1:nsu) = rbbbs(1:nsu)
    zsu(1:nsu) = zbbbs(1:nsu)
    ! no arrays associated with csu
    csu(:) = 0.d0

    ilimt          = limitr
    rlimt(1:ilimt) = rlim(1:ilimt)
    zlimt(1:ilimt) = zlim(1:ilimt)


    deallocate(fpol, pres, ffprim, pprime, psirz, qpsi)
    deallocate(rbbbs, zbbbs, rlim, zlim)
    deallocate(pressw,pwprim,dmion,rhovn)

2000 format(6a8,3i4)
2020 format(5e16.9)
2022 format(2i5)
2024 format(i5,e16.9,i5)

  end subroutine read_geqdsk


!!$!***************************************************************
!!$!
!!$!   Wrapper for lesq6
!!$!
!!$!***************************************************************
!!$
!!$  subroutine lesq6_wrapper(ilesq6, mdeg, x, y, xout, yout)
!!$! arguments
!!$    integer(4), intent(in) :: ilesq6, mdeg
!!$    real(8), dimension(:), intent(in)  :: x, y, xout
!!$    real(8), dimension(:), intent(out) :: yout
!!$! local variables
!!$    integer(4) :: nmax, icon, i, j, nout
!!$    real(8), dimension(:), allocatable :: ww, wc
!!$
!!$    nmax = size(x)
!!$    nout = size(yout)
!!$
!!$    allocate(ww(nmax), wc(nmax))
!!$    ww(:) = 1.d0
!!$    call lesq6( ilesq6, x, y, nmax, mdeg, ww, wc, icon )
!!$    if( icon /= 0 ) then
!!$       write(6,'(a)') '*** lesq6 error ***'
!!$       stop
!!$    endif
!!$    yout(1) = wc(1)
!!$    do i = 2, nout
!!$       yout(i) = 0.0_8
!!$       do j = 1, mdeg+1
!!$          yout(i) = yout(i) + wc(j) * xout(i)**(j-1)
!!$       enddo
!!$    enddo
!!$    deallocate(ww, wc)
!!$
!!$  end subroutine lesq6_wrapper
!!$
!!$!=======================================================================
!!$  subroutine lesq6( lsm, x, y, n, m, w, c, icon )
!!$!                   I    I  I  I  I  I  O  O
!!$!
!!$!     Function :
!!$!        Calculate the Coefficient of Least Square Polynomial Function
!!$!      Approximation
!!$!
!!$!     Arguement
!!$!        LSM  : Fitting Mode
!!$!             :   IAND(LSM,1).NE.0 : Derivative at Center = 0
!!$!             :   IAND(LSM,2).NE.0 : Fix Center Value ( Y(1) )
!!$!             :   IAND(LSM,4).NE.0 : Fix Edge Value   ( Y(N) )
!!$!        X    : X Data
!!$!        Y    : Y Data
!!$!        N    : Number of Data Points
!!$!        M    : Polynomial Function Degree
!!$!        W    : Weight of Each Data
!!$!        C    : Coefficient of Polynomial Function
!!$!        ICON : Return Code
!!$!=======================================================================
!!$! arguments
!!$    integer, intent(in)  :: lsm, n, m
!!$    real(8), intent(in)  :: x(n), y(n), w(n)
!!$    integer, intent(out) :: icon
!!$    real(8), intent(out) :: c(m+1)
!!$
!!$! local variables
!!$    integer    i, irc, j, k, m1
!!$    real(8)    a1, b1, xa(m+1,m+2), xc(n,m+1)
!!$!        XA   : Work Area Used for Coefficient Matrix of Normal Equation
!!$!        XC   : Work Area Used for Saving the X**(j)
!!$
!!$!-----------------------------------------------------------------------
!!$
!!$!----  Initialize and Check Data
!!$
!!$    icon = 1
!!$    if( n .lt. m+1 ) return
!!$    m1 = m + 1
!!$
!!$!----  Set Xi**(j),j=0,M
!!$
!!$    do i=1,n
!!$       xc(i,1) = 1.d0
!!$    enddo
!!$    do j=2,m1
!!$       do i=1,n
!!$          xc(i,j) = xc(i,j-1)*x(i)
!!$       enddo
!!$    enddo
!!$
!!$!----  Calculate the Coefficient Matrix of Normal Equation
!!$
!!$!----           n
!!$!----    XAij = Σ(Xk**(i) * Xk**(j) * Wk)
!!$!----           k=1
!!$
!!$    do i=1,m1
!!$       do j=i,m1
!!$          a1 = 0.0d0
!!$          do k=1,n
!!$             a1 = a1 + xc(k,i)*xc(k,j)*w(k)
!!$          enddo
!!$          xa(i,j) = a1
!!$          xa(j,i) = a1
!!$       enddo
!!$    enddo
!!$
!!$!----  Calculate the Constant Vector of Normal Equation
!!$
!!$!----         n
!!$!----    Bj = Σ(Xk**(j) * Yk * Wk)
!!$!----         k=1
!!$
!!$    do i=1,m1
!!$       b1 = 0.0d0
!!$       do k=1,n
!!$          b1 = b1 + xc(k,i)*y(k)*w(k)
!!$       enddo
!!$       xa(i,m+2) = b1
!!$    enddo
!!$!    write(ft06,*)'XA= '
!!$!    write(ft06,'(1X,1P5E11.3)')((XA(I,J),J=1,5),I=1,4)
!!$
!!$!----  Derivative at Center = 0
!!$
!!$    if( iand(lsm,1).ne.0 ) then
!!$       xa(2,1) = 0.0d0
!!$       xa(2,2) = 1.0d0
!!$       do i=3,m+1
!!$          xa(2,i) = 0.0d0
!!$       enddo
!!$       xa(2,m+2) = 0.0d0
!!$    endif
!!$
!!$!----  Fix Center Value
!!$
!!$    if( iand(lsm,2).ne.0 ) then
!!$       xa(1,1) = 1.0d0
!!$       do i=2,m+1
!!$          xa(1,i) = 0.0d0
!!$       enddo
!!$       xa(1,m+2) = y(1)
!!$    endif
!!$
!!$!----  Fix Edge Value
!!$
!!$    if( iand(lsm,4).ne.0 ) then
!!$       do i=1,m+1
!!$          xa(m+1,i) = 1.0d0
!!$       enddo
!!$       xa(m+1,m+2) = y(n)
!!$    endif
!!$
!!$!----  Solve the Equation by a Process of Gauss-Jordan Elimination
!!$
!!$    call gauss6( xa, m1, irc )
!!$    if( irc.ne.0 ) return
!!$
!!$!----  Set Coefficient
!!$
!!$    do i=1,m1
!!$       c(i) = xa(i,m1+1)
!!$    enddo
!!$    icon = 0
!!$!-----------------------------------------------------------------------
!!$  end subroutine lesq6
!!$
!!$
!!$!=======================================================================
!!$  subroutine gauss6( a, m, ic )
!!$!                    U  I  O
!!$!
!!$!     Function :
!!$!        Solve the Equation by a Process of Gauss-Jordan Elimination
!!$!
!!$!         |A11 A12 ... A1m|   | X1 |   |A1,m+1|
!!$!         |A21 A22 ... A2m| * | X2 | = |A2,m+1|
!!$!         | |   |   |   | |   | |  |   |  |   |
!!$!         |Am1 Am2 ... Amm|   | Xm |   |Am,m+1|
!!$!
!!$!     Arguement
!!$!        A    : Coefficient Matrix and Constant Vector
!!$!             : A(*,1:M) : Coefficient Matrix
!!$!             : A(*,M+1) : Inp : Constant Vector
!!$!             :          : Out : Solution to Simultaneous Equation
!!$!        M    : Number of Simultaneous Equation
!!$!        IC   : Return Code
!!$!=======================================================================
!!$! arguments
!!$    integer, intent(in)    :: m
!!$    integer, intent(out)   :: ic
!!$    real(8), intent(inout) :: a(m,m+1)
!!$
!!$! local variables
!!$    integer    i, ip, j, k, kk
!!$    real(8)    aa, epsl, pv, w
!!$
!!$!-----------------------------------------------------------------------
!!$
!!$!----  Initialize Arguement and Data
!!$
!!$    ic = 1
!!$! modified 1/1 lines by kamata 2003/10/16
!!$!   epsl = 1.0d-14
!!$    epsl = 1.0d-28
!!$
!!$!----  Loop of the Number of Row
!!$
!!$    do k=1,m
!!$
!!$!------  Select Pivot
!!$
!!$       pv = 0.0d0
!!$       do i=k,m
!!$          if( pv.lt.dabs(a(i,k)) ) then
!!$             pv = dabs(a(i,k))
!!$             ip = i
!!$          endif
!!$       enddo
!!$       if( pv.lt.epsl ) then
!!$          write(6,*) 'ERROR AT GAUSS : PIVOT'
!!$          write(6,'(1P6E12.5)') (a(i,k),i=k,m)
!!$          return
!!$       endif
!!$
!!$!------  Exchange the Rows
!!$
!!$       do j=k,m+1
!!$          w       = a(k,j)
!!$          a(k,j)  = a(ip,j)
!!$          a(ip,j) = w
!!$       enddo
!!$
!!$!------  Process for Row Selected as Pivot
!!$
!!$       aa = 1.0d0 / a(k,k)
!!$       kk = k + 1
!!$       do j=kk,m+1
!!$          a(k,j) = a(k,j) * aa
!!$       enddo
!!$
!!$!------  Process for Other Rows
!!$
!!$       do i=1,m
!!$          if(i.eq.k) cycle !  pivot
!!$          aa = a(i,k)
!!$          do j=kk,m+1
!!$             a(i,j) = a(i,j) - aa*a(k,j)
!!$          enddo
!!$       enddo
!!$
!!$    enddo
!!$
!!$    ic = 0
!!$!-----------------------------------------------------------------------
!!$  end subroutine gauss6

end module eqread_mod

