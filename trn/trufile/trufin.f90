MODULE trufin
! ------------------------------------------------------------------------
!  ***  time interpolating to the time of every step
!  ***  subsitute experimetal data into calculation variables
!  ***  store values for graphic output
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,nsum,nrum,ntum, nrmax,ntxmax, tmu
  USE trlib, ONLY: lin_itp

  IMPLICIT NONE

  PUBLIC 

  ! intermediate variables for graphic output ----------------------------
  REAL(rkind) :: &
       ! golobal variables
       tmug,rrug,raug,phiaug,bbug,rkapug,rdltug,ripug,wthug,wtotug,zeffug, &
       ! source
       pnbug,pecug,pibwug,picug,plhug,pohmug,pradug

  REAL(rkind),DIMENSION(1:nrum) :: & 
       ! profile variables
       qpug,bpug,zeffrug,wrotug, &
       ! current density
       jtotug,jnbug,jecug,jicug,jlhug,jbsug 

  REAL(rkind),DIMENSION(2,1:nrum) :: & ! electron and ion
       ! power deposition profile
       qnbug,qecug,qibwug,qicug,qlhug,qfusug,qwallug, &
       ! particle source
       snbug

  REAL(rkind),DIMENSION(1:nrum) :: & ! profile variables
       ! power deposition profile
       qohmug,qradug, &
       ! particle source
       swallug,       &
       ! geometric quantity
       pvolug,psurug,rmjrhoug,rmnrhoug,ar1rhoug,ar2rhoug,rkprhoug, &
       dvrhoug,arrhoug,abrhoug,ttrhoug

  REAL(rkind),DIMENSION(1:nsum,1:nrum) :: &
       rnug,rnfug,rtug
  ! ----------------------------------------------------------------------

  REAL(rkind),DIMENSION(1:ntum) :: temp

CONTAINS

  SUBROUTINE tr_uf_init(mdluf)
    USE trcomm,ONLY: ntmax,dt,tlmax,ntlmax,tmu,time_slc

    INTEGER(ikind),INTENT(IN) :: mdluf

    INTEGER(ikind) :: nt,ntxinit
    REAL(rkind)    :: tmax

    SELECT CASE(mdluf)
    CASE(1)

       ! get the initial time from the exp. data at an arbitrary moment
       CALL tr_uf_time_slice(time_slc,tmu,ntxmax,ntxinit)

    CASE(2)
       nt = 0
       ntmax_set: DO
          nt   = nt + 1
          tmax = nt * dt
          IF(tmax > tlmax)THEN
             ntlmax = nt - 1
             EXIT
          END IF
       END DO ntmax_set
!       ntmax = ntlmax
       WRITE(6,*) ' - Set NTLMAX to fit with the exp.data. NTLMAX= ',ntlmax
       WRITE(6,*) ! spacing

    END SELECT
    
    RETURN
  END SUBROUTINE tr_uf_init

! **********************************************************************

! ----------------------------------------------------------------------
! Substitution routines of experimental data into calculation variables
!
!   - tr_ufin_global      : global variables (RR,RA,BB, ...etc.)
!   - tr_ufin_density     : density profiles including fast ions
!   - tr_ufin_rotation    : toroidal velocity
!   - tr_ufin_temperature : temperature profiles
!   - tr_ufin_field       : quantities associated with poloidal flux
!   - tr_ufin_source      : source profile (density and energy)
!   - tr_ufin_geometry    : geometric quantities
!
! < input >
!   time : the time at which the experimental data are substituted [sec]
!   tmid = 0 : the exp. data are not substituted
!        = 1 : all the exp. data are substituted
!        = 2 : the exp. data are substituted depending on 'id_neq'
!  mdlid : (only for 'tr_ufin_field')
!          the switch corresponding to 'mdlijq' in trcomm.f90
!        = 1 : create from jtot using RIP as boundary condition
!        = 2 : create from qp using RIP as boundary condition  
!        = 3 : create from jtot not using RIP                  
!        = 4 : create from qp not using RIP                    
!
! < output >
!   ierr : error identifier
! ----------------------------------------------------------------------

  SUBROUTINE tr_ufin_global(time,tmid,ierr)
! ----------------------------------------------------------------------
!   global variables
! ----------------------------------------------------------------------
    USE trcomm,ONLY: &
             rr, ra, rip, bb, rkap, rdlt, &
         tmu,rru,rau,ripu,bbu,rkapu,rdltu,phiau,zeffu,wthu,wtotu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr

    ierr = 0
!!$    CALL TIMESPL(time,rrug,tmu,rru,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,raug,tmu,rau,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,bbug,tmu,bbu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,rkapug,tmu,rkapu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,rdltug,tmu,rdltu,ntxmax,ntum,ierr)
!!$!    CALL TIMESPL(time,phiaug,tmu,phiau,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,zeffug,tmu,zeffu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,wthug,tmu,wthu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,wtotug,tmu,wtotu,ntxmax,ntum,ierr)

    CALL lin_itp(time,rrug,tmu,rru,ntxmax,ntum)
    CALL lin_itp(time,raug,tmu,rau,ntxmax,ntum)
    CALL lin_itp(time,bbug,tmu,bbu,ntxmax,ntum)
    CALL lin_itp(time,rkapug,tmu,rkapu,ntxmax,ntum)
    CALL lin_itp(time,rdltug,tmu,rdltu,ntxmax,ntum)
!    CALL lin_itp(time,phiaug,tmu,phiau,ntxmax,ntum)
    CALL lin_itp(time,zeffug,tmu,zeffu,ntxmax,ntum)
    CALL lin_itp(time,wthug,tmu,wthu,ntxmax,ntum)
    CALL lin_itp(time,wtotug,tmu,wtotu,ntxmax,ntum)

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_global: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    ! substitution for variables used in calculation
    IF(tmid > 0)THEN
       RR   = rrug
       RA   = raug
       BB   = ABS(bbug)
       RKAP = rkapug
       RDLT = rdltug
    END IF

    RETURN
  END SUBROUTINE tr_ufin_global


  SUBROUTINE tr_ufin_density(time,tmid,ierr)
! ----------------------------------------------------------------------
!   density and ZEFF
! ----------------------------------------------------------------------
    USE trcomm,ONLY: id_neq,id_neqnr,nsa_nsu,nsa_nsfu, &
                     tmu,rn,rnu,rnfu,z_eff,zeffru
    USE truf0d,ONLY: idnm,idnfast
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind) :: nr,nsa,neq,nsu,nsfu,id

    ierr = 0

    DO nr = 1, nrmax+1
       DO nsu = 1, nsum
          IF(idnm(nsu))THEN
             temp(1:ntum) = rnu(nsu,1:ntum,nr)
!             CALL TIMESPL(time,rnug(nsu,nr),tmu,temp,ntxmax,ntum,ierr)
             CALL lin_itp(time,rnug(nsu,nr),tmu,temp,ntxmax,ntum)
          END IF
          IF(idnfast(nsu))THEN
             temp(1:ntum) = rnu(nsu,1:ntum,nr)
!             CALL TIMESPL(time,rnfug(nsu,nr),tmu,temp,ntxmax,ntum,ierr)
             CALL lin_itp(time,rnfug(nsu,nr),tmu,temp,ntxmax,ntum)
          END IF
       END DO

!       CALL TIMESPL(time,zeffrug(nr),tmu,zeffru(:,nr),ntxmax,ntum,ierr)
       temp(1:ntum) = zeffru(1:ntum,nr)
       CALL lin_itp(time,zeffrug(nr),tmu,zeffru(:,nr),ntxmax,ntum)
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_density: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF


    ! substitution for variables used in calculation
    IF(tmid > 0)THEN
       DO nsu = 1, nsum
          IF(idnm(nsu))THEN
             nsa = nsa_nsu(nsu)
             IF(nsa == 0) CYCLE

             neq = 3*nsa - 1
             id  = id_neq(neq)
             IF(tmid == 1) id = 0 ! initial profile

             SELECT CASE(id)
             CASE(0) ! input values over all radial points
                rn(nsa,0:nrmax) = rnug(nsu,1:nrmax+1)
             CASE(2) ! input value only at plasma surface
                rn(nsa,nrmax) = rnug(nsu,nrmax+1)
             END SELECT
          END IF
       END DO
       DO nsfu = 1, nsum
          IF(idnfast(nsfu))THEN
             nsa = nsa_nsfu(nsfu)
             IF(nsa == 0) CYCLE

             neq = 3*nsa - 1
             id  = id_neq(neq)
             IF(tmid == 1) id = 0 ! initial profile

             SELECT CASE(id)
             CASE(0)   ! input values over all radial points
                rn(nsa,0:nrmax) = rnfug(nsfu,1:nrmax+1)
             CASE(1) ! input value only at plasma surface
                rn(nsa,nrmax) = rnfug(nsfu,nrmax+1)
             CASE(11) ! input values over the fixed region
                DO nr = 0, nrmax
                   IF(id_neqnr(neq,nr)==2) rn(nsa,nr) = rnfug(nsfu,nr+1)
                END DO
             END SELECT
          END IF
       END DO
    END IF

!    z_eff(0:nrmax) = zeffrug(1:nrmax+1)
    RETURN
  END SUBROUTINE tr_ufin_density


  SUBROUTINE tr_ufin_rotation(time,tmid,ierr)
! ----------------------------------------------------------------------
!   toroidal rotation
!   * interim way
! ----------------------------------------------------------------------
    USE trcomm, ONLY: tmu,wrot,wrotu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr

    INTEGER(ikind) :: nr

    ierr = 0

    DO nr = 1, nrmax+1
       temp(1:ntum) = wrotu(1:ntum,nr)
!       CALL TIMESPL(time,wrotug(nr),tmu,temp,ntxmax,ntum,ierr)
       CALL lin_itp(time,wrotug(nr),tmu,temp,ntxmax,ntum)
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_rotation: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    IF(tmid > 0)THEN
       wrot(0:nrmax) = wrotug(1:nrmax+1)
    END IF

    RETURN
  END SUBROUTINE tr_ufin_rotation


  SUBROUTINE tr_ufin_temperature(time,tmid,ierr)
! -----------------------------------------------------------------------
!   temperature
! -----------------------------------------------------------------------
    USE trcomm, ONLY: id_neq,id_neqnr,nsa_nsu,rt,rtu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr

    INTEGER(ikind) :: nr,nsu,nsa,neq,id

    ierr = 0

    DO nr = 1, nrmax+1
       DO nsu = 1, nsum
          temp(1:ntum) = rtu(nsu,1:ntum,nr)
!          CALL TIMESPL(time,rtug(nsu,nr),tmu,temp,ntxmax,ntum,ierr)
          CALL lin_itp(time,rtug(nsu,nr),tmu,temp,ntxmax,ntum)
       END DO
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_temperature: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    IF(tmid > 0)THEN
       DO nsu = 1, nsum
          nsa = nsa_nsu(nsu)
          IF(nsa == 0) CYCLE
          
          neq = 3*nsa + 1
          id  = id_neq(neq)
          IF(tmid == 1) id = 0 ! initial profile

          SELECT CASE(id)
          CASE(0) ! input values over all radial points
             rt(nsa,0:nrmax) = rtug(nsu,1:nrmax+1)
          CASE(1) ! input value only at plasma surface
             rt(nsa,nrmax) = rtug(nsu,nrmax+1)
          CASE(11) ! input values over the fixed region
             DO nr = 0, nrmax
                IF(id_neqnr(neq,nr)==2) rt(nsa,nr) = rtug(nsu,nr+1)
             END DO
          END SELECT
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_ufin_temperature


  SUBROUTINE tr_ufin_field(time,tmid,mdlid,ierr)
! ----------------------------------------------------------------------
!   profile variables associated with poloidal magnetic field
! ----------------------------------------------------------------------
    USE trcomm,ONLY: rkev,id_neq,                                  &
             rip, qp, bp, jtot,jcd_nb,jcd_ec,jcd_ic,jcd_lh,jbs_nc, &
         tmu,ripu,qpu,bpu,jtotu,jnbu,jecu,jicu,jlhu,jbsu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid,mdlid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind) :: nr,id

    ierr = 0

!!$    CALL TIMESPL(time,ripug,tmu,ripu,ntxmax,ntum,ierr)
    CALL lin_itp(time,ripug,tmu,ripu,ntxmax,ntum)

    DO nr = 1, nrmax+1
!!$       CALL TIMESPL(time,qpug(nr),tmu,qpu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,bpug(nr),tmu,bpu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,jtotug(nr),tmu,jtotu(:,nr),ntxmax,ntum,ierr)
!!$
!!$       CALL TIMESPL(time,jnbug(nr),tmu,jnbu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,jecug(nr),tmu,jecu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,jicug(nr),tmu,jicu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,jlhug(nr),tmu,jlhu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,jbsug(nr),tmu,jbsu(:,nr),ntxmax,ntum,ierr)

       temp(1:ntum) = qpu(1:ntum,nr)
       CALL lin_itp(time,qpug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = bpu(1:ntum,nr)
       CALL lin_itp(time,bpug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = jtotu(1:ntum,nr)
       CALL lin_itp(time,jtotug(nr),tmu,temp,ntxmax,ntum)

       temp(1:ntum) = jnbu(1:ntum,nr)
       CALL lin_itp(time,jnbug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = jecu(1:ntum,nr)
       CALL lin_itp(time,jecug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = jicu(1:ntum,nr)
       CALL lin_itp(time,jicug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = jlhu(1:ntum,nr)
       CALL lin_itp(time,jlhug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = jbsu(1:ntum,nr)
       CALL lin_itp(time,jbsug(nr),tmu,temp,ntxmax,ntum)
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_profile: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    ! substitution for variables used in calculation
    IF(tmid > 0)THEN
       id = id_neq(1)
       IF(tmid == 1) id = 0 ! initial profile

       IF(id == 0)THEN ! magnetic diffusion equation is not solved.
          SELECT CASE(mdlid)
          CASE(1)
             rip  = - ripug * 1.d-6
             jtot(0:nrmax) = - jtotug(1:nrmax+1)
          CASE(2)
             rip  = - ripug * 1.d-6
             qp(0:nrmax)   = qpug(1:nrmax+1)
          CASE(3)
             jtot(0:nrmax) = - jtotug(1:nrmax+1)
          CASE(4)
             qp(0:nrmax)   = qpug(1:nrmax+1)
          END SELECT

       ELSE IF(id == 2)THEN ! magnetic diffusion equation is solved.
          SELECT CASE(mdlid)
          CASE(1,2)
             rip  = - ripug * 1.d-6
          CASE(3)
             jtot(0:nrmax) = - jtotug(1:nrmax+1)
          CASE(4)
             qp(nrmax)   = qpug(nrmax+1)
          END SELECT
       END IF

       jcd_nb(0:nrmax) = - jnbug(1:nrmax+1)
       jcd_ec(0:nrmax) = - jecug(1:nrmax+1)
       jcd_ic(0:nrmax) = - jicug(1:nrmax+1)
       jcd_lh(0:nrmax) = - jlhug(1:nrmax+1)
       jbs_nc(0:nrmax) = - jbsug(1:nrmax+1)
    END IF

    RETURN
  END SUBROUTINE tr_ufin_field


  SUBROUTINE tr_ufin_source(time,tmid,ierr)
! ----------------------------------------------------------------------
!   power input and source/sink terms
! ----------------------------------------------------------------------
    USE trcomm,ONLY: tmu,pnbu,pecu,pibwu,picu,plhu,pohmu,pradu,    &
         jtotu,jnbu,jecu,jicu,jlhu,jbsu,qnbu,qecu,qibwu,qicu,qlhu, &
         qfusu,qohmu,qradu,qwallu,snbu,swallu, &
         pnb,pec,pibw,pic,plh,poh,prl,pwl,pnf,snb,spl,swl,  &
         pnb_t,pec_t,pibw_t,pic_t,plh_t,poh_t,prl_t

    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: tmid
    REAL(rkind)   ,INTENT(IN)  :: time
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind) :: nr

    ierr = 0

!!$    CALL TIMESPL(time,pnbug,tmu,pnbu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,pecug,tmu,pecu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,pibwug,tmu,pibwu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,picug,tmu,picu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,plhug,tmu,plhu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,pohmug,tmu,pohmu,ntxmax,ntum,ierr)
!!$    CALL TIMESPL(time,pradug,tmu,pradu,ntxmax,ntum,ierr)

    CALL lin_itp(time,pnbug,tmu,pnbu,ntxmax,ntum)
    CALL lin_itp(time,pecug,tmu,pecu,ntxmax,ntum)
    CALL lin_itp(time,pibwug,tmu,pibwu,ntxmax,ntum)
    CALL lin_itp(time,picug,tmu,picu,ntxmax,ntum)
    CALL lin_itp(time,plhug,tmu,plhu,ntxmax,ntum)
    CALL lin_itp(time,pohmug,tmu,pohmu,ntxmax,ntum)
    CALL lin_itp(time,pradug,tmu,pradu,ntxmax,ntum)

    DO nr = 1, nrmax+1
!!$       CALL TIMESPL(time,qnbug(1,nr),tmu,qnbu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qnbug(2,nr),tmu,qnbu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qecug(1,nr),tmu,qecu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qecug(2,nr),tmu,qecu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qibwug(1,nr),tmu,qibwu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qibwug(2,nr),tmu,qibwu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qicug(1,nr),tmu,qicu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qicug(2,nr),tmu,qicu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qlhug(1,nr),tmu,qlhu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qlhug(2,nr),tmu,qlhu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qohmug(nr),tmu,qohmu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qradug(nr),tmu,qradu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qfusug(1,nr),tmu,qfusu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,qfusug(2,nr),tmu,qfusu(2,:,nr),ntxmax,ntum,ierr)
!!$
!!$       CALL TIMESPL(time,snbug(1,nr),tmu,snbu(1,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,snbug(2,nr),tmu,snbu(2,:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,swallug(nr),tmu,swallu(:,nr),ntxmax,ntum,ierr)

       temp(1:ntum) = qnbu(1,1:ntum,nr)
       CALL lin_itp(time,qnbug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qnbu(2,1:ntum,nr)
       CALL lin_itp(time,qnbug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qecu(1,1:ntum,nr)
       CALL lin_itp(time,qecug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qecu(2,1:ntum,nr)
       CALL lin_itp(time,qecug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qibwu(1,1:ntum,nr)
       CALL lin_itp(time,qibwug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qibwu(2,1:ntum,nr)
       CALL lin_itp(time,qibwug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qicu(1,1:ntum,nr)
       CALL lin_itp(time,qicug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qicu(2,1:ntum,nr)
       CALL lin_itp(time,qicug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qlhu(1,1:ntum,nr)
       CALL lin_itp(time,qlhug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qlhu(2,1:ntum,nr)
       CALL lin_itp(time,qlhug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qohmu(1:ntum,nr)
       CALL lin_itp(time,qohmug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qradu(1:ntum,nr)
       CALL lin_itp(time,qradug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qwallu(1,1:ntum,nr)
       CALL lin_itp(time,qwallug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qwallu(2,1:ntum,nr)
       CALL lin_itp(time,qwallug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qfusu(1,1:ntum,nr)
       CALL lin_itp(time,qfusug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = qfusu(2,1:ntum,nr)
       CALL lin_itp(time,qfusug(2,nr),tmu,temp,ntxmax,ntum)

       temp(1:ntum) = snbu(1,1:ntum,nr)
       CALL lin_itp(time,snbug(1,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = snbu(2,1:ntum,nr)
       CALL lin_itp(time,snbug(2,nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = swallu(1:ntum,nr)
       CALL lin_itp(time,swallug(nr),tmu,temp,ntxmax,ntum)
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_source: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    ! substitution for variables used in calculation
    IF(tmid > 0)THEN
       ! global variables
       pnb_t  = pnbug *1.d-6
       pec_t  = pecug *1.d-6
       pibw_t = pibwug *1.d-6
       pic_t  = picug *1.d-6
       plh_t  = plhug *1.d-6
       poh_t  = pohmug *1.d-6
       prl_t  = pradug *1.d-6

       ! profiles
       pnb(1,0:nrmax)  = qnbug(1,1:nrmax+1)
       pnb(2,0:nrmax)  = qnbug(2,1:nrmax+1)
       pec(1,0:nrmax)  = qecug(1,1:nrmax+1)
       pec(2,0:nrmax)  = qecug(2,1:nrmax+1)
       pibw(1,0:nrmax) = qibwug(1,1:nrmax+1)
       pibw(2,0:nrmax) = qibwug(2,1:nrmax+1)
       pic(1,0:nrmax)  = qicug(1,1:nrmax+1)
       pic(2,0:nrmax)  = qicug(2,1:nrmax+1)
       plh(1,0:nrmax)  = qlhug(1,1:nrmax+1)
       plh(2,0:nrmax)  = qlhug(2,1:nrmax+1)
       poh(1,0:nrmax)  = qohmug(1:nrmax+1)
       prl(1,0:nrmax)  = qradug(1:nrmax+1)
       pwl(1,0:nrmax)  = qwallug(1,1:nrmax+1)
       pwl(2,0:nrmax)  = qwallug(2,1:nrmax+1)
       pnf(1,0:nrmax)  = qfusug(1,1:nrmax+1)
       pnf(2,0:nrmax)  = qfusug(2,1:nrmax+1)

       snb(1,0:nrmax)  = snbug(1,1:nrmax+1)
       snb(2,0:nrmax)  = snbug(2,1:nrmax+1)
       ! for the time being
       swl(1,0:nrmax)  = swallug(1:nrmax+1)
    END IF

    RETURN
  END SUBROUTINE tr_ufin_source


  SUBROUTINE tr_ufin_geometry(time,tmid,ierr)
! ----------------------------------------------------------------------
!   metric quantity
! ----------------------------------------------------------------------
    USE trcomm,ONLY: tmu,pvolu,psuru,rmjrhou,rmnrhou,ar1rhou,ar2rhou, &
                     rkprhou,dvrhou,arrhou,abrhou,ttrhou,             &
                     pvolrho,psurrho,rmjrho,rmnrho,ar1rho,ar2rho,     &
                     rkprho,dvrho
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN) :: tmid
    REAL(rkind)   ,INTENT(IN) :: time
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind) :: nr

    ierr = 0

    DO nr = 1, nrmax+1
!!$       CALL TIMESPL(time,pvolug(nr),tmu,pvolu(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,psurug(nr),tmu,psuru(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,rmjrhoug(nr),tmu,rmjrhou(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,rmnrhoug(nr),tmu,rmnrhou(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,ar1rhoug(nr),tmu,ar1rhou(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,ar2rhoug(nr),tmu,ar2rhou(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,rkprhoug(nr),tmu,rkprhou(:,nr),ntxmax,ntum,ierr)
!!$       CALL TIMESPL(time,dvrhoug(nr),tmu,dvrhou(:,nr),ntxmax,ntum,ierr)

       temp(1:ntum) = pvolu(1:ntum,nr)
       CALL lin_itp(time,pvolug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = psuru(1:ntum,nr)
       CALL lin_itp(time,psurug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = rmjrhou(1:ntum,nr)
       CALL lin_itp(time,rmjrhoug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = rmnrhou(1:ntum,nr)
       CALL lin_itp(time,rmnrhoug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = ar1rhou(1:ntum,nr)
       CALL lin_itp(time,ar1rhoug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = ar2rhou(1:ntum,nr)
       CALL lin_itp(time,ar2rhoug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = rkprhou(1:ntum,nr)
       CALL lin_itp(time,rkprhoug(nr),tmu,temp,ntxmax,ntum)
       temp(1:ntum) = dvrhou(1:ntum,nr)
       CALL lin_itp(time,dvrhoug(nr),tmu,temp,ntxmax,ntum)
    END DO

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX uf_ufin_metric: Time interpolation error. IERR=', ierr
       WRITE(6,*) 'XX skipped to interpolate time data.'
       RETURN
    END IF

    ! substitution for variables used in calculation
    IF(tmid > 0)THEN
       pvolrho(0:nrmax) = pvolug(1:nrmax+1)
       psurrho(0:nrmax) = psurug(1:nrmax+1)
       rmjrho(0:nrmax)  = rmjrhoug(1:nrmax+1)
       rmnrho(0:nrmax)  = rmnrhoug(1:nrmax+1)
       ar1rho(0:nrmax)  = ar1rhoug(1:nrmax+1)
       ar2rho(0:nrmax)  = ar2rhoug(1:nrmax+1)
       rkprho(0:nrmax)  = rkprhoug(1:nrmax+1)
       dvrho(0:nrmax)   = dvrhoug(1:nrmax+1)
    END IF

    RETURN
  END SUBROUTINE tr_ufin_geometry

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_uf_time_slice(time_slc,tl,ntxumax,ntsl)
! ------------------------------------------------------------------------
!   acquire the profile of a 2D variable at the certain time point
! ------------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN)  :: ntxumax
    INTEGER(ikind),INTENT(OUT) :: ntsl
    REAL(rkind),                INTENT(INOUT) :: time_slc
    REAL(rkind),DIMENSION(ntum),INTENT(IN)    :: tl

    INTEGER(ikind) :: ioerr, ntx, ntx_min
    REAL(rkind) :: tl_min,tl_min_old


    IF(ntxumax.NE.1)THEN
       DO
          IF(time_slc.LT.tl(1) .OR. time_slc.GT.tl(ntxumax))THEN
             WRITE(6,'(A,F9.5,A,F9.5,A)')             &
             &  '# Input arbitrary time between: ',   &
             &   tl(1),' sec. - ',tl(ntxumax),' sec.'
             READ(5,*,IOSTAT=ioerr) time_slc
             IF(ioerr.NE.0 .OR. time_slc.LT.tl(1)) CYCLE
          END IF
          EXIT
       END DO

       IF(time_slc.GT.tl(ntxumax))THEN
          time_slc = tl(ntxumax)
          WRITE(6,'(A,F9.5,A,A,F9.5,A)')                    &
          &    ' Designated time: ',time_slc,' sec.',       &
          &    ' has been replaced by ',tl(ntxumax),' sec.'
          ntsl = ntxumax
          RETURN
       ELSE IF(time_slc==tl(ntxumax))THEN
          ntsl = ntxumax
          RETURN
       END IF

       tl_min = tl(ntxumax)
       DO ntx = 1, ntxumax
          IF(ABS(tl(ntx)-time_slc) .LE. 1.d-5)THEN
             ntsl = ntx
             EXIT
          END IF

          tl_min_old = tl_min
          tl_min     = MIN(ABS(tl(ntx)-time_slc), tl_min)
          IF(tl_min_old==tl_min .OR. ntx==ntxumax)THEN
             ntx_min  = ntx - 1
             IF(ntx==ntxumax) ntx_min = ntxumax
             WRITE(6,'(A,F9.5,A,A,F9.5,A)')                    &
             &    ' Designated time: ',time_slc,' sec.',       &
             &    ' has been replaced by ',tl(ntx_min),' sec.'

             time_slc = tl(ntx_min)
             ntsl     = ntx_min
             EXIT
          END IF             
       END DO
    END IF
       
    RETURN
  END SUBROUTINE tr_uf_time_slice

END MODULE trufin
