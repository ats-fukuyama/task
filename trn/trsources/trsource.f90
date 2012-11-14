MODULE trsource

  USE trcomm,ONLY: ikind,rkind,neqmax,nrmax

  PRIVATE
  PUBLIC tr_source1,tr_source2

CONTAINS

  SUBROUTINE tr_source1
! ------------------------------------------------------------------------
!   *** calculation source terms which do not need non-linear calculations
! ------------------------------------------------------------------------
    USE trcomm, ONLY: t,time_slc,mdluf,mdlsrc, &
                      poh,pnb,pec,pic,plh,pibw,pnf,prl,snb,spl,swl
    USE trufin, ONLY: tr_ufin_source
    IMPLICIT NONE
    REAL(rkind)    :: time
    INTEGER(ikind) :: ierr

    poh(1:2,0:nrmax)  = 0.d0
    pnb(1:2,0:nrmax)  = 0.d0
    pec(1:2,0:nrmax)  = 0.d0
    pic(1:2,0:nrmax)  = 0.d0
    plh(1:2,0:nrmax)  = 0.d0
    pibw(1:2,0:nrmax) = 0.d0
    pnf(1:2,0:nrmax)  = 0.d0
    prl(1:2,0:nrmax)  = 0.d0
    snb(1:2,0:nrmax)  = 0.d0
    spl(1:2,0:nrmax)  = 0.d0
    swl(1:2,0:nrmax)  = 0.d0


    IF(mdluf==1)THEN
       time = time_slc
    ELSE IF(mdluf==2)THEN
       time = t
    END IF

    SELECT CASE(mdlsrc)
    CASE(1) ! simple model
       CALL tr_src_simple

    CASE(2) ! conventional heating calculation
       CALL tr_src_simple

    CASE(6) ! only initial condition
       IF(t /= 0.d0) RETURN
       ! experimental data
       CALL tr_ufin_source(time,0,ierr)

    CASE(7) ! successive input
       CALL tr_ufin_source(time,0,ierr)

    CASE(8:9)
       ! other module or code
    END SELECT


    RETURN
  END SUBROUTINE tr_source1


  SUBROUTINE tr_source2
! ------------------------------------------------------------------------
!   *** calculation source terms which need non-linear calculations
!   *** substitution all contribution to calculation variable 'str'
! ------------------------------------------------------------------------
    USE trcomm,ONLY: rkev,nsa_neq,nva_neq,mdlsrc,str,htr, &
         poh,pnb,pec,pic,plh,pibw,pnf,prl,snb,spl,swl,    &
         jbs_nc,jex_nc,jcd_nb,jcd_ec,jcd_ic,jcd_lh
    IMPLICIT NONE

    INTEGER(ikind) :: neq,nsa,nva

    str(1:neqmax,0:nrmax)=0.D0
    htr(1:neqmax,0:nrmax)=0.D0 ! external driven current


    CALL tr_src_energy_ohm  ! calculate ohmic heating

    SELECT CASE(mdlsrc)
    CASE(0:1) ! conventional heating
       ! ----- energy -----

       ! ----- particle -----

    CASE(6)

    CASE(7)
       
    END SELECT



    ! substitution each contribution of sources
    ! (from both tr_source1 and tr_source2)
    DO neq = 1, neqmax
       nva = nva_neq(neq)
       nsa = nsa_neq(neq)

       SELECT CASE(nva)
       CASE(0) ! magnetic
          htr(1,0:nrmax) =  jbs_nc(0:nrmax) + jex_nc(0:nrmax) &
                          + jcd_nb(0:nrmax) + jcd_ec(0:nrmax) &
                          + jcd_ic(0:nrmax) + jcd_lh(0:nrmax)
       CASE(1) ! density
          IF(nsa == 1)THEN
             str(neq,0:nrmax) = ( snb(1,0:nrmax)  &
                                + spl(1,0:nrmax)  &
                                + swl(1,0:nrmax))*1.d-20
          ELSE IF(nsa == 2)THEN
             str(neq,0:nrmax) = ( snb(2,0:nrmax)  &
                                + spl(2,0:nrmax)  &
                                + swl(2,0:nrmax))*1.d-20
          END IF
       CASE(2) ! toroidal rotation
       CASE(3) ! energy
          IF(nsa == 1)THEN ! energy, electron
             str(neq,0:nrmax) = ( poh(1,0:nrmax)  &
                                + pnb(1,0:nrmax)  &
                                + pec(1,0:nrmax)  &
                                + pic(1,0:nrmax)  &
                                + plh(1,0:nrmax)  &
                                + pibw(1,0:nrmax) &
                                + pnf(1,0:nrmax)  &
                                - prl(1,0:nrmax))/(rkev*1.d20)
                                
          ELSE IF(nsa == 2)THEN ! energy, ion
             str(neq,0:nrmax) = ( pnb(2,0:nrmax)  &
                                + pec(2,0:nrmax)  &
                                + pic(2,0:nrmax)  &
                                + plh(2,0:nrmax)  &
                                + pibw(2,0:nrmax) &
                                + pnf(2,0:nrmax)  &
                                - prl(2,0:nrmax))/(rkev*1.d20)
          END IF
       END SELECT
    END DO

!    str(2:neqmax,0:nrmax) = str_simple(2:neqmax,0:nrmax)*1.d6

    RETURN
  END SUBROUTINE tr_source2

! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_src_simple
! ------------------------------------------------------------------------
!   simple source model for particle and energy equations
!   for the time being
! ------------------------------------------------------------------------
    USE trcomm, ONLY: rkev,nrmax,nsamax,neqmax,nva_neq,ph0,phs,rhog,ra, &
         dvrho,str_simple,nrd1,nrd2,nrd3,nrd4, &
         pnb,pnb_tot,pnb_r0,pnb_rw,pnb_eng,snb
    IMPLICIT NONE
    INTEGER(ikind) :: nr, neq

    REAL(rkind) :: sum_nb,pnb0,rhom,dvrhom,drhog

    str_simple = 0.d0
    DO nr = 0, nrmax
       DO neq=1,neqmax
          IF(nva_neq(neq) == 3) THEN
             str_simple(neq,nr) = phs+(ph0-phs)*(1.D0-(rhog(nr)/ra)**2)*1.d6
!             str_simple(neq,nr) = phs+(ph0-phs)*(1.D0-(rg(nr)/ra)**2)
!             str_simple(neq,nr) = 0.d0
          END IF
       END DO
    END DO

    ! simple RF : Gaussian profile : only energy
    sum_nb = 0.d0
    DO nr = 1, nrmax
       rhom   = 0.5d0*(rhog(nr)+rhog(nr-1))
       dvrhom = 0.5d0*(dvrho(nr)+dvrho(nr-1))
       drhog  = rhog(nr) - rhog(nr-1)

       sum_nb = sum_nb &
              + DEXP( -((ra*rhom-pnb_r0)/pnb_rw)**2 ) *dvrhom*drhog
    ENDDO
    pnb0 = pnb_tot*1.d6 / sum_nb
    DO nr = 0, nrmax
       pnb(1,nr) = 0.5d0*pnb0*DEXP(-((ra*rhog(nr)-pnb_r0)/pnb_rw)**2)
       pnb(2,nr) = 0.5d0*pnb0*DEXP(-((ra*rhog(nr)-pnb_r0)/pnb_rw)**2)
    ENDDO

!    nrd1(0:nrmax) = eta(0:nrmax)
!    nrd2(0:nrmax) = ezoh(0:nrmax)
!    nrd3(0:nrmax) = poh(0:nrmax)
!    nrd4(0:nrmax) = joh(0:nrmax)

    RETURN
  END SUBROUTINE tr_src_simple


  SUBROUTINE tr_src_energy_ohm
! ------------------------------------------------------------------------
!   ohmic heating
! ------------------------------------------------------------------------
    USE trcomm, ONLY: poh,joh,ezoh,eta
    IMPLICIT NONE

    INTEGER(ikind) :: nr

    ! ohmic heating [W/m^3]
    DO nr = 0, nrmax
       ezoh(nr)  = eta(nr)*joh(nr)
       poh(1,nr) = ezoh(nr)*joh(nr)    
    END DO

    RETURN
  END SUBROUTINE tr_src_energy_ohm

END MODULE trsource
