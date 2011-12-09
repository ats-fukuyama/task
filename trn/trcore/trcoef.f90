MODULE trcoef

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_coef

  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: ADDW,AKDW,AVDW,AVKDW
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: ADDWD,ADDWP,AKDWD,AKDWP

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_coef
    USE trcomm, ONLY: rkind,ikind,mdlprv,nrmax,nsamax,dtr,str,mdld
    USE trglf23, ONLY: tr_glf23
    IMPLICIT NONE

    SELECT CASE(mdld)
    CASE(60,61)
       ALLOCATE(vtrglf(nrmax,3*nsamax))
       ALLOCATE(dtrglf(nrmax,3*nsamax,3*nsamax))
       CALL tr_glf23(dtrglf,vtrglf)
       DTR(1:nrmax,1:3*nsamax,1:3*nsamax) &
            =DTR(1:nrmax,1:3*nsamax,1:3*nsamax) &
            +DTR(1:nrmax,1:3*nsamax,1:3*nsamax)
    END SELECT
    

    CALL calc_dtr
    CALL calc_vtr
    CALL calc_ctr
    CALL calc_htr
    CALL calc_str
    IF(mdlprv /= 0) CALL Pereverzev_method

    RETURN
  END SUBROUTINE tr_coef

! ----- calculate diffusion coefficient -----

  SUBROUTINE calc_dtr
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,mdld,ltcr,d0,d1,rg,rt,dtr, &
         nsa_neq,lt_save
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ, nsa
    REAL(rkind) :: LT
    REAL(rkind),DIMENSION(neqmax,nrmax):: DTR_DIAG

    SELECT CASE(MDLD)
       ! --- FLAT PROFILE ---
    CASE(0)
       DO NR = 1, NRMAX
          DO NEQ = 1, NEQMAX
             nsa=nsa_neq(neq)
             DTR_DIAG(NEQ,NR) = D0
             if(nsa /= 0) lt_save(nsa,nr)=0.d0
          END DO
       END DO
      
       ! --- Typical stiff model (Pereverzev) ---
    CASE(1)
       DO NR = 1, NRMAX
          DO NEQ = 1, NEQMAX
             nsa=nsa_neq(neq)
             IF(nsa == 0) THEN
                DTR_DIAG(NEQ,NR) = D0
             ELSE
                IF(RT(nsa,NR)-RT(nsa,NR-1) > 0) THEN
                   LT = 0.D0
                ELSE
                   LT = - (RT(nsa,NR)-RT(nsa,NR-1))/(RG(NR)-RG(NR-1))
                END IF
                lt_save(nr,nsa)=lt
                IF(LT > LTCR)THEN
                   DTR_DIAG(NEQ,NR) = D0 + D1 - D1*LTCR/LT
                ELSE
                   DTR_DIAG(NEQ,NR) = D0
                END IF
             END IF
          END DO
       END DO
    END SELECT
         
    DO NR = 1, NRMAX
       DO NEQ = 1, NEQMAX
          dtr(NEQ,1:NEQMAX,NR) = 0.D0
          dtr(NEQ,NEQ,NR) = DTR_DIAG(NEQ,NR)
       END DO
    END DO
    NR = 0
       DO NEQ = 1, NEQMAX
          dtr(NEQ,1:NEQMAX,NR) = dtr(NEQ,1:NEQMAX,NR+1)
       END DO
    RETURN
  END SUBROUTINE calc_dtr

! ----- calculate convenction velocity -----

  SUBROUTINE calc_vtr
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,vtr
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    NR=0
       DO NEQ = 1, NEQMAX
          vtr(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    DO NR = 1, NRMAX
       DO NEQ = 1, NEQMAX
          vtr(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    END DO
    RETURN
  END SUBROUTINE calc_vtr

! ----- calculate energy transfer -----

  SUBROUTINE calc_ctr
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,ctr
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    DO NR = 0, NRMAX
       DO NEQ = 1, NEQMAX
          ctr(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    END DO
    RETURN
  END SUBROUTINE calc_ctr

! ----- calculate current source -----

  SUBROUTINE calc_htr
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,rg,ra,htr
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    DO NR = 0, NRMAX
       DO NEQ = 1, NEQMAX
          htr(neq,nr) = 0.D0
       END DO
    END DO
    RETURN
  END SUBROUTINE calc_htr

! ----- calculate heat source -----

  SUBROUTINE calc_str
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,ph0,phs,rg,ra,str
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    DO NR = 0, NRMAX
       DO NEQ = 1, NEQMAX
          str(neq,nr) = phs+(ph0-phs)*(1.D0-(RG(NR)/RA)**2)
       END DO
    END DO
    RETURN
  END SUBROUTINE calc_str

! -----Pereverzev method -----

  SUBROUTINE Pereverzev_method
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,ph0,rg,rt_prev, &
         mdlprv,dprv1,dprv2,dtr,vtr,nsa_neq
    IMPLICIT NONE
    REAL(rkind) :: VDB,DDB,LT,drt,rtave
    INTEGER(ikind) :: NR,NEQ,NSA
      
    DO NR = 1, NRMAX
       DO NEQ = 1, NEQMAX
          nsa=nsa_neq(neq)
          IF(nsa == 0) THEN
             DDB=0.D0
             VDB=0.D0
          ELSE
             drt=rt_prev(nsa,nr)-rt_prev(nsa,nr-1)
             rtave=0.5D0*(rt_prev(nsa,nr)+rt_prev(nsa,nr-1))
             IF(drt > 0.D0) THEN
                lt = 0.D0
             ELSE
                lt = - drt/(rg(nr)-rg(nr-1))
             END IF
             SELECT CASE(MDLPRV)
             CASE(1)
                DDB = DPRV1
             CASE(2)
                DDB = DPRV2*dtr(NEQ,NEQ,NR)
             CASE(3)
                DDB = DPRV2*dtr(NEQ,NEQ,NR)+DPRV1
             END SELECT
             VDB = DDB * lt / rtave
          ENDIF

          dtr(NEQ,NEQ,NR) = dtr(NEQ,NEQ,NR) + DDB
          vtr(NEQ,NEQ,NR) = vtr(NEQ,NEQ,NR) + VDB
       END DO
    END DO
    NR = 0
       DO NEQ = 1, NEQMAX
          dtr(NEQ,1:NEQMAX,NR) = dtr(NEQ,1:NEQMAX,NR+1)
          vtr(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    RETURN
  END SUBROUTINE Pereverzev_method
END MODULE trcoef
