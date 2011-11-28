MODULE trcoef

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_coef

  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: ADDW,AKDW,AVDW,AVKDW
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: ADDWD,ADDWP,AKDWD,AKDWP

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_coef
    USE trcomm, ONLY: rkind,ikind,mdlprv,nrmax,nsamax,dfa,pha,mdld
    USE trglf, ONLY: tr_glf23
    IMPLICIT NONE

    SELECT CASE(mdld)
    CASE(60,61)
       ALLOCATE(ADDW(nrmax,nsamax))
       ALLOCATE(AKDW(nrmax,nsamax))
       ALLOCATE(AVDW(nrmax,nsamax))
       ALLOCATE(AVKDW(nrmax,nsamax))
       ALLOCATE(ADDWD(nrmax,nsamax,nsamax))
       ALLOCATE(ADDWP(nrmax,nsamax,nsamax))
       ALLOCATE(AKDWD(nrmax,nsamax,nsamax))
       ALLOCATE(AKDWP(nrmax,nsamax,nsamax))
       CALL tr_glf23(ADDW,ADDWD,ADDWP,AKDW,AKDWD,AKDWP,AVDW,AVKDW)
    END SELECT
    

    CALL tr_calc_dfa
    CALL tr_calc_vca
    CALL tr_calc_exa
    CALL tr_calc_pha
    IF(mdlprv /= 0) CALL Pereverzev_method

    RETURN
  END SUBROUTINE tr_coef

! ----- calculate diffusion coefficient -----

  SUBROUTINE tr_calc_dfa
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,mdld,ltcr,d0,d1,rg,rt,dfa, &
         nsa_neq,lt_save
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ, nsa
    REAL(rkind) :: LT
    REAL(rkind),DIMENSION(neqmax,nrmax):: DFG

    SELECT CASE(MDLD)
       ! --- FLAT PROFILE ---
    CASE(0)
       DO NR = 1, NRMAX
          DO NEQ = 1, NEQMAX
             DFG(NEQ,NR) = D0
          END DO
          lt_save(nr)=0.d0
       END DO
      
       ! --- Typical stiff model (Pereverzev) ---
    CASE(1)
       DO NR = 1, NRMAX
          DO NEQ = 1, NEQMAX
             nsa=nsa_neq(neq)
             IF(nsa == 0) THEN
                DFG(NEQ,NR) = D0
                lt_save(nr)=0.D0
             ELSE
                IF(RT(nsa,NR)-RT(nsa,NR-1) > 0) THEN
                   LT = 0.D0
                ELSE
                   LT = - (RT(nsa,NR)-RT(nsa,NR-1))/(RG(NR)-RG(NR-1))
                END IF
                lt_save(nr)=lt
                IF(LT > LTCR)THEN
                   DFG(NEQ,NR) = D0 + D1 - D1*LTCR/LT
                ELSE
                   DFG(NEQ,NR) = D0
                END IF
             END IF
          END DO
       END DO
    END SELECT
         
    DO NR = 1, NRMAX
       DO NEQ = 1, NEQMAX
          DFa(NEQ,1:NEQMAX,NR) = 0.D0
          DFa(NEQ,NEQ,NR) = DFG(NEQ,NR)
       END DO
    END DO
    NR = 0
       DO NEQ = 1, NEQMAX
          DFa(NEQ,1:NEQMAX,NR) = DFa(NEQ,1:NEQMAX,NR+1)
       END DO
    RETURN
  END SUBROUTINE tr_calc_dfa

! ----- calculate convenction velocity -----

  SUBROUTINE tr_calc_vca
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,vca
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    NR=0
       DO NEQ = 1, NEQMAX
          VCa(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    DO NR = 1, NRMAX
       DO NEQ = 1, NEQMAX
          VCa(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_calc_vca

! ----- calculate energy transfer -----

  SUBROUTINE tr_calc_exa
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,exa
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    DO NR = 0, NRMAX
       DO NEQ = 1, NEQMAX
          EXa(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_calc_exa

! ----- calculate heat source -----

  SUBROUTINE tr_calc_pha
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,ph0,phs,rg,ra,pha
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ

    DO NR = 0, NRMAX
       DO NEQ = 1, NEQMAX
          PHa(neq,nr) = phs+(ph0-phs)*(1.D0-(RG(NR)/RA)**2)
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_calc_pha

! -----Pereverzev method -----

  SUBROUTINE Pereverzev_method
    USE trcomm, ONLY: ikind,rkind,nrmax,neqmax,ph0,rg,rt_prev, &
         mdlprv,dprv1,dprv2,dfa,vca,nsa_neq
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
                DDB = DPRV2*DFa(NEQ,NEQ,NR)
             CASE(3)
                DDB = DPRV2*DFa(NEQ,NEQ,NR)+DPRV1
             END SELECT
             VDB = DDB * lt / rtave
          ENDIF

          DFa(NEQ,NEQ,NR) = DFa(NEQ,NEQ,NR) + DDB
          VCa(NEQ,NEQ,NR) = VCa(NEQ,NEQ,NR) + VDB
       END DO
    END DO
    NR = 0
       DO NEQ = 1, NEQMAX
          DFa(NEQ,1:NEQMAX,NR) = DFa(NEQ,1:NEQMAX,NR+1)
          VCa(NEQ,1:NEQMAX,NR) = 0.D0
       END DO
    RETURN
  END SUBROUTINE Pereverzev_method
END MODULE trcoef
