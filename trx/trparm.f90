! trparm.f90

MODULE trparm

  PRIVATE
  PUBLIC tr_parm
  PUBLIC tr_nlin

CONTAINS

!     ***********************************************************

!           PARAMETER INPUT

!     ***********************************************************

  SUBROUTINE tr_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

      USE TRCOMM
      USE libkio
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER,INTENT(OUT):: IERR

    1 CALL TASK_PARM(MODE,'TR',KIN,tr_nlin,trplst,IERR)
      IF(IERR.NE.0) RETURN

      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
    END SUBROUTINE tr_parm

!     ****** INPUT NAMELIST ******

    SUBROUTINE tr_nlin(NID,IST,IERR)

      USE trcomm
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NID
      INTEGER,INTENT(OUT):: IST, IERR

      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
                    PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    PROFJ1,PROFJ2,ALP, &
                    DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
                    EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA, &
                    TPRST,CNN,CDW,CWEB,CALF,CNB,CSPRS, &
                    model_chi_tb,model_dp_tb,model_vk_tb,model_vp_tb, &
                    model_chi_nc,model_dp_nc,model_vk_nc,model_vp_nc, &
                    model_eta,model_bs,model_tpfrac, & 
                    factor_chi_tb,factor_dp_tb,factor_vk_tb,factor_vp_tb, &
                    factor_chi_nc,factor_dp_nc,factor_vk_nc,factor_vp_nc, &
                    factor_eta,factor_bs, &
                    MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDEDGE,MDLIMP, &
                    MDLXP,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05, &
                    PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG,PNBCD,MDLNB, &
                    PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC, &
                    PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH, &
                    PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC, &
                    PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL, &
                    MDLPR,SYNCABS,SYNCSELF, &
                    KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS, &
                    MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI, &
                    NSMAX,NSZMAX,NSNMAX, &
                    KUFDIR,KUFDEV,KUFDCG, &
                    TIME_INT,MODEP,MDNI,MDLJQ,MDLPCK, &
                    MDLPSC,NPSCMAX,PSCTOT,PSCR0,PSCRW,NSPSC, &
                    PBSCD,MDLCD

      IF(NID.GE.0) THEN
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
      ENDIF
      IST=0
      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
    END SUBROUTINE tr_nlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE trplst

      WRITE(6,*) '# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA'
      WRITE(6,*) '(PA,PZ,PN,PNS,PT,PTS:NSM)'
      WRITE(6,*) 'PNC,PNFE,PNNU,PNNUS'
      WRITE(6,*) 'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'
      WRITE(6,*) 'PROFJ1,PROFJ2,ALP'
      WRITE(6,*) 'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'
      WRITE(6,*) 'EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA'
      WRITE(6,*) 'TPRST,CNN,CDW,CWEB,CALF,CNB,CSPRS'
      WRITE(6,*) 'model_chi_tb,model_dp_tb,model_vk_tb,model_vp_tb'
      WRITE(6,*) 'model_chi_nc,model_dp_nc,model_vk_nc,model_vp_nc'
      WRITE(6,*) 'model_eta,model_bs,model_tpfrac'
      WRITE(6,*) 'factor_chi_tb,factor_dp_tb,factor_vk_tb,factor_vp_tb'
      WRITE(6,*) 'factor_chi_nc,factor_dp_nc,factor_vk_nc,factor_vp_nc'
      WRITE(6,*) 'factor_eta,factor_bs'
      WRITE(6,*) 'MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDEDGE,MDLIMP'
      WRITE(6,*) 'MDLXP,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05'
      WRITE(6,*) 'PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG,PNBCD,MDLNB'
      WRITE(6,*) 'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'
      WRITE(6,*) 'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'
      WRITE(6,*) 'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'
      WRITE(6,*) 'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL'
      WRITE(6,*) 'MDLPR,SYNCABS,SYNCSELF'
      WRITE(6,*) 'KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS'
      WRITE(6,*) 'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEQI'
      WRITE(6,*) 'NSMAX,NSZMAX,NSNMAX'
      WRITE(6,*) 'KUFDIR,KUFDEV,KUFDCG'
      WRITE(6,*) 'TIME_INT,MODEP,MDNI,MDLJQ,MDLPCK'
      WRITE(6,*) 'MDLPSC,NPSCMAX,PSCTOT,PSCR0,PSCRW,NSPSC'
      RETURN
    END SUBROUTINE trplst

!     ***** CHECK INPUT PARAMETERS *****

    SUBROUTINE trchek(IERR)

      USE TRCOMM, ONLY : NGTSTP, NRMAX, NTM, NTMAX
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: IERR


      IERR=0

      IF(NRMAX.LT.1) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX =',NRMAX
         IERR=1
      ENDIF

      IF(NTMAX.LT.0.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF

      RETURN
    END SUBROUTINE trchek

END MODULE trparm
