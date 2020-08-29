! wmparm.f90

MODULE wmparm

  PRIVATE
  PUBLIC wm_parm,wm_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE wm_parm(mode,kin,ierr)
    
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

    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN):: kin
    INTEGER,INTENT(OUT):: ierr
    
1   CALL TASK_PARM(mode,'WM',kin,wm_nlin,wm_plst,ierr)
    IF(ierr.NE.0) RETURN

    CALl wm_check(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1
    IF(ierr.NE.0) ierr=ierr+100

    RETURN
  END SUBROUTINE wm_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wm_nlin(nid,ist,ierr)

    USE wmcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr

    NAMELIST /WM/ &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           NSMAX,NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,PPN0,PTN0,RF_PL, &
           MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
           RHOGMN,RHOGMX,MODEFR,MODEFW,IDEBUG, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,KNAMTR,KNAMEQ2, &
           MODELP,MODELV,NCMIN,NCMAX,PMAX,EMAX, &
           NRMAX,NTHMAX,NHHMAX,NPHMAX,factor_nth,factor_nhh, &
           NRMAX_DP,NTHMAX_DP,NPMAX_DP,NSAMAX_DP,RHON_MIN,RHON_MAX, &
           MODELP,MODELV,NCMIN,NCMAX,NS_NSA_DP,EPSRT,LMAXRT, &
           NSUMAX,NSWMAX,B0_FACT, &
           RF,RFI,RD,PRFIN,BETAJ,NTH0,NPH0,NHC, &
           NAMAX,AJ,AEWGT,AEWGZ,APH,THJ1,THJ2,PHJ1,PHJ2,ANTANG, &
           NPRINT,NGRAPH,MODELJ,MODELA,MODELM,MDLWMK,PNA,PNAL,PTA,ZEFF, &
           FRMIN,FRMAX,FIMIN,FIMAX,FI0,NGFMAX,NGXMAX,NGYMAX, &
           SCMIN,SCMAX,NSCMAX,LISTEG,FRINI,FIINI,DLTNW,EPSNW,LMAXNW,LISTNW, &
           MODENW,NCONT,ILN1,IBL1,ICL1,ILN2,IBL2,ICL2,WAEMIN,WAEMAX,nthgmax

    ierr=0

    READ(nid,WM,IOSTAT=ist,ERR=9800,END=9900)
    IF(ist.NE.0) THEN
       WRITE(6,'(A,I5)') 'XX wm_nlin: READ ERROR: IOSTAT=',ist
       ierr=ist
       RETURN
    END IF
    RETURN

9800 IERR=8
    WRITE(6,'(A)') 'XX wm_nlin: READ ERROR'
    RETURN
9900 IERR=9
    WRITE(6,'(A)') 'XX wm_nlin: READ END OF FILE'
    RETURN
  END SUBROUTINE wm_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wm_plst

    WRITE(6,'(A)') &
         '# &WM: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,', &
         '       NSMAX,NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,', &
         '       PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,', &
         '       RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,', &
         '       PPN0,PTN0,RF_PL,', &
         '       MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,', &
         '       RHOGMN,RHOGMX,MODEFR,MODEFW,IDEBUG,', &
         '       KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,', &
         '       MODELP,MODELV,NCMIN,NCMAX,PMAX,EMAX,', &
         '       NPMAX,NTHMAX,NRMAX,NSAMAX,RHON_MIN,RHON_MAX,', &
         '       NHHMAX,NPHMAX,factor_nth,factor_nhh,', &
         '       NSUMAX,NSWMAX,B0_FACT,', &
         '       RF,RFI,RD,PRFIN,BETAJ,NTH0,NPH0,NHC,', &
         '       NAMAX,AJ,AEWGT,AEWGZ,APH,THJ1,THJ2,PHJ1,PHJ2,ANTANG,', &
         '       NPRINT,NGRAPH,MODELJ,MODELA,MODELM,MDLWMK,', &
         '       PNA,PNAL,PTA,ZEFF,', &
         '       FRMIN,FRMAX,FIMIN,FIMAX,FI0,NGFMAX,NGXMAX,NGYMAX,', &
         '       SCMIN,SCMAX,NSCMAX,LISTEG,FRINI,FIINI,', &
         '       DLTNW,EPSNW,LMAXNW,LISTNW,', &
         '       MODENW,NCONT,ILN1,IBL1,ICL1,ILN2,IBL2,ICL2,WAEMIN,WAEMAX', &
         '       nthgmax'
    RETURN
  END SUBROUTINE wm_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE wm_check(ierr)

    USE wmcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    IF(NSMAX.LT.0.OR.NSMAX.GT.NSM) THEN
       WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NSMAX'
       WRITE(6,*) '                  NSMAX,NSM =',NSMAX,NSM
       IERR=1
    ENDIF

    IF((NAMAX.LT.1).OR.(NAMAX.GT.NAM)) THEN
       WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NAMAX'
       WRITE(6,*) '                  NAMAX,NAM =',NAMAX,NAM         
       IERR=1
    ENDIF

    IF(MODELJ.LT.0) THEN
       WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL MODELJ'
       WRITE(6,*) '                  MODELJ =',MODELJ
       IERR=1
    ENDIF

    RETURN
  END SUBROUTINE wm_check

!     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE wm_broadcast

    USE wmcomm_parm
    USE plparm,ONLY: pl_broadcast
    USE dpparm,ONLY: dp_broadcast
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata
    INTEGER:: NS

    !----- PL input parameters -----     

    CALL pl_broadcast

!----- DP input parameters -----

    CALL dp_broadcast

! --- WM specific input parameters ---

    idata( 1)=NRMAX
    idata( 2)=NTHMAX
    idata( 3)=NHHMAX
    idata( 4)=NPHMAX
    idata( 5)=NSUMAX
    idata( 6)=NSWMAX
    idata( 7)=NTH0
    idata( 8)=NPH0
    idata( 9)=NHC
    idata(10)=NAMAX
    idata(11)=NPRINT
    idata(12)=NGRAPH
    idata(13)=MODELJ
    idata(14)=MODELA
    idata(15)=MODELM
    idata(16)=NGFMAX
    idata(17)=NGXMAX
    idata(18)=NGYMAX
    idata(19)=NSCMAX
    idata(20)=LISTEG
    idata(21)=LMAXNW
    idata(22)=LISTNW
    idata(23)=MODENW
    idata(24)=NCONT
    idata(25)=ILN1
    idata(26)=IBL1
    idata(27)=ICL1
    idata(28)=ILN2
    idata(29)=IBL2
    idata(30)=ICL2
    idata(31)=NTHGMAX
    idata(32)=MDLWMK
    idata(33)=factor_nth
    idata(34)=factor_nhh

    CALL mtx_broadcast_integer(idata,34)

    NRMAX=idata( 1)
    NTHMAX=idata( 2)
    NHHMAX=idata( 3)
    NPHMAX=idata( 4)
    NSUMAX=idata( 5)
    NSWMAX=idata( 6)
    NTH0=idata( 7)
    NPH0=idata( 8)
    NHC=idata( 9)
    NAMAX=idata(10)
    NPRINT=idata(11)
    NGRAPH=idata(12)
    MODELJ=idata(13)
    MODELA=idata(14)
    MODELM=idata(15)
    NGFMAX=idata(16)
    NGXMAX=idata(17)
    NGYMAX=idata(18)
    NSCMAX=idata(19)
    LISTEG=idata(20)
    LMAXNW=idata(21)
    LISTNW=idata(22)
    MODENW=idata(23)
    NCONT=idata(24)
    ILN1=idata(25)
    IBL1=idata(26)
    ICL1=idata(27)
    ILN2=idata(28)
    IBL2=idata(29)
    ICL2=idata(30)
    NTHGMAX=idata(31)
    MDLWMK=idata(32)
    factor_nth=idata(33)
    factor_nhh=idata(34)

    rdata( 1)=RF
    rdata( 2)=RFI
    rdata( 3)=RD
    rdata( 4)=PRFIN
    rdata( 5)=0.D0
    rdata( 6)=PNA
    rdata( 7)=PNAL
    rdata( 8)=PTA
    rdata( 9)=ZEFF
    rdata(10)=0.D0
    rdata(11)=0.D0
    rdata(12)=B0_FACT
    rdata(13)=FRMIN
    rdata(14)=FRMAX
    rdata(15)=FIMIN
    rdata(16)=FIMAX
    rdata(17)=FI0
    rdata(18)=SCMIN
    rdata(19)=SCMAX
    rdata(20)=FRINI
    rdata(21)=FIINI
    rdata(22)=DLTNW
    rdata(23)=WAEMIN
    rdata(24)=WAEMAX

    CALL mtx_broadcast_real8(rdata,24)

    RF=rdata( 1)
    RFI=rdata( 2)
    RD=rdata( 3)
    PRFIN=rdata( 4)
!    BETAJ=rdata( 5)
    PNA=rdata( 6)
    PNAL=rdata( 7)
    PTA=rdata( 8)
    ZEFF=rdata( 9)
!    RMAX=rdata(10)
!    RMIN=rdata(11)
    B0_FACT=rdata(12)
    FRMIN=rdata(13)
    FRMAX=rdata(14)
    FIMIN=rdata(15)
    FIMAX=rdata(16)
    FI0=rdata(17)
    SCMIN=rdata(18)
    SCMAX=rdata(19)
    FRINI=rdata(20)
    FIINI=rdata(21)
    DLTNW=rdata(22)
    WAEMIN=rdata(23)
    WAEMAX=rdata(24)
    factor_nth=rdata(25)
    factor_nhh=rdata(26)

    CALL mtx_broadcast_real8(AJ,NAMAX)
    CALL mtx_broadcast_real8(AEWGT,NAMAX)
    CALL mtx_broadcast_real8(AEWGZ,NAMAX)
    CALL mtx_broadcast_real8(APH,NAMAX)
    CALL mtx_broadcast_real8(THJ1,NAMAX)
    CALL mtx_broadcast_real8(THJ2,NAMAX)
    CALL mtx_broadcast_real8(PHJ1,NAMAX)
    CALL mtx_broadcast_real8(PHJ2,NAMAX)
    CALL mtx_broadcast_real8(ANTANG,NAMAX)
    CALL mtx_broadcast_real8(BETAJ,NAMAX)

  END SUBROUTINE wm_broadcast
END MODULE wmparm
