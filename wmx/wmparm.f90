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
           NRMAX,NTHMAX,NPMAX,NSAMAX,RMIN,RMAX, &
           NHHMAX,NPHMAX,factor_nth,factor_nhh, &
           NSUMAX,NSWMAX,B0_FACT, &
           RF,RFI,RD,PRFIN,BETAJ,NTH0,NPH0,NHC, &
           NAMAX,AJ,AEWGT,AEWGZ,APH,THJ1,THJ2,PHJ1,PHJ2,ANTANG, &
           NPRINT,NGRAPH,MODELJ,MODELA,MODELM,PNA,PNAL,PTA,ZEFF, &
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
         '       NPMAX,NTHMAX,NRMAX,NSAMAX,RMIN,RMAX,', &
         '       NHHMAX,NPHMAX,factor_nth,factor_nhh,', &
         '       NSUMAX,NSWMAX,B0_FACT,', &
         '       RF,RFI,RD,PRFIN,BETAJ,NTH0,NPH0,NHC,', &
         '       NAMAX,AJ,AEWGT,AEWGZ,APH,THJ1,THJ2,PHJ1,PHJ2,ANTANG,', &
         '       NPRINT,NGRAPH,MODELJ,MODELA,MODELM,PNA,PNAL,PTA,ZEFF,', &
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
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata
    INTEGER:: NS

!----- PL input parameters -----     

    idata( 1)=NSMAX
    idata( 2)=MODELG
    idata( 3)=MODELN
    idata( 4)=MODELQ
    idata( 5)=IDEBUG
    idata( 6)=MODEFR
    idata( 7)=MODEFW
    idata( 8)=MODEL_PROF
    idata( 9)=MODEL_PROF

    CALL mtx_broadcast_integer(idata,9)

    NSMAX =idata( 1)
    MODELG=idata( 2)
    MODELN=idata( 3)
    MODELQ=idata( 4)
    IDEBUG=idata( 5)
    MODEFR=idata( 6)
    MODEFW=idata( 7)
    MODEL_PROF=idata( 8)
    MODEL_PROF=idata( 9)

    rdata( 1)=RR
    rdata( 2)=RA
    rdata( 3)=RB
    rdata( 4)=RKAP
    rdata( 5)=RDLT
    rdata( 6)=BB
    rdata( 7)=Q0
    rdata( 8)=QA
    rdata( 9)=RIP
    rdata(10)=PROFJ
    rdata(11)=RHOMIN
    rdata(12)=QMIN
    rdata(13)=RHOEDG
    rdata(14)=RHOGMN
    rdata(15)=RHOGMX
    rdata(16)=PPN0
    rdata(17)=PTN0
    rdata(18)=RF_PL

    CALL mtx_broadcast_real8(rdata,22)

    RR    =rdata( 1)
    RA    =rdata( 2)
    RB    =rdata( 3)
    RKAP  =rdata( 4)
    RDLT  =rdata( 5)
    BB    =rdata( 6)
    Q0    =rdata( 7)
    QA    =rdata( 8)
    RIP   =rdata( 9)
    PROFJ =rdata(10)
    RHOMIN=rdata(11)
    QMIN  =rdata(12)
    RHOEDG=rdata(13)
    RHOGMN=rdata(14)
    RHOGMX=rdata(15)
    PPN0=rdata(16)
    PTN0=rdata(17)
    RF_PL=rdata(18)

    CALL mtx_broadcast_integer(NPA,NSMAX)
    CALL mtx_broadcast_real8(PA,NSMAX)
    CALL mtx_broadcast_real8(PZ,NSMAX)
    CALL mtx_broadcast_real8(PN,NSMAX)
    CALL mtx_broadcast_real8(PNS,NSMAX)
    CALL mtx_broadcast_real8(PTPR,NSMAX)
    CALL mtx_broadcast_real8(PTPP,NSMAX)
    CALL mtx_broadcast_real8(PTS,NSMAX)
    CALL mtx_broadcast_real8(PU,NSMAX)
    CALL mtx_broadcast_real8(PUS,NSMAX)
    CALL mtx_broadcast_real8(RHOITB,NSMAX)
    CALL mtx_broadcast_real8(PNITB,NSMAX)
    CALL mtx_broadcast_real8(PTITB,NSMAX)
    CALL mtx_broadcast_real8(PUITB,NSMAX)
    CALL mtx_broadcast_real8(PROFN1,NSMAX)
    CALL mtx_broadcast_real8(PROFN2,NSMAX)
    CALL mtx_broadcast_real8(PROFT1,NSMAX)
    CALL mtx_broadcast_real8(PROFT2,NSMAX)
    CALL mtx_broadcast_real8(PROFU1,NSMAX)
    CALL mtx_broadcast_real8(PROFU2,NSMAX)
    CALL mtx_broadcast_real8(PZCL,NSMAX)

    CALL mtx_broadcast_integer(ID_NS,NSMAX)
    DO NS=1,NSMAX
       CALL mtx_broadcast_character(KID_NS(NS),2)
    END DO

    CALL mtx_broadcast_character(KNAMEQ,80)
    CALL mtx_broadcast_character(KNAMWR,80)
    CALL mtx_broadcast_character(KNAMFP,80)
    CALL mtx_broadcast_character(KNAMWM,80)
    CALL mtx_broadcast_character(KNAMPF,80)
    CALL mtx_broadcast_character(KNAMFO,80)
    CALL mtx_broadcast_character(KNAMTR,80)
    CALL mtx_broadcast_character(KNAMEQ2,80)

!----- DP input parameters -----

    CALL mtx_broadcast_integer(MODELP,NSMAX)
    CALL mtx_broadcast_integer(MODELV,NSMAX)
    CALL mtx_broadcast_integer(NCMIN,NSMAX)
    CALL mtx_broadcast_integer(NCMAX,NSMAX)
    CALL mtx_broadcast_real8(PMAX,NSMAX)
    CALL mtx_broadcast_real8(EMAX,NSMAX)

! --- WM specific input parameters ---

    idata( 1)=NRMAX
    idata( 2)=NTHMAX
    idata( 3)=NHHMAX
    idata( 4)=NPHMAX
    idata( 5)=0
    idata( 6)=0
    idata( 7)=NPMAX
    idata( 8)=NSAMAX
    idata( 9)=NSUMAX
    idata(10)=NTH0
    idata(11)=NPH0
    idata(12)=NHC
    idata(13)=NAMAX
    idata(14)=NPRINT
    idata(15)=NGRAPH
    idata(16)=MODELJ
    idata(17)=MODELA
    idata(18)=MODELM
    idata(19)=NGFMAX
    idata(20)=NGXMAX
    idata(21)=NGYMAX
    idata(22)=LMAXNW
    idata(23)=LISTNW
    idata(24)=MODENW
    idata(25)=NCONT
    idata(26)=ILN1
    idata(27)=IBL1
    idata(28)=ICL1
    idata(29)=ILN2
    idata(30)=IBL2
    idata(31)=ICL2
    idata(32)=nthgmax

    CALL mtx_broadcast_integer(idata,32)

    NRMAX=idata( 1)
    NTHMAX=idata( 2)
    NHHMAX=idata( 3)
    NPHMAX=idata( 4)
!    NTHMAX_F=idata( 5)
!    NHHMAX_F=idata( 6)
    NPMAX=idata( 7)
    NSAMAX=idata( 8)
    NSUMAX=idata( 9)
    NTH0=idata(10)
    NPH0=idata(11)
    NHC=idata(12)
    NAMAX=idata(13)
    NPRINT=idata(14)
    NGRAPH=idata(15)
    MODELJ=idata(16)
    MODELA=idata(17)
    MODELM=idata(18)
    NGFMAX=idata(19)
    NGXMAX=idata(20)
    NGYMAX=idata(21)
    LMAXNW=idata(22)
    LISTNW=idata(23)
    MODENW=idata(24)
    NCONT=idata(25)
    ILN1=idata(26)
    IBL1=idata(27)
    ICL1=idata(28)
    ILN2=idata(29)
    IBL2=idata(30)
    ICL2=idata(31)
    nthgmax=idata(32)

    rdata( 1)=RF
    rdata( 2)=RFI
    rdata( 3)=RD
    rdata( 4)=PRFIN
    rdata( 5)=BETAJ
    rdata( 6)=PNA
    rdata( 7)=PNAL
    rdata( 8)=PTA
    rdata( 9)=ZEFF
    rdata(10)=RMAX
    rdata(11)=RMIN
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
    rdata(25)=factor_nth
    rdata(26)=factor_nhh

    CALL mtx_broadcast_real8(rdata,26)

    RF=rdata( 1)
    RFI=rdata( 2)
    RD=rdata( 3)
    PRFIN=rdata( 4)
    BETAJ=rdata( 5)
    PNA=rdata( 6)
    PNAL=rdata( 7)
    PTA=rdata( 8)
    ZEFF=rdata( 9)
    RMAX=rdata(10)
    RMIN=rdata(11)
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

  END SUBROUTINE wm_broadcast
END MODULE wmparm
