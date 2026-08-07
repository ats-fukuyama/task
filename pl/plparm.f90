! plparm.f90

MODULE plparm

  PRIVATE
  PUBLIC pl_parm
  PUBLIC pl_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

    SUBROUTINE pl_parm(MODE,KIN,IERR)

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

      USE libkio
      USE plcomm, ONLY: nsmax,pt,ptpr,ptpp,pu,pupr,pupp,rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN):: mode
      CHARACTER(LEN=*),INTENT(IN)::  kin
      INTEGER,INTENT(OUT):: ierr
      REAL(rkind):: pt_save(nsmax),pu_save(nsmax)
      INTEGER:: ns

1     CONTINUE
      DO ns=1,nsmax
         pt_save(ns)=pt(ns)
         pu_save(ns)=pu(ns)
         pt(ns)=0.d0
         pu(ns)=0.d0
      END DO
      CALL TASK_PARM(MODE,'PL',KIN,plnlin,plplst,IERR)
      IF(IERR.NE.0) THEN
         DO ns=1,nsmax
            pt(ns)=pt_save(ns)
            pu(ns)=pu_save(ns)
         END DO
         RETURN
      END IF

      DO ns=1,nsmax
         IF(pt(ns).EQ.0.d0) THEN
            pt(ns)=(ptpr(ns)+2.D0*ptpp(ns))/3.D0
         ELSE
            ptpr(ns)=pt(ns)
            ptpp(ns)=pt(ns)
         END IF
         IF(pu(ns).EQ.0.D0) THEN
            pu(ns)=pupr(ns)
         ELSE
            pupr(ns)=pu(ns)
            pupp(ns)=0.D0
         END IF
      END DO
      
      CALl plcheck(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
    END subroutine pl_parm

!     ****** INPUT NAMELIST ******

    SUBROUTINE plnlin(NID,IST,IERR)

      use plcomm_parm

      implicit none
      integer,intent(in) :: NID
      integer,intent(out) :: IST,IERR
      INTEGER:: NS

      NAMELIST /PL/ &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOIL,ZCOIL,BCOIL,NCOILMAX, &
           NSMAX,NPA,PA,PZ, &
           PN,PNS,PNM,PT,PTS,PTM,PTPR,PTPP,PU,PUS,PUM,PUPR,PUPP,PNUC,PZCL, &
           PROFN1,PROFN2,PROFN3,PROFT1,PROFT2,PROFT3,PROFU1,PROFU2,PROFU3, &
           nsfmax,nszmax,nsnmax,nstmax, &
           ID_NS,KID_NS, &
           RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
           PPN0,PTN0,RF_PL,BAXIS_SCALED, &
           r_corner,z_corner, &
           br_corner,bz_corner,bt_corner, &
           pn_corner,ptpr_corner,ptpp_corner, &
           profn_travis_g,profn_travis_h,profn_travis_p,profn_travis_q, &
           profn_travis_w,proft_travis_g,proft_travis_h,proft_travis_p, &
           proft_travis_q,proft_travis_w, &
           MODELG,MODELB,MODELQ,model_coll,MODEL_NPROF, &
           model_prof,model_prof_time,model_sigv, &
           model_eqdsk_psi, &
           RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
           knam_profg_TOTAL,knam_profm_TOTAL,knam_coil, &
           MODEFR,MODEFW,IDEBUG,mdlplw

      READ(NID,PL,IOSTAT=IST,ERR=9800,END=9900)

      IF(MODEL_PROF.EQ.0) THEN
         DO NS=2,NSMAX
            PROFN1(NS)=PROFN1(1)
            PROFN2(NS)=PROFN2(1)
            PROFN3(NS)=PROFN3(1)
            PROFT1(NS)=PROFT1(1)
            PROFT2(NS)=PROFT2(1)
            PROFT3(NS)=PROFT3(1)
            PROFU1(NS)=PROFU1(1)
            PROFU2(NS)=PROFU2(1)
            PROFU3(NS)=PROFU3(1)
         END DO
      END IF

      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
    END SUBROUTINE plnlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE plplst

      implicit none

      WRITE(6,*) '&PL : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'
      WRITE(6,*) 'RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOI,ZCOIL,BCOIL,NCOILMAX,'
      WRITE(6,*) 'NSMAX,PA,PZ,PN,PNS,PNM,PT,PTS,PTPR,PTPP,PU,PUS,PUPR,PUPP,'
      WRITE(6,*) 'PROFN1,PROFN2,PROFN3,PROFT1,PROFT2,PROFT3,PROFU1,PROFU2,PROFT3'
      WRITE(6,*) 'PNUC,PZCL,ID_NS,KID_NS,'
      WRITE(6,*) 'r_corner,z_corner,br_corner,bz_corner,bt_corner,'
      WRITE(6,*) 'pn_corner,ptpr_corner,ptpp_corner,'
      WRITE(6,*) 'profn_travis_g,profn_travis_h,profn_travis_p,'
      WRITE(6,*) 'profn_travis_q,profn_travis_w,proft_travis_g,'
      WRITE(6,*) 'proft_travis_h,proft_travis_p,proft_travis_q,'
      WRITE(6,*) 'proft_travis_w,'
      WRITE(6,*) 'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'
      WRITE(6,*) 'PPN0,PTN0,RFCL,BAXIS_SCALED,'
      WRITE(6,*) 'MODELG,MODELB,MODELQ,'
      WRITE(6,*) 'model_coll,MODEL_NPROF,RHOGMN,RHOGMX,'
      WRITE(6,*) 'model_prof,model_prof_time,model_sigv,'
      WRITE(6,*) 'model_eqdsk_psi,'
      WRITE(6,*) 'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2,'
      WRITE(6,*) 'knam_profg_TOTAL,knam_profm_TOTAL,knam_coil,'
      WRITE(6,*) 'MODEFW,MODEFR,IDEBUG,mdlplw'
      RETURN
    END SUBROUTINE plplst

!     ****** CHECK INPUT PARAMETER ******

    SUBROUTINE plcheck(IERR)

      use plcomm
      implicit none
      integer:: IERR,NS
      real(rkind):: RHOE0,RHOES

      IERR=0

      IF(MODELG.NE.3) THEN
         IF((MODELG.LT.0).OR.(MODELG.GT.11)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELG: MODELG=',MODELG
            IERR=1
         ENDIF
         IF((model_prof.LT.0).OR.(model_prof.GT.41)) THEN
            WRITE(6,*) 'XX plcheck: INVALID model_prof: model_prof=',model_prof
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ELSE
         IF((model_prof.LT.0).OR.(model_prof.GT.41)) THEN
            WRITE(6,*) 'XX plcheck: INVALID model_prof: model_prof=',model_prof
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1).AND.(MODELG.NE.3)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ENDIF

      RHOE0=0.D0
      RHOES=0.D0
      DO NS=1,NSMAX
         RHOE0=RHOE0+PZ(NS)*PN(NS)
         RHOES=RHOES+PZ(NS)*PNS(NS)
      ENDDO
      IF(ABS(RHOE0).GT.1.D-10) THEN
         WRITE(6,*) 'XX PLPARM: CHARGE NEUTRALITY ERROR AT CENTER'
         IERR=1
      ENDIF
      IF(ABS(RHOES).GT.1.D-10) THEN
         WRITE(6,*) 'XX PLPARM: CHARGE NEUTRALITY ERROR AT SURFACE'
         IERR=1
      ENDIF

      RETURN
    END SUBROUTINE plcheck

    ! --- Broadcast pl input parameters ---

  SUBROUTINE pl_broadcast

    USE plcomm_parm
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata
    INTEGER:: NS

    idata( 1)=NSMAX
    idata( 2)=MODELG
    idata( 3)=MODELB
    idata( 4)=model_prof
    idata( 5)=MODELQ
    idata( 6)=IDEBUG
    idata( 7)=MODEFR
    idata( 8)=MODEFW
    idata( 9)=mdlplw
    idata(10)=model_coll
    idata(11)=model_prof_time
    idata(12)=MODEL_NPROF
    idata(13)=NCOILMAX
    idata(14)=model_sigv

    CALL mtx_broadcast_integer(idata,14)
    
    NSMAX=idata( 1)
    MODELG=idata( 2)
    MODELB=idata( 3)
    model_prof=idata( 4)
    MODELQ=idata( 5)
    IDEBUG=idata( 6)
    MODEFR=idata( 7)
    MODEFW=idata( 8)
    mdlplw=idata( 9)
    model_coll=idata(10)
    model_prof_time=idata(11)
    MODEL_NPROF=idata(12)
    NCOILMAX=idata(13)
    model_sigv=idata(14)

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
    rdata(11)=RMIR
    rdata(12)=ZBB
    rdata(13)=Hpitch1
    rdata(14)=Hpitch2
    rdata(15)=RRCH
    rdata(16)=RHOMIN
    rdata(17)=QMIN
    rdata(18)=RHOEDG
    rdata(19)=RHOGMN
    rdata(20)=RHOGMX
    rdata(21)=PPN0
    rdata(22)=PTN0
    rdata(23)=RF_PL
    rdata(24)=profn_travis_g
    rdata(25)=profn_travis_h
    rdata(26)=profn_travis_p
    rdata(27)=profn_travis_q
    rdata(28)=profn_travis_w
    rdata(29)=proft_travis_g
    rdata(30)=proft_travis_h
    rdata(31)=proft_travis_p
    rdata(32)=proft_travis_q
    rdata(33)=proft_travis_w
    rdata(34)=BAXIS_SCALED
    
    CALL mtx_broadcast_real8(rdata,34)
    
    RR=rdata( 1)
    RA=rdata( 2)
    RB=rdata( 3)
    RKAP=rdata( 4)
    RDLT=rdata( 5)
    BB=rdata( 6)
    Q0=rdata( 7)
    QA=rdata( 8)
    RIP=rdata( 9)
    PROFJ=rdata(10)
    RMIR=rdata(11)
    ZBB=rdata(12)
    Hpitch1=rdata(13)
    Hpitch2=rdata(14)
    RRCH=rdata(15)
    RHOMIN=rdata(16)
    QMIN=rdata(17)
    RHOEDG=rdata(18)
    RHOGMN=rdata(19)
    RHOGMX=rdata(20)
    PPN0=rdata(21)
    PTN0=rdata(22)
    RF_PL=rdata(23)
    profn_travis_g=rdata(24)
    profn_travis_h=rdata(25)
    profn_travis_p=rdata(26)
    profn_travis_q=rdata(27)
    profn_travis_w=rdata(28)
    proft_travis_g=rdata(29)
    proft_travis_h=rdata(30)
    proft_travis_p=rdata(31)
    proft_travis_q=rdata(32)
    proft_travis_w=rdata(33)
    BAXIS_SCALED=rdata(34)

    CALL mtx_broadcast_integer(NPA,NSMAX)
    CALL mtx_broadcast_integer(ID_NS,NSMAX)
    DO NS=1,NSMAX
       CALL mtx_broadcast_character(KID_NS(NS),4)
    END DO
    CALL mtx_broadcast_real8(PA,NSMAX)
    CALL mtx_broadcast_real8(PZ,NSMAX)
    CALL mtx_broadcast_real8(PN,NSMAX)
    CALL mtx_broadcast_real8(PNS,NSMAX)
    CALL mtx_broadcast_real8(PTPR,NSMAX)
    CALL mtx_broadcast_real8(PTPP,NSMAX)
    CALL mtx_broadcast_real8(PTS,NSMAX)
    CALL mtx_broadcast_real8(PU,NSMAX)
    CALL mtx_broadcast_real8(PUS,NSMAX)
    CALL mtx_broadcast_real8(PUPR,NSMAX)
    CALL mtx_broadcast_real8(PUPP,NSMAX)
    CALL mtx_broadcast_real8(RHOITB,NSMAX)
    CALL mtx_broadcast_real8(PNITB,NSMAX)
    CALL mtx_broadcast_real8(PTITB,NSMAX)
    CALL mtx_broadcast_real8(PUITB,NSMAX)
    CALL mtx_broadcast_real8(PROFN1,NSMAX)
    CALL mtx_broadcast_real8(PROFN2,NSMAX)
    CALL mtx_broadcast_real8(PROFN3,NSMAX)
    CALL mtx_broadcast_real8(PROFT1,NSMAX)
    CALL mtx_broadcast_real8(PROFT2,NSMAX)
    CALL mtx_broadcast_real8(PROFT3,NSMAX)
    CALL mtx_broadcast_real8(PROFU1,NSMAX)
    CALL mtx_broadcast_real8(PROFU2,NSMAX)
    CALL mtx_broadcast_real8(PROFU3,NSMAX)
    CALL mtx_broadcast_real8(PNUC,NSMAX)
    CALL mtx_broadcast_real8(PZCL,NSMAX)

    CALL mtx_broadcast_real8(RCOIL,NCOILMAX)
    CALL mtx_broadcast_real8(ZCOIL,NCOILMAX)
    CALL mtx_broadcast_real8(BCOIL,NCOILMAX)

    CALL mtx_broadcast_real8(r_corner,3)
    CALL mtx_broadcast_real8(z_corner,3)
    CALL mtx_broadcast_real8(br_corner,3)
    CALL mtx_broadcast_real8(bz_corner,3)
    CALL mtx_broadcast_real8(bt_corner,3)
    DO NS=1,NSMAX
       CALL mtx_broadcast_real8(pn_corner(1:3,NS),3)
       CALL mtx_broadcast_real8(ptpr_corner(1:3,NS),3)
       CALL mtx_broadcast_real8(ptpp_corner(1:3,NS),3)
    END DO

    CALL mtx_broadcast_character(KNAMEQ,80)
    CALL mtx_broadcast_character(KNAMWR,80)
    CALL mtx_broadcast_character(KNAMFP,80)
    CALL mtx_broadcast_character(KNAMWM,80)
    CALL mtx_broadcast_character(KNAMPF,80)
    CALL mtx_broadcast_character(KNAMFO,80)
    CALL mtx_broadcast_character(KNAMTR,80)
    CALL mtx_broadcast_character(KNAMEQ2,80)
    CALL mtx_broadcast_character(knam_coil,80)
    CALL mtx_broadcast_character(knam_profg_TOTAL,128)
    CALL mtx_broadcast_character(knam_profm_TOTAL,128)
  END SUBROUTINE pl_broadcast

 END MODULE plparm
