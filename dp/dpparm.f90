! dpparm.f90

MODULE DPPARM

  PRIVATE
  PUBLIC dp_parm,dp_chek,dpprep,dpprep_local,dp_broadcast,dp_view

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE DP_PARM(MODE,KIN,IERR)

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
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR
    EXTERNAL EQCHEK

    IERR=0

1   CALL TASK_PARM(MODE,'DP',KIN,DPNLIN,DPPLST,IERR)
    IF(IERR.NE.0) RETURN

    CALL EQCHEK(IERR)
    CALL DP_CHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    CALL DPPREP_LOCAL(IERR)
    RETURN
  END SUBROUTINE DP_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE DPNLIN(NID,IST,IERR)

    USE dpcomm_parm_local,DPX=>DP
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR
    INTEGER:: NS

    NAMELIST /DP/ &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOIL,ZCOIL,BCOIL,NCOILMAX, &
           NSMAX,NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PUPR,PUPP,PNUC,PZCL, &
           ID_NS,KID_NS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
           PPN0,PTN0,RF_PL,BAXIS_SCALED, &
           r_corner,z_corner, &
           br_corner,bz_corner,bt_corner, &
           pn_corner,ptpr_corner,ptpp_corner, &
           profn_travis_g,profn_travis_h,profn_travis_p,profn_travis_q, &
           profn_travis_w,proft_travis_g,proft_travis_h,proft_travis_p, &
           proft_travis_q,proft_travis_w, &
           MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
           RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
           MODEFR,MODEFW,IDEBUG,mdlplw, &
           MODELP,MODELV,NCMIN,NCMAX, &
           RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,RK0,RKANG0, &
           MODEL_ES,EPSRT,LMAXRT, &
           NS_NSA_DP,PMAX,EMAX,RHON_MIN,RHON_MAX, &
           NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_DP, &
           RF1,RFI1,RKX1,RKY1,RKZ1,RX1,RY1,RZ1,RK1, &
           RF2,RFI2,RKX2,RKY2,RKZ2,RX2,RY2,RZ2,RK2, &
           NGXMAX,NGYMAX,NGPMAX, &
           EPSDP,EPSRF,NORMF,NORMK,NFLOUT

    READ(NID,DP,IOSTAT=IST,ERR=9800,END=9900)
    IF(MODEL_PROF.EQ.0) THEN
       DO NS=1,NSMAX
          PROFN1(NS)=PROFN1(1)
          PROFN2(NS)=PROFN2(1)
          PROFT1(NS)=PROFT1(1)
          PROFT2(NS)=PROFT2(1)
          PROFU1(NS)=PROFU1(1)
          PROFU2(NS)=PROFU2(1)
       END DO
    END IF
    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE DPNLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE DPPLST

    WRITE(6,601)
    RETURN

  601 FORMAT(' ','# &DP : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOI,ZCOIL,BCOIL,NCOILMAX,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,'/ &
             9X,'PU,PUS,PUPR,PUPP,PNUC,PZCL,'/ &
             9X,'ID_NS,KID_NS,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'r_corner,z_corner,br_corner,bz_corner,bt_corner,'/ &
             9X,'pn_corner,ptpr_corner,ptpp_corner,'/ &
             9X,'profn_travis_g,profn_travis_h,profn_travis_p,'/ &
             9X,'profn_travis_q,profn_travis_w,proft_travis_g,'/ &
             9X,'proft_travis_h,proft_travis_p,proft_travis_q,'/ &
             9X,'proft_travis_w,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,BAXIS_SCALED,'/ &
             9X,'MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG,mdlplw'/ &
             9X,'MODELP,MODELV,NCMIN,NCMAX,'/ &
             9X,'RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,RK0,RKANG0,'/ &
             9X,'MODEL_ES,EPSRT,LMAXRT,'/ &
             9X,'NS_NSA_DP,PMAX,EMAX,ROHN_MIN,ROHN_MAX,'/ &
             9X,'NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_DP,'/ &
             9X,'RF1,RFI1,RKX1,RKY1,RKZ1,RX1,RY1,RZ1,RK1,'/ &
             9X,'RF2,RFI2,RKX2,RKY2,RKZ2,RX2,RY2,RZ2,RK2,'/ &
             9X,'NGXMAX,NGYMAX,NGPMAX,'/ &
             9X,'EPSDP,EPSRF,NORMF,NORMK,NFLOUT')
  END SUBROUTINE DPPLST

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE DP_CHEK(IERR)

    USE dpcomm_parm_local
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR

    IERR=0

    RETURN
  END SUBROUTINE DP_CHEK

!     ***** Setup velocity distribution function *****

  SUBROUTINE DPPREP(NTHMAX_DP_1,NRMAX_DP_1,RMIN_1,RMAX_1,IERR)

    USE dpcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTHMAX_DP_1,NRMAX_DP_1
    REAL(rkind),INTENT(IN):: RMIN_1,RMAX_1
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS

    NTHMAX_DP=NTHMAX_DP_1
    NRMAX_DP=NRMAX_DP_1
    DO NS=1,NSMAX
       RHON_MIN(NS)=RMIN_1
       RHON_MAX(NS)=RMAX_1
    END DO
    CALL DPPREP_LOCAL(IERR)
    RETURN
  END SUBROUTINE DPPREP

  SUBROUTINE DPPREP_LOCAL(IERR)

    USE dpcomm
    USE dpfpin
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,IND

    IERR=0
    IND=0
    DO NS=1,NSMAX
       IF(MODELV(NS).EQ.2.OR. &
          MODELV(NS).EQ.4.OR. &
          MODELV(NS).EQ.9) IND=1
    END DO

    IF(IND.EQ.1) THEN
       CALL DPLDFP(IERR)
       IF(IERR.NE.0) RETURN
    END IF

    DO NS=1,NSMAX
       IF(MODELV(NS).EQ.1) THEN
          CALL DPLDFM(NS,0,IERR)
          IF(IERR.NE.0) RETURN
       END IF
       IF(MODELV(NS).EQ.3) THEN
          CALL DPLDFM(NS,1,IERR)
          IF(IERR.NE.0) RETURN
       END IF
    END DO
    RETURN
  END SUBROUTINE DPPREP_LOCAL

!     ****** SHOW PARAMETERS ******

  SUBROUTINE DP_VIEW

    USE dpcomm_parm_local
    IMPLICIT NONE
    INTEGER:: NS

    WRITE(6,100)
    DO NS=1,NSMAX
       WRITE(6,110) &
            NS,MODELP(NS),MODELV(NS),NCMIN(NS),NCMAX(NS)
    ENDDO
    WRITE(6,101)
    DO NS=1,NSMAX
       WRITE(6,111) &
            NS,PMAX(NS),EMAX(NS),RHON_MIN(NS),RHON_MAX(NS)
    ENDDO
    WRITE(6,612) ' NPMAX_DP   ',NPMAX_DP ,'NTHMAX_DP  ',NTHMAX_DP
    WRITE(6,612) ' NRMAX_DP   ',NRMAX_DP ,'NSAMAX_DP  ',NSAMAX_DP
    WRITE(6,602) ' NORMF ',NORMF ,'NORMK ',NORMK
    WRITE(6,612) ' MODEL_ES   ',MODEL_ES
    WRITE(6,601) ' EPSRT ',EPSRT ,'EPSDP ',EPSDP
    WRITE(6,601) ' EPSRF ',EPSRF
    WRITE(6,602) ' LMAXRT',LMAXRT,'NFLOUT',NFLOUT

    RETURN

100 FORMAT('NS    MODELP  MODELV  NCMIN   NCMAX')
101 FORMAT('NS    PMAX        EMAX        RHON_MIN    RHON_MAX')
110 FORMAT(I2,' ',4I8)                               
111 FORMAT(I2,' ',1P4E12.4)                               
601 FORMAT(A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
602 FORMAT(A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
611 FORMAT(A12,'=',1PE11.3 :2X,A12,'=',1PE11.3: &
            2X,A12,'=',1PE11.3)
612 FORMAT(A12,'=',I7,4X   :2X,A12,'=',I7,4X: &
            2X,A12,'=',I7)
  END SUBROUTINE DP_VIEW

  ! --- Broadcast dp input parameters
  
  SUBROUTINE dp_broadcast

    USE dpcomm_parm
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata

    idata(1)=NPMAX_DP
    idata(2)=NTHMAX_DP
    idata(3)=NRMAX_DP
    idata(4)=NSAMAX_DP
    idata(5)=LMAXRT
    idata(6)=MODEL_ES

    CALL mtx_broadcast_integer(idata,6)

    NPMAX_DP=idata(1)
    NTHMAX_DP=idata(2)
    NRMAX_DP=idata(3)
    NSAMAX_DP=idata(4)
    LMAXRT=idata(5)
    MODEL_ES=idata(6)

    rdata(1)=EPSRT

    CALL mtx_broadcast_real8(rdata,1)

    EPSRT=rdata(1)
    
    CALL mtx_broadcast_integer(NS_NSA_DP,NSAMAX_DP)
    CALL mtx_broadcast_integer(MODELP,NSAMAX_DP)
    CALL mtx_broadcast_integer(MODELV,NSAMAX_DP)
    CALL mtx_broadcast_integer(NCMIN,NSAMAX_DP)
    CALL mtx_broadcast_integer(NCMAX,NSAMAX_DP)

    CALL mtx_broadcast_real8(PMAX,NSAMAX_DP)
    CALL mtx_broadcast_real8(EMAX,NSAMAX_DP)
    CALL mtx_broadcast_real8(RHON_MIN,NSAMAX_DP)
    CALL mtx_broadcast_real8(RHON_MAX,NSAMAX_DP)

    RETURN
  END SUBROUTINE dp_broadcast

  ! --- Broadcast local dp input parameters
  
  SUBROUTINE dp_broadcas_localt

    USE dpcomm_parm_local
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata

    idata(1)=NGXMAX
    idata(2)=NGYMAX
    idata(3)=NGPMAX
    idata(4)=NORMF
    idata(5)=NORMK
    idata(6)=NFLOUT

    CALL mtx_broadcast_integer(idata,6)

    NGXMAX=idata(1)
    NGYMAX=idata(2)
    NGPMAX=idata(3)
    NORMF=idata(4)
    NORMK=idata(5)
    NFLOUT=idata(6)

    rdata( 1)=RF0
    rdata( 2)=RFI0
    rdata( 3)=RKX0
    rdata( 4)=RKY0
    rdata( 5)=RKZ0
    rdata( 6)=RX0
    rdata( 7)=RY0
    rdata( 8)=RZ0
    rdata( 9)=RK0
    rdata(10)=RKANG0
    rdata(11)=RF1
    rdata(12)=RFI1
    rdata(13)=RKX1
    rdata(14)=RKY1
    rdata(15)=RKZ1
    rdata(16)=RX1
    rdata(17)=RY1
    rdata(18)=RZ1
    rdata(19)=RK1
    rdata(20)=RF2
    rdata(21)=RFI2
    rdata(22)=RKX2
    rdata(23)=RKY2
    rdata(24)=RKZ2
    rdata(25)=RX2
    rdata(26)=RY2
    rdata(27)=RZ2
    rdata(28)=RK2
    rdata(29)=EPSDP
    rdata(30)=EPSRF

    CALL mtx_broadcast_real8(rdata,30)

    RF0=rdata( 1)
    RFI0=rdata( 2)
    RKX0=rdata( 3)
    RKY0=rdata( 4)
    RKZ0=rdata( 5)
    RX0=rdata( 6)
    RY0=rdata( 7)
    RZ0=rdata( 8)
    RK0=rdata( 9)
    RKANG0=rdata(10)
    RF1=rdata(11)
    RFI1=rdata(12)
    RKX1=rdata(13)
    RKY1=rdata(14)
    RKZ1=rdata(15)
    RX1=rdata(16)
    RY1=rdata(17)
    RZ1=rdata(18)
    RK1=rdata(19)
    RF2=rdata(20)
    RFI2=rdata(21)
    RKX2=rdata(22)
    RKY2=rdata(23)
    RKZ2=rdata(24)
    RX2=rdata(25)
    RY2=rdata(26)
    RZ2=rdata(27)
    RK2=rdata(28)
    EPSDP=rdata(29)
    EPSRF=rdata(30)

    idata(1)=MODEL_ES
    idata(2)=LMAXRT
    idata(3)=NPMAX_DP
    idata(4)=NTHMAX_DP
    idata(5)=NRMAX_DP
    idata(6)=NSAMAX_DP
    CALL mtx_broadcast_integer(idata,6)
    MODEL_ES=idata(1)
    LMAXRT=idata(2)
    NPMAX_DP=idata(3)
    NTHMAX_DP=idata(4)
    NRMAX_DP=idata(5)
    NSAMAX_DP=idata(6)

    rdata(1)=EPSRT
    CALL mtx_broadcast_real8(rdata,1)
    EPSRT=rdata(1)

    RETURN
  END SUBROUTINE dp_broadcas_localt
END MODULE DPPARM
