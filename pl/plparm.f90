! plparm.f90

MODULE plparm

  PRIVATE
  PUBLIC pl_parm,pl_view,pl_broadcast

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
      IMPLICIT NONE
      INTEGER,INTENT(IN):: mode
      CHARACTER(LEN=*),INTENT(IN)::  kin
      INTEGER,INTENT(OUT):: ierr

    1 CALL TASK_PARM(MODE,'PL',KIN,plnlin,plplst,IERR)
      IF(IERR.NE.0) RETURN

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
           MODEFR,MODEFW,IDEBUG,mdlplw

      READ(NID,PL,IOSTAT=IST,ERR=9800,END=9900)

      IF(MODEL_PROF.EQ.0) THEN
         DO NS=2,NSMAX
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
    END SUBROUTINE plnlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE plplst

      implicit none
      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &PL : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
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
             9X 'proft_travis_w,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,BAXIS_SCALED,'/ &
             9X,'MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG,mdlplw')
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
         IF((MODELN.LT.0).OR.(MODELN.GT.31)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ELSE
         IF((MODELN.LT.0).OR.(MODELN.GT.31)) THEN
            WRITE(6,*) 'XX plcheck: INVALID MODELN: MODELN=',MODELN
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
    idata( 4)=MODELN
    idata( 5)=MODELQ
    idata( 6)=IDEBUG
    idata( 7)=MODEFR
    idata( 8)=MODEFW
    idata( 9)=mdlplw
    idata(10)=MODEL_PROF
    idata(11)=MODEL_NPROF
    idata(12)=NCOILMAX

    CALL mtx_broadcast_integer(idata,12)
    
    NSMAX=idata( 1)
    MODELG=idata( 2)
    MODELB=idata( 3)
    MODELN=idata( 4)
    MODELQ=idata( 5)
    IDEBUG=idata( 6)
    MODEFR=idata( 7)
    MODEFW=idata( 8)
    mdlplw=idata( 9)
    MODEL_PROF=idata(10)
    MODEL_NPROF=idata(11)
    NCOILMAX=idata(12)

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
    CALL mtx_broadcast_real8(PROFT1,NSMAX)
    CALL mtx_broadcast_real8(PROFT2,NSMAX)
    CALL mtx_broadcast_real8(PROFU1,NSMAX)
    CALL mtx_broadcast_real8(PROFU2,NSMAX)
    CALL mtx_broadcast_real8(PNUC,NSMAX)
    CALL mtx_broadcast_real8(PZCL,NSMAX)
    CALL mtx_broadcast_integer(NPA,NSMAX)
    CALL mtx_broadcast_integer(ID_NS,NSMAX)
    DO NS=1,NSMAX
       CALL mtx_broadcast_character(KID_NS(NS),2)
    END DO

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
  END SUBROUTINE pl_broadcast

!     ****** SHOW PARAMETERS ******

    SUBROUTINE pl_view

      use plcomm_parm
      implicit none
      integer:: NS,i,NCOIL

      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    , &
                   'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'RKAP  ',RKAP  ,'RDLT  ',RDLT  , &
                   'Q0    ',Q0    ,'QA    ',QA
      WRITE(6,601) 'RIP   ',RIP   ,'PROFJ ',PROFJ

      IF(MODELG.EQ.0) THEN
         WRITE(6,601) 'RMIR  ',RMIR  ,'ZBB   ',ZBB
      ENDIF
      IF(MODELG.EQ.11) THEN
         WRITE(6,601) 'Hpitc1',Hpitch1,'Hpitc2',Hpitch2, &
                      'RRCH  ',RRCH
      END IF
      IF(MODELG.EQ.0.AND.MOD(MODELB,2).EQ.1) THEN
         WRITE(6,'(A)') '     NCOIL   RCOIL       ZCOIL       BCOIL'
         DO NCOIL=1,NCOILMAX
            WRITE(6,'(I10,1P3E12.4)') &
                 NCOIL,RCOIL(NCOIL),ZCOIL(NCOIL),BCOIL(NCOIL)
         END DO
      END IF

      WRITE(6,601) 'RHOEDG',RHOEDG,'RHOGMN',RHOGMN, &
                   'RHOGMX',RHOGMX
      WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB, &
                   'MODELN',MODELN,'MODELQ',MODELQ
      WRITE(6,604) 'MODEFR',MODEFR,'MODEFW',MODEFW, &
                   'mdlplw',mdlplw
      WRITE(6,'(A,I5)') ' MODEL_PROF  =',MODEL_PROF
      WRITE(6,'(A,I5)') ' MODEL_NPROF =',MODEL_NPROF

      WRITE(6,100)
      DO NS=1,NSMAX
         WRITE(6,110) NS,NPA(NS),PA(NS),PZ(NS),PN(NS),PNS(NS)
      ENDDO
      WRITE(6,120)
      DO NS=1,NSMAX
         WRITE(6,130) NS,PTPR(NS),PTPP(NS),PTS(NS),PU(NS),PUS(NS)
      ENDDO
      WRITE(6,131)
      DO NS=1,NSMAX
         WRITE(6,132) NS,PUPR(NS),PUPP(NS),PNUC(NS),PZCL(NS)
      ENDDO
      WRITE(6,140)
      DO NS=1,NSMAX
         WRITE(6,150) NS,RHOITB(NS),PNITB(NS),PTITB(NS),PUITB(NS)
      END DO

      WRITE(6,160)
      DO NS=1,NSMAX
         WRITE(6,170) NS,PROFN1(NS),PROFN2(NS),PROFT1(NS),PROFT2(NS), &
                         PROFU1(NS),PROFU2(NS)
      END DO

      WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  , &
                   'RF_PL ',RF_PL

      IF(MODELG.EQ.11.OR.MODELG.EQ.13) THEN
         WRITE(6,606) 'r_corner:  ',(r_corner(i),i=1,3)
         WRITE(6,606) 'z_corner:  ',(z_corner(i),i=1,3)
         WRITE(6,606) 'br_corner: ',(br_corner(i),i=1,3)
         WRITE(6,606) 'bz_corner: ',(bz_corner(i),i=1,3)
         WRITE(6,606) 'bt_corner: ',(bt_corner(i),i=1,3)
         DO ns=1,NSMAX
            WRITE(6,'(A,I3)') 'ns=',ns
            WRITE(6,606) 'pn_corner: ',(pn_corner(i,ns),i=1,3)
            WRITE(6,606) 'ptpr_corner:',(ptpr_corner(i,ns),i=1,3)
            WRITE(6,606) 'ptpp_corner:',(ptpp_corner(i,ns),i=1,3)
         END DO
      END IF

      RETURN

  100 FORMAT(' ','NS    NPA         PA          PZ          ', &
                       'PN          PNS')
  110 FORMAT(' ',I2,' ',I5,7X,1P4E12.4)
  120 FORMAT(' ','NS    PTPR        PTPP        PTS         ', &
                       'PU          PUS')
  130 FORMAT(' ',I2,' ',1P5E12.4)                               
  131 FORMAT(' ','NS    PUPR        PUPP        PNUC        PZCL')
  132 FORMAT(' ',I2,' ',1P4E12.4)                               
  140 FORMAT(' ','NS    RHOITB      PNITB       PTITB       PUITB')
  150 FORMAT(' ',I2,' ',1P4E12.4)                               
  160 FORMAT(' ','NS    PROFN1      PROFN2      PROFT1      ', &
                       'PROFT2      PROFU1      PROFU2')
  170 FORMAT(' ',I2,' ',1P6E12.4)                               
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
             2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
             2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  606 FORMAT(' ',A12,1P3E12.4)
    END SUBROUTINE pl_view
 END MODULE plparm
