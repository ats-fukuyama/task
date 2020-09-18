! wmgsub.f90

MODULE wmgsub

  PRIVATE
  PUBLIC wm_gsub,wm_greq,wm_greqg,wm_gfwr,wm_gxeq,wm_gxeqp,wm_g3dr,wm_gprm, &
       wm_circle

CONTAINS

! ****** DRAW LINES OF 1D GRAPH ******

  SUBROUTINE wm_gsub(NX,GX,GXMIN,GXMAX,GY,NY,GP,KTITL)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NX,NY
    REAL(4),INTENT(IN):: GX(NX),GXMIN,GXMAX,GY(NX,NY),GP(4)
    CHARACTER(LEN=6):: KTITL*6
    INTEGER:: I
    REAL(4):: GYMIN,GYMAX,GYMIN1,GYMAX1,GSYMIN,GSYMAX,GSX,GSY
    INTEGER:: ICL(6),IPAT(6)
    DATA ICL/7,6,5,4,3,2/,IPAT/0,2,4,6,3,1/

    IF(NY.EQ.0) RETURN

    GYMIN= 1.E32
    GYMAX=-1.E32
    DO I=1,NY
       CALL GMNMX1(GY(1,I),1,NX,1,GYMIN1,GYMAX1)
       GYMIN=MIN(GYMIN,GYMIN1)
       GYMAX=MAX(GYMAX,GYMAX1)
    ENDDO

!     CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
    GSX=0.1

    IF(GYMIN.GT.0.0.AND.GYMAX.GT.0.0) THEN
       GYMIN=0.0
    ELSE IF(GYMIN.LT.0.0.AND.GYMAX.LT.0.0) THEN
       GYMAX=0.0
    ENDIF
    CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)

    IF(ABS(GSYMAX-GSYMIN).LT.1.E-15) GOTO 9000

    CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),GXMIN,GXMAX,GSYMIN,GSYMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
    CALL GVALUE(GXMIN,GSX*5,0.0,0.0,NGULEN(GSX*5))
    CALL GSCALE(0.0,0.0,0.0,GSY,0.1,9)
    CALL GVALUE(0.0,0.0,0.0,GSY*2,NGULEN(GSY*2))

    DO I=1,NY
       CALL SETLIN(0,2,ICL(I))
       CALL GPLOTP(GX,GY(1,I),1,NX,1,0,0,IPAT(I))
    ENDDO
    CALL SETLIN(0,2,7)

9000 CONTINUE
    CALL MOVE(GP(1),GP(4)+0.2)
    CALL TEXT(KTITL,6)

    RETURN
  END SUBROUTINE wm_gsub

!     ****** DRAW 2D POLOIDAL GRAPH ******

  SUBROUTINE wm_greq(K2,K3,K4)

    USE wmcomm
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2,K3,K4
    REAL(4),ALLOCATABLE:: GY(:,:)
    COMPLEX(rkind):: CFL
    INTEGER:: NA3,NG3,NHH,NTH,NR

    ALLOCATE(GY(nrmax+1,nthmax))

    IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
       IF(K3.EQ.'R') THEN
          NA3=1
          NG3=1
       ELSEIF(K3.EQ.'T') THEN
          NA3=1
          NG3=2
       ELSEIF(K3.EQ.'Z') THEN
          NA3=1
          NG3=3
       ELSEIF(K3.EQ.'S') THEN
          NA3=2
          NG3=1
       ELSEIF(K3.EQ.'H') THEN
          NA3=2
          NG3=2
       ELSEIF(K3.EQ.'B') THEN
          NA3=2
          NG3=3
       ELSEIF(K3.EQ.'+') THEN
          NA3=3
          NG3=1
       ELSEIF(K3.EQ.'-') THEN
          NA3=3
          NG3=2
       ELSEIF(K3.EQ.'P') THEN
          NA3=3
          NG3=3
       ELSE
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN wm_greq'
          GOTO 9000
       ENDIF
    ELSEIF(K2.EQ.'P') THEN
       IF(K3.EQ.'E') THEN
          NG3=1
       ELSEIF(K3.EQ.'D') THEN
          NG3=2
       ELSEIF(K3.EQ.'T') THEN
          NG3=3
       ELSEIF(K3.EQ.'1') THEN
          NG3=1
       ELSEIF(K3.EQ.'2') THEN
          NG3=2
       ELSEIF(K3.EQ.'3') THEN
          NG3=3
       ELSEIF(K3.EQ.'4') THEN
          NG3=4
       ELSEIF(K3.EQ.'5') THEN
          NG3=5
       ELSEIF(K3.EQ.'6') THEN
          NG3=6
       ELSE
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN wm_greq'
          GOTO 9000
       ENDIF
    ELSE IF(K2.EQ.'J') THEN
       NG3=0
    ELSE IF(K2.EQ.'F') THEN
       IF(K3.EQ.'S') THEN
          NG3=1
       ELSEIF(K3.EQ.'B') THEN
          NG3=2
       ELSEIF(K3.EQ.'Q') THEN
          NG3=3
       ELSEIF(K3.EQ.'2') THEN
          NG3=4
       ELSEIF(K3.EQ.'3') THEN
          NG3=5
       ELSEIF(K3.EQ.'J') THEN
          NG3=6
       ELSEIF(K3.EQ.'R') THEN
          NG3=7
       ELSEIF(K3.EQ.'Z') THEN
          NG3=8
       ELSEIF(K3.EQ.'4') THEN
          NG3=9
       ELSEIF(K3.EQ.'5') THEN
          NG3=10
       ELSEIF(K3.EQ.'6') THEN
          NG3=11
       ELSEIF(K3.EQ.'7') THEN
          NG3=12
       ELSE
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN wm_greq'
          GOTO 9000
       ENDIF
    ELSE
       WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #2 IN wm_greq'
       GOTO 9000
    ENDIF

2   IF(NHHMAX.EQ.1) THEN
       NHH=1
    ELSE
       WRITE(6,*) '## INPUT NHH : 1..',NHHMAX
       READ(5,*,ERR=2,END=9000) NHH
    END IF
    IF(NHH.LT.1.OR.NHH.GT.NHHMAX) THEN
       WRITE(6,*) 'XX ILLEGAL NHH'
       GOTO 2
    ENDIF

    DO NTH=1,MDSIZ
       DO NR=1,NRMAX+1
          IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
             IF(K2.EQ.'E') THEN
                IF(NA3.EQ.1) THEN
                   CFL=CEFLD(NG3,NTH,NHH,NR)
                ELSEIF(NA3.EQ.2) THEN
                   CFL=CEN(NG3,NTH,NHH,NR)
                ELSEIF(NA3.EQ.3) THEN
                   CFL=CEP(NG3,NTH,NHH,NR)
                ENDIF
             ELSEIF(K2.EQ.'B') THEN
                CFL=CBFLD(NG3,NTH,NHH,NR)
             ENDIF

             IF(K4.EQ.'R') THEN
                GY(NR,NTH)=GUCLIP(DBLE(CFL))
             ELSEIF(K4.EQ.'I') THEN
                GY(NR,NTH)=GUCLIP(DIMAG(CFL))
             ELSEIF(K4.EQ.'A') THEN
                GY(NR,NTH)=GUCLIP(ABS(CFL))
             ELSE
                WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN wm_grth'
                GOTO 9000
             ENDIF
          ELSEIF(K2.EQ.'P') THEN
             IF(NR.LE.NRMAX) THEN
                GY(NR,NTH)=GUCLIP(PABS(NTH,NHH,NR,NG3))
             ELSE
                GY(NR,NTH)=0.0
             ENDIF
          ELSE IF (K2.EQ.'J') THEN
             IF(NR.LE.NRMAX) THEN
                GY(NR,NTH)=GUCLIP(PCUR(NTH,NHH,NR))
             ELSE
                GY(NR,NTH)=0.0
             ENDIF
          ELSE IF (K2.EQ.'F') THEN
             IF(NG3.EQ.1) THEN
                GY(NR,NTH) =GUCLIP(PSIP(NR))
             ELSEIF(NG3.EQ.2) THEN
                GY(NR,NTH) =GUCLIP(BPST(NTH,NHH,NR))
             ELSEIF(NG3.EQ.3) THEN
                GY(NR,NTH) =GUCLIP(BFLD(3,NTH,NHH,NR)/BFLD(2,NTH,NHH,NR))
             ELSEIF(NG3.EQ.4) THEN
                GY(NR,NTH) =GUCLIP(BFLD(2,NTH,NHH,NR))
             ELSEIF(NG3.EQ.5) THEN
                GY(NR,NTH) =GUCLIP(BFLD(3,NTH,NHH,NR))
             ELSEIF(NG3.EQ.6) THEN
                GY(NR,NTH) =GUCLIP(RJ(NTH,NHH,NR))
             ELSEIF(NG3.EQ.7) THEN
                GY(NR,NTH) =GUCLIP(RPS(NTH,NR))
             ELSEIF(NG3.EQ.8) THEN
                GY(NR,NTH) =GUCLIP(ZPS(NTH,NR))
             ELSEIF(NG3.EQ.9) THEN
                GY(NR,NTH) =GUCLIP(BTP(NTH,NR))
             ELSEIF(NG3.EQ.10) THEN
                GY(NR,NTH) =GUCLIP(BPT(NTH,NR))
             ELSEIF(NG3.EQ.11) THEN
                GY(NR,NTH) =GUCLIP(BPR(NTH,NR))
             ELSEIF(NG3.EQ.12) THEN
                GY(NR,NTH) =GUCLIP(BPZ(NTH,NR))
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    IF(NGRAPH.EQ.0) THEN
       CALL wm_gfwr(GY,NHH,K2,K3,K4)
       GOTO 9000
    ELSEIF(NGRAPH.EQ.1) THEN
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL wm_gxeq(GY,NHH,K2,K3)
    ELSE IF(NGRAPH.EQ.2) THEN
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL wm_gxeqp(GY,NHH,K2,K3)
    ELSE IF(NGRAPH.EQ.3) THEN
       IF(NTHMAX.LE.2) THEN
          WRITE(6,*) 'XX NTHMAX:',NTHMAX,'  CONDITION:NTHMAX >= 4'
          GO TO 9000
       ENDIF
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL wm_g3da(GY,NHH,K2,K3,0)
    ELSE IF(NGRAPH.EQ.4) THEN
       IF(NTHMAX.LE.2) THEN
          WRITE(6,*) 'XX NTHMAX:',NTHMAX,'  CONDITION:NTHMAX >= 4'
          GO TO 9000
       ENDIF
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL wm_g3da(GY,NHH,K2,K3,4)
    ELSE
       WRITE(6,*) 'XX wm_greq: UNDEFINED NGRAPH: NGRAPH=',NGRAPH
       GO TO 9000
    ENDIF

    CALL MOVE(20.0,17.5)
    IF(K2.EQ.'E') CALL TEXT('E ',2)
    IF(K2.EQ.'B') CALL TEXT('B ',2)
    IF(K2.EQ.'P') CALL TEXT('P ',2)
    IF(K2.EQ.'J') CALL TEXT('J ',2)
    CALL MOVE(20.0,17.1)
    IF(K2.NE.'P') THEN
       IF(K3.EQ.'R') CALL TEXT('r ',2)
       IF(K3.EQ.'T') CALL TEXT('theta ',6)
       IF(K3.EQ.'Z') CALL TEXT('phi ',4)
       IF(K3.EQ.'S') CALL TEXT('s ',2)
       IF(K3.EQ.'H') CALL TEXT('h ',2)
       IF(K3.EQ.'B') CALL TEXT('b ',2)
       IF(K3.EQ.'+') CALL TEXT('+ ',2)
       IF(K3.EQ.'-') CALL TEXT('- ',2)
       IF(K3.EQ.'P') CALL TEXT('P ',2)
    ELSE
       IF(K3.EQ.'E') CALL TEXT('e ',2)
       IF(K3.EQ.'D') CALL TEXT('D ',2)
       IF(K3.EQ.'T') CALL TEXT('T ',2)
       IF(K3.EQ.'1') CALL TEXT('1 ',2)
       IF(K3.EQ.'2') CALL TEXT('2 ',2)
       IF(K3.EQ.'3') CALL TEXT('3 ',2)
       IF(K3.EQ.'4') CALL TEXT('4 ',2)
       IF(K3.EQ.'5') CALL TEXT('5 ',2)
       IF(K3.EQ.'6') CALL TEXT('6 ',2)
    ENDIF
    CALL MOVE(20.0,16.7)
    IF(K4.EQ.'R') CALL TEXT('Real ',5)
    IF(K4.EQ.'I') CALL TEXT('Imag ',5)
    IF(K4.EQ.'A') CALL TEXT('Abs ',4)

    CALL wm_gprm('C',K3,0,0,0,0)

    CALL PAGEE

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_greq

!     ****** DRAW 1D RADIAL GRAPHS ******

  SUBROUTINE wm_greqg(K2,K3,K4)

    USE wmcomm
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2,K3,K4
    REAL(4),ALLOCATABLE:: GY(:,:)
    CHARACTER(LEN=6)::  RTITL(8)
    INTEGER:: NHH,NTH,NR
    REAL(4):: GP(4,4)
    DATA GP/ 3.0, 10.8,  9.5, 16.5, &
             3.0, 10.8,  1.0,  8.0, &
            13.8, 21.6,  9.5, 16.5, &
            13.8, 21.6,  1.0,  8.0/

    ALLOCATE(GY(nrmax+1,nthmax))

2   CONTINUE
    IF(NHHMAX.EQ.1) THEN
       NHH=1
    ELSE
       WRITE(6,*) '## INPUT NHH : 1..',NHHMAX
       READ(5,*,ERR=2,END=9000) NHH
    ENDIF
    IF(NHH.LT.1.OR.NHH.GT.NHHMAX) THEN
       WRITE(6,*) 'XX ILLEGAL NHH'
       GOTO 2
    ENDIF

    CALL PAGES
    CALL SETCHS(0.3,0.0)
    RTITL(1)='RG11  '
    RTITL(2)='RG12  '
    RTITL(3)='RG13  '
    RTITL(4)='RG22  '
    RTITL(5)='RG23  '
    RTITL(6)='RG33  '
    RTITL(7)='RJ    '
    RTITL(8)='PSI   '
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG11(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(1),GP(1,1))
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG12(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(2),GP(1,2))
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG13(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(3),GP(1,3))
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG22(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(4),GP(1,4))
    CALL wm_gprm('C',K3,0,0,0,0)
    CALL PAGEE

    CALL PAGES
    CALL SETCHS(0.3,0.0)
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG23(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(5),GP(1,1))
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RG33(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(6),GP(1,2))
    DO NTH=1,NTHMAX
       DO NR=1,NRMAX+1
          GY(NR,NTH)=GUCLIP(RJ(NTH,NHH,NR))
       ENDDO
    ENDDO
    CALL wm_ggr(GY,NTHMAX,RTITL(7),GP(1,3))
    DO NR=1,NRMAX+1
       GY(NR,1)=GUCLIP(PSIP(NR))
    ENDDO
    CALL wm_ggr(GY,1,RTITL(8),GP(1,4))
    CALL wm_gprm('C',K3,0,0,0,0)
    CALL PAGEE     
    CALL MOVE(20.0,17.5)
    IF(K2.EQ.'G') CALL TEXT('RG',2)

    CALL wm_gprm('C',K3,0,0,0,0)

    CALL PAGEE

9000 RETURN
  END SUBROUTINE wm_greqg

!     ****** WRITE GRAPHIC DATA IN FILE ******

  SUBROUTINE wm_gfwr(GGL,NHH,K2,K3,K4)

    USE wmcomm
    IMPLICIT NONE
    REAL(4),INTENT(IN):: GGL(nrmax,nthmax)
    INTEGER,INTENT(IN):: NHH
    CHARACTER(LEN=1),INTENT(IN):: K2,K3,K4
    REAL(4),ALLOCATABLE:: GRL(:,:),GZL(:,:)
    INTEGER:: NR,NTH,NFD
    
    ALLOCATE(GRL(nrmax,nthmax),GZL(nrmax,nthmax))

    DO NR=1,NRMAX+1
       DO NTH=1,NTHMAX
          GRL(NR,NTH)=GUCLIP(RPS(NTH,NR))
          GZL(NR,NTH)=GUCLIP(ZPS(NTH,NR))
       ENDDO
    ENDDO

    NFD=23
    WRITE(NFD,'(3A1,I8)') K2,K3,K4,NHH
    WRITE(NFD,'(2I8)') NRMAX+1,NTHMAX
    WRITE(NFD,'(1P3E15.7)') ((GRL(NR,NTH),GZL(NR,NTH),GGL(NR,NTH), &
                              NR=1,NRMAX+1),NTH=1,NTHMAX)

    WRITE(6,*) '## wm_gfwr: Data written in fort.23'
    RETURN
  END SUBROUTINE wm_gfwr

!     ****** DRAW COUNTOUR IN MAGNETIC SURFACE COORDINATES ******

  SUBROUTINE wm_gxeq(GGL,NHH,K2,K3)

    USE wmcomm
    USE wmprof
    IMPLICIT NONE
    REAL(4),INTENT(IN):: GGL(nrmax+1,nthmax)
    INTEGER,INTENT(IN):: NHH
    CHARACTER(LEN=1):: K2,K3
    REAL(4),ALLOCATABLE:: GBY(:,:),GFL(:,:),GRL(:,:),GZL(:,:),GRS(:),GZS(:)
    REAL(4),ALLOCATABLE:: GTHR(:,:),GTCO(:,:),GTRC(:,:),GTLC(:,:)
    REAL(rkind),ALLOCATABLE:: THR(:,:),TCO(:,:),TRC(:,:),TLC(:,:)
    REAL(rkind):: RN(nsmax),RTPR(nsmax),RTPP(nsmax),RU(nsmax)
    REAL(rkind):: WP(nsmax),WC(nsmax)
    INTEGER:: NS,NR,NTH,NTHF,NHHF,NTHGS,NTHP,NTHL,IRORG,NC,NRLCFS,NSTEP
    INTEGER:: NTHG,NSU
    REAL(rkind):: RCOS,RSIN,FACT,VAL,WF,RHOL,DTH,DPH
    REAL(rkind):: BYL,BST,BSP,RKTH,RKPH,RKPR,RNPR
    REAL(rkind):: RNL,AM,AE,FACTORHR,FACTORRC,FACTORLC
    REAL(rkind):: BABS,BSUPTH,BSUPPH,BABSP,BSUPTHP,BSUPPHP
    REAL(rkind):: BCF,DTHG
    REAL(4):: GRLEN,GZLEN,GPR,GPZ,GRORG
    REAL(4):: GGRMIN,GGRMAX,GGRSTP,GGZMIN,GGZMAX,GGZSTP
    REAL(4):: GBCF,GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP

    ALLOCATE(THR(nrmax+1,nthmax_f),TCO(nrmax+1,nthmax_f))
    ALLOCATE(TRC(nrmax+1,nthmax_f),TLC(nrmax+1,nthmax_f))
    ALLOCATE(GBY(nrmax+1,nthmax_f))
    ALLOCATE(GFL(nrmax+1,nthmax_g))
    ALLOCATE(GRL(nrmax+1,nthmax_g),GZL(nrmax+1,nthmax_g))
    ALLOCATE(GTHR(nrmax+1,nthmax_f),GTCO(nrmax+1,nthmax_f))
    ALLOCATE(GTRC(nrmax+1,nthmax_f),GTLC(nrmax+1,nthmax_f))
    ALLOCATE(GRS(nsumax+1),GZS(nsumax+1))

    SELECT CASE(MODELG)
    CASE(0:2)
       DTHG=2.D0*PI/nthmax_g
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             RCOS=COS(DTHG*(NTH-1))
             RSIN=SIN(DTHG*(NTH-1))
             IF(MODELG.EQ.0) THEN
                RPSG(NTH,NR)  =      XR(NR)*RCOS
             ELSE
                RPSG(NTH,NR)  = RR + XR(NR)*RCOS
             ENDIF
             ZPSG(NTH,NR)  =         XR(NR)*RSIN
          ENDDO
       END DO
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             GRL(NR,NTH)=GUCLIP(RPSG(NTH,NR))
             GZL(NR,NTH)=GUCLIP(ZPSG(NTH,NR))
          ENDDO
       ENDDO
    CASE(3)
       ! RPSG and ZPSG loaded in wmsetg_tor
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             GRL(NR,NTH)=GUCLIP(RPSG(NTH,NR))
             GZL(NR,NTH)=GUCLIP(ZPSG(NTH,NR))
          ENDDO
       ENDDO
    CASE(4,6)
       nthmax_g=NTHMAX
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             NHHF=(NHH-1)*factor_nhh +1
             NTHF=(NTH-1)*factor_nth +1
             GRL(NR,NTH)=GUCLIP(RPST(NTHF,NHHF,NR))
             GZL(NR,NTH)=GUCLIP(ZPST(NTHF,NHHF,NR))
          ENDDO
       ENDDO
    END SELECT
    
    NTHGS=nthmax_g/NTHMAX
    DO NR=1,NRMAX+1
       DO NTH=1,NTHMAX
          NTHP=NTH+1
          IF(NTHP.GT.NTHMAX) NTHP=1
          DO NTHG=1,NTHGS
             NTHL=(NTH-1)*NTHGS+NTHG
             FACT=DBLE(NTHG-1)/DBLE(NTHGS)
             VAL=(1.D0-FACT)*DBLE(GGL(NR,NTH))+FACT*DBLE(GGL(NR,NTHP))
             GFL(NR,NTHL)=GUCLIP(VAL)
          ENDDO
       ENDDO
    ENDDO

!     ***** CALCULATE RESONANCE AND CUTOFF *****

    IF(K2.EQ.'P') THEN

       WF=2*PI*RF*1.0D6
       DO NR=1,NRMAX+1
          rhol=xrho(nr)
          CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
          dth=2.d0*pi/nthmax
          dph=2.d0*pi/nphmax
          DO NTH=1,NTHMAX
             NTHP=NTH+1
             IF(NTHP.GT.NTHMAX) NTHP=1
             CALL WMCMAG(NR,2*NTH-1, 2*NHH-1,BABS, BSUPTH, BSUPPH )
             CALL WMCMAG(NR,2*NTHP-1,2*NHH-1,BABSP,BSUPTHP,BSUPPHP)
             DO NTHG=1,NTHGS
                NTHL=(NTH-1)*NTHGS+NTHG
                FACT=DBLE(NTHG-1)/DBLE(NTHGS)
                BYL=(1.D0-FACT)*BABS  +FACT*BABSP
                BST=(1.D0-FACT)*BSUPTH+FACT*BSUPTHP
                BSP=(1.D0-FACT)*BSUPPH+FACT*BSUPPHP
                RKTH=NTH0*BST/BYL
                RKPH=NPH0*BSP/BYL
                RKPR=RKTH+RKPH
                RNPR=RKPR/(WF/VC)
                GBY(NR,NTHL)=GUCLIP(BYL)
                DO NS=1,NSMAX
                   RNL=RN(NS)*1.0D20
                   AM=PA(NS)*AMP
                   AE=PZ(NS)*AEE
                   WP(NS)=(AE**2)*RNL/(AM*EPS0)
                   WC(NS)=AE*BYL/AM
                ENDDO

                FACTORHR=1.D0
                FACTORRC=1.D0
                FACTORLC=1.D0
                DO NS=1,NSMAX
                   FACTORHR=FACTORHR*(1.D0-(WC(NS)/WF)**2)
                   FACTORRC=FACTORRC*(1.D0+ WC(NS)/WF)
                   FACTORLC=FACTORLC*(1.D0- WC(NS)/WF)
                ENDDO
                TCO(NR,NTHL)=1.D0
                THR(NR,NTHL)=FACTORHR
                TRC(NR,NTHL)=FACTORRC*(1.D0-RNPR**2)
                TLC(NR,NTHL)=FACTORLC*(1.D0-RNPR**2)
                DO NS=1,NSMAX
                   TCO(NR,NTHL)=TCO(NR,NTHL) &
                        -(WP(NS)/WF**2)
                   THR(NR,NTHL)=THR(NR,NTHL) &
                        -(WP(NS)/WF**2)/(1.D0-(WC(NS)/WF)**2) &
                        *FACTORHR
                   TRC(NR,NTHL)=TRC(NR,NTHL) &
                        -(WP(NS)/WF**2)/(1.D0+WC(NS)/WF) &
                        *FACTORRC
                   TLC(NR,NTHL)=TLC(NR,NTHL) &
                        -(WP(NS)/WF**2)/(1.D0-WC(NS)/WF) &
                        *FACTORLC
                ENDDO
                GTCO(NR,NTHL)=GUCLIP(TCO(NR,NTHL))
                GTHR(NR,NTHL)=GUCLIP(THR(NR,NTHL))
                GTRC(NR,NTHL)=GUCLIP(TRC(NR,NTHL))
                GTLC(NR,NTHL)=GUCLIP(TLC(NR,NTHL))
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!     ***** DRAW FRAME *****

    GRLEN=GUCLIP(RGMAX-RGMIN)
    GZLEN=GUCLIP(ZGMAX-ZGMIN)
    IF(GRLEN.GT.GZLEN) THEN
       GPR=15.0
       GPZ=15.0*GZLEN/GRLEN
    ELSE
       GPR=15.0*GRLEN/GZLEN
       GPZ=15.0
    ENDIF
    CALL GQSCAL(GUCLIP(RGMIN),GUCLIP(RGMAX),GGRMIN,GGRMAX,GGRSTP)
    CALL GQSCAL(GUCLIP(ZGMIN),GUCLIP(ZGMAX),GGZMIN,GGZMAX,GGZSTP)
    IRORG=(NINT(GGRMIN/GGRSTP)+3)/2*2
    GRORG=IRORG*GGRSTP

    CALL SETLIN(0,2,7)
    CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ, &
                GUCLIP(RGMIN),GUCLIP(RGMAX), &
                GUCLIP(ZGMIN),GUCLIP(ZGMAX))
    CALL GFRAME
    CALL GSCALE(GRORG,GGRSTP,0.0,0.0,0.1,9)
    CALL GVALUE(GRORG,GGRSTP*2,0.0,0.0,NGULEN(2*GGRSTP))
    CALL GSCALE(0.0,0.0,0.0,GGZSTP,0.1,9)
    CALL GVALUE(0.0,0.0,0.0,GGZSTP*2,NGULEN(2*GGZSTP))

!     ***** DRAW RESONANCE AND CUTOFF *****

    IF(K2.EQ.'P') THEN

!     ****** DRAW CYCROTRON REASONANCE SURFACE (yellow, dot-dashed) ****** 

       DO NS=1,NSMAX
          DO NC=1,4
             BCF=PA(NS)*AMP*WF/(NC*ABS(PZ(NS))*AEE)
             GBCF=GUCLIP(BCF)
             CALL SETRGB(1.0,1.0,(NC-1)*0.1)
             CALL CONTQ5(GBY,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                         GBCF,1.0,1,2,4,KACONT)
          ENDDO
       ENDDO

!     ****** DRAW PLASMA CUTOFF SURFACE (medium blue, two-dot-dashed) ******

       CALL SETRGB(0.0,1.0,1.0)
       CALL CONTQ5(GTCO,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.0,1000.0,1,2,6,KACONT)

!     ****** DRAW RIGHT CUT OFF SURFACE (light blue, two-dot-dashed) ******

       CALL SETRGB(0.5,1.0,1.0)
       CALL CONTQ5(GTRC,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.0,1000.0,1,2,6,KACONT)

!     ****** DRAW LEFT CUT OFF SURFACE (light purple, two-dots-dashed) ******

       CALL SETRGB(1.0,0.5,1.0)
       CALL CONTQ5(GTLC,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.0,1000.0,1,2,6,KACONT)

!     ****** DRAW HYBRID RESONANCE SURFACE (purple, long-dashed) ******

       CALL SETRGB(1.0,0.0,1.0)
       CALL CONTQ5(GTHR,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.0,1000.0,1,2,3,KACONT)

    ENDIF

!     ****** DRAW MAIN DATA ****** 

!    CALL GMNMX2(GFL,nrmax+1,1,NRMAX+1,1,1,nthmax_g,1,GFMIN,GFMAX)
    DO NR=1,NRMAX
       IF (xrho(NR).LT. 1.D0) THEN
          NRLCFS=NR
       END IF
    ENDDO
    CALL GMNMX2(GFL,nrmax+1,1,NRLCFS,1,1,nthmax_g,1,GFMIN,GFMAX)

    IF(ABS(GFMIN).GT.GFMAX) GFMAX=ABS(GFMIN)
    GFMIN=-ABS(GFMAX)
    CALL GQSCAL(GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP)
    GGFSTP=0.5*GGFSTP
    NSTEP=INT((GGFMAX-GGFMIN)/GGFSTP)+1
    IF (NSTEP .LT. NCONT )THEN
       NSTEP=NCONT
       GGFSTP=(GGFMAX-GGFMIN)/REAL(NSTEP)
    ENDIF

!Chonda      IF(GFMIN*GFMAX.GT.0.0) THEN
    IF(DBLE(GFMIN)*DBLE(GFMAX).GT.0.d0) THEN ! To avoid floating overflow
       CALL SETLIN(-1,-1,6)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   GGFMIN,GGFSTP,NSTEP,2,0,KACONT)
    ELSE
       CALL SETLIN(-1,-1,6)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.5*GGFSTP, GGFSTP,NSTEP,2,0,KACONT)
       CALL SETLIN(-1,-1,5)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                  -0.5*GGFSTP,-GGFSTP,NSTEP,2,2,KACONT)
    ENDIF

!     ****** DRAW PLASMA AND WALL SURFACE ****** 

    DO NSU=1,NSUMAX
       GRS(NSU)=GUCLIP(RSU(NSU,NHH))
       GZS(NSU)=GUCLIP(ZSU(NSU,NHH))
    ENDDO
    GRS(NSUMAX+1)=GRS(1)
    GZS(NSUMAX+1)=GZS(1)
    CALL SETLIN(-1,-1,7)
    CALL GPLOTP(GRS,GZS,1,NSUMAX+1,1,0,0,0)
    DO NSU=1,NSUMAX
       GRS(NSU)=GUCLIP(RSW(NSU,NHH))
       GZS(NSU)=GUCLIP(ZSW(NSU,NHH))
    ENDDO
    GRS(NSUMAX+1)=GRS(1)
    GZS(NSUMAX+1)=GZS(1)
    CALL SETLIN(-1,-1,7)
    CALL GPLOTP(GRS,GZS,1,NSUMAX+1,1,0,0,0)

    CALL MOVE(17.5,16.0)
    CALL TEXT('MAX :',5)
    CALL NUMBR(GFMAX,'(1PE9.2)',9)
    CALL MOVE(17.5,15.5)
    CALL TEXT('MIN :',5)
    CALL NUMBR(GFMIN,'(1PE9.2)',9)
    CALL MOVE(17.5,15.0)
    CALL TEXT('STEP:',5)
    CALL NUMBR(GGFSTP,'(1PE9.2)',9)

    CALL MOVE(17.5,14.5)
    CALL TEXT('RAXS:',5)
    CALL NUMBR(REAL(RAXIS),'(1PE9.2)',9)
    CALL MOVE(17.5,14.0)
    CALL TEXT('ZAXS:',5)
    CALL NUMBR(REAL(ZAXIS),'(1PE9.2)',9)

    RETURN
  END SUBROUTINE wm_gxeq

!     ****** PAINT CONTOUR IN MAGNETIC SURFACE COORDINATES ******

  SUBROUTINE wm_gxeqp(GGL,NHH,K2,K3)

    USE wmcomm
    IMPLICIT NONE
    REAL(4),INTENT(IN):: GGL(nrmax,nthmax)
    INTEGER,INTENT(IN):: NHH
    CHARACTER(LEN=1):: K2,K3
    INTEGER,PARAMETER:: NRGBA=5
    INTEGER,PARAMETER:: NSTEPM=101
    REAL(4),ALLOCATABLE:: GBY(:,:),GFL(:,:),GRL(:,:),GZL(:,:),GRS(:),GZS(:)
    REAL(4),ALLOCATABLE:: GDL(:),GRGBL(:,:)
    REAL(4):: GRGBA(3,NRGBA),GLA(NRGBA)
    DATA GRGBA/0.0,0.0,1.0, &
               0.0,1.0,1.0, &
               1.0,1.0,1.0, &
               1.0,1.0,0.0, &
               1.0,0.0,0.0/
    DATA GLA/0.0,0.40,0.5,0.60,1.0/
    INTEGER:: NR,NTH,NHHF,NTHF,NTHGS,NTHP,NTHPF,NTHG,NTHL,NRLCFS
    INTEGER:: ISTEP,I,ITORG,NNHD,NSW,IRORG,NSTEP,NSU,NHHD
    REAL(rkind):: FACT,VAL,DTHG,RCOS,RSIN
    REAL(4):: GZA,GDZ,GFACT,GRLEN,GZLEN,GPR,GPZ
    REAL(4):: GGRMIN,GGRMAX,GGRSTP,GGZMIN,GGZMAX,GGZSTP,GRORG
    REAL(4):: GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP
    REAL(rkind):: BICF
    REAL(4):: GBICF
    
    ALLOCATE(GBY(nrmax,nthmax_f))
    ALLOCATE(GFL(nrmax,nthmax_f),GRL(nrmax,nthmax_f),GZL(nrmax,nthmax_f))
    ALLOCATE(GRS(nsumax+1),GZS(nsumax+1))
    ALLOCATE(GDL(NSTEPM),GRGBL(3,0:NSTEPM))
      
    SELECT CASE(MODELG)
    CASE(0:2)
       DTHG=2.D0*PI/nthmax_g
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             RCOS=COS(DTHG*(NTH-1))
             RSIN=SIN(DTHG*(NTH-1))
             IF(MODELG.EQ.0) THEN
                RPSG(NTH,NR)  =      XR(NR)*RCOS
             ELSE
                RPSG(NTH,NR)  = RR + XR(NR)*RCOS
             ENDIF
             ZPSG(NTH,NR)  =         XR(NR)*RSIN
          ENDDO
       END DO
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             GRL(NR,NTH)=GUCLIP(RPSG(NTH,NR))
             GZL(NR,NTH)=GUCLIP(ZPSG(NTH,NR))
          ENDDO
       ENDDO
    CASE(3)
       ! RPSG and ZPSG loaded in wmsetg_tor
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             GRL(NR,NTH)=GUCLIP(RPSG(NTH,NR))
             GZL(NR,NTH)=GUCLIP(ZPSG(NTH,NR))
          ENDDO
       ENDDO
    CASE(4,6)
       nthmax_g=NTHMAX
       DO NR=1,NRMAX+1
          DO NTH=1,nthmax_g
             NHHF=(NHH-1)*factor_nhh +1
             NTHF=(NTH-1)*factor_nth +1
             GRL(NR,NTH)=GUCLIP(RPST(NTHF,NHHF,NR))
             GZL(NR,NTH)=GUCLIP(ZPST(NTHF,NHHF,NR))
          ENDDO
       ENDDO
    END SELECT
    
    NTHGS=nthmax_g/NTHMAX
    DO NR=1,NRMAX+1
       DO NTH=1,NTHMAX
          NTHP=NTH+1
          NHHF=(NHH-1)*FACTOR_NHH+1
          NTHF=(NTH-1)*FACTOR_NTH+1
          NTHPF=(NTHP-1)*FACTOR_NTH+1
          IF(NTHP.GT.NTHMAX) NTHP=1
          DO NTHG=1,NTHGS
             NTHL=(NTH-1)*NTHGS+NTHG
             FACT=DBLE(NTHG-1)/DBLE(NTHGS)
             VAL=(1.D0-FACT)*BPST(NTHF,NHHF,NR)+FACT*BPST(NTHPF,NHHF,NR)
             GBY(NR,NTHL)=GUCLIP(VAL)
             VAL=(1.D0-FACT)*DBLE(GGL(NR,NTH))+FACT*DBLE(GGL(NR,NTHP))
             GFL(NR,NTHL)=GUCLIP(VAL)
          ENDDO
       ENDDO
    ENDDO

    CALL GMNMX2(GFL,nrmax+1,1,NRMAX+1,1,1,nthmax_g,1,GFMIN,GFMAX)
    DO NR=1,NRMAX
       IF(xrho(NR).LT.1.D0) THEN
          NRLCFS=NR
       END IF
    ENDDO
    CALL GMNMX2(GFL,nrmax+1,1,NRLCFS,1,1,nthmax_g,1,GFMIN,GFMAX)
    IF(ABS(GFMIN).GT.GFMAX) GFMAX=ABS(GFMIN)
    GFMIN=-ABS(GFMAX)

    CALL GQSCAL(GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP)

    ISTEP=NCONT

!      IF(GFMIN*GFMAX.LT.0.D0) THEN

    GZA=MAX(ABS(GGFMAX),ABS(GGFMIN))
    GDZ=2*GZA/ISTEP
    DO I=1,ISTEP
       GDL(I)=GDZ*(I-0.5)-GZA
    ENDDO

    DO I=0,ISTEP
       GFACT=REAL(I)/REAL(ISTEP)
       CALL GUSRGB(GFACT,GRGBL(1,I),NRGBA,GLA,GRGBA)
    ENDDO
!      ELSE
!         GDZ=(GGFMAX-GGFMIN)/ISTEP
!         GZA=GGFMIN
!         DO I=1,ISTEP
!            GDL(I)=GDZ*I+GZA
!         ENDDO
!         DO I=0,ISTEP
!            GFACT=REAL(I)/REAL(ISTEP)
!            CALL GUSRGB(GFACT,GRGBL(1,I),NRGBB,GLB,GRGBB)
!         ENDDO

!      ENDIF

    GRLEN=GUCLIP(RGMAX-RGMIN)
    GZLEN=GUCLIP(ZGMAX-ZGMIN)
    IF(GRLEN.GT.GZLEN) THEN
       GPR=15.0
       GPZ=15.0*GZLEN/GRLEN
    ELSE
       GPR=15.0*GRLEN/GZLEN
       GPZ=15.0
    ENDIF
    CALL GQSCAL(GUCLIP(RGMIN),GUCLIP(RGMAX),GGRMIN,GGRMAX,GGRSTP)
    CALL GQSCAL(GUCLIP(ZGMIN),GUCLIP(ZGMAX),GGZMIN,GGZMAX,GGZSTP)
    IRORG=(NINT(GGRMIN/GGRSTP)+3)/2*2
    GRORG=IRORG*GGRSTP

    CALL SETLIN(0,2,7)
    CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ, &
                GUCLIP(RGMIN),GUCLIP(RGMAX), &
                GUCLIP(ZGMIN),GUCLIP(ZGMAX))
    CALL GFRAME
    CALL GSCALE(GRORG,GGRSTP,0.0,0.0,0.1,9)
    CALL GVALUE(GRORG,GGRSTP*2,0.0,0.0,NGULEN(2*GGRSTP))
    CALL GSCALE(0.0,0.0,0.0,GGZSTP,0.1,9)
    CALL GVALUE(0.0,0.0,0.0,GGZSTP*2,NGULEN(2*GGZSTP))

    CALL CONTF6(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                GDL,GRGBL,ISTEP,2)

    IF(GFMIN*GFMAX.GT.0.0) THEN
       CALL SETLIN(0,0,6)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   GGFMIN,GGFSTP,NSTEP,2,0,KACONT)
    ELSE
!       CALL SETLIN(0,0,6)
!       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
!                   0.25*GGFSTP, GGFSTP,NSTEP,2,0,KACONT)
       CALL SETRGB(0.8,0.8,0.0)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   0.15*GGFSTP, GGFSTP,1,2,0,KACONT)
!       CALL SETLIN(0,0,5)
!       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
!                  -0.25*GGFSTP,-GGFSTP,NSTEP,2,2,KACONT)
       CALL SETRGB(0.0,0.8,0.8)
       CALL CONTQ5(GFL,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                  -0.15*GGFSTP,-GGFSTP,1,2,0,KACONT)
    ENDIF

    CALL RGBBAR(2.3+GPR,2.8+GPR,2.0,2.0+GPZ,GRGBL,ISTEP+1,1)

    IF(K2.EQ.'P'.AND.K3.EQ.'3') THEN

       BICF=2.D0*PI*AMP*RF*1.D6/AEE
       GBICF=GUCLIP(BICF)
       CALL SETRGB(0.0,1.0,0.0)
       CALL CONTQ5(GBY,GRL,GZL,nrmax+1,NRMAX+1,nthmax_g, &
                   GBICF,1.0,1,2,7,KACONT)

    ENDIF

    DO NSU=1,NSUMAX
       GRS(NSU)=GUCLIP(RSU(NSU,NHH))
       GZS(NSU)=GUCLIP(ZSU(NSU,NHH))
    ENDDO
    GRS(NSUMAX+1)=GRS(1)
    GZS(NSUMAX+1)=GZS(1)
    CALL SETRGB(0.0,1.0,0.0)
    CALL GPLOTP(GRS,GZS,1,NSUMAX+1,1,0,0,0)

    IF(MODELG.EQ.4.OR.MODELG.EQ.6) THEN
       NHHD=NHH
    ELSE
       NHHD=1
    ENDIF
    DO NSW=1,NSWMAX
       GRS(NSW)=GUCLIP(RSW(NSW,NHHD))
       GZS(NSW)=GUCLIP(ZSW(NSW,NHHD))
    ENDDO
    CALL SETRGB(0.0,0.0,0.0)
    CALL GPLOTP(GRS,GZS,1,NSWMAX,1,0,0,0)

    CALL SETLIN(0,2,7)
    CALL MOVE(18.5,16.0)
    CALL TEXT('MAX :',5)
    CALL MOVE(18.5,15.5)
    CALL NUMBR(GFMAX,'(1PE9.2)',9)
    CALL MOVE(18.5,15.0)
    CALL TEXT('MIN :',5)
    CALL MOVE(18.5,14.5)
    CALL NUMBR(GFMIN,'(1PE9.2)',9)
    CALL MOVE(18.5,14.0)
    CALL TEXT('STEP:',5)
    CALL MOVE(18.5,13.5)
    CALL NUMBR(GDZ,'(1PE9.2)',9)

    RETURN
  END SUBROUTINE wm_gxeqp

!     ****** DRAW BIRD-EYE VIEW IN MAGNETIC SURFACE COORDINATES ******

  SUBROUTINE wm_g3da(GFL,NHH,K2,K3,IND)

    USE wmcomm
    IMPLICIT NONE
    REAL(4),INTENT(IN):: GFL(nrmax,nthmax)
    INTEGER,INTENT(IN):: NHH,IND
    CHARACTER(LEN=1):: K2,K3
    INTEGER,PARAMETER:: NDIVM=101

    REAL(4),ALLOCATABLE:: GFM(:,:)
    REAL(rkind),ALLOCATABLE:: RBS(:,:),ZBS(:,:),FBS(:,:)
    REAL(rkind),ALLOCATABLE:: RSWF(:),ZSWF(:),RSWB(:),ZSWB(:)
    REAL(rkind),ALLOCATABLE:: UR(:,:,:,:),UZ(:,:,:,:),UF(:,:,:,:)
    REAL(rkind),ALLOCATABLE:: FX(:,:),FY(:,:),FXY(:,:)
    REAL(rkind),ALLOCATABLE:: U1(:,:),F1(:),U2(:,:),F2(:)
    REAL(rkind),ALLOCATABLE:: XMP(:),YMP(:)
    REAL(rkind):: A(2,2),B(2)

    INTEGER:: NR,NTH,NHHF,NTHF,NHHD,NSW,NRMN,NDIV,I,NSWMAX1,NSWMAX2,MN
    INTEGER:: NYP,IXFST,NXP,IERR,ICNVG,ITER
    INTEGER:: NY,NSU
    REAL(rkind):: DTH,RMIN,RMAX,ZMIN,ZMAX,RMIN1,RMAX1,ZMIN1,ZMAX1
    REAL(rkind):: XLNG,YLNG,DXM,DYM,EPS,RHO,THETA,RHOP,XMPP,YMPP
    REAL(rkind):: DRHO,DTHETA,ZUPL,ZDWL,XDIF,YDIF,FMP,DFX,DFY
    REAL(4):: GZMIN,GZMAX,GZMIN1,GZMAX1,GMAX,GXMIN,GXMAX,GYMIN,GYMAX
    REAL(4):: GXL,GYL,GRL,GZL,GPHI,GTHETA,GRADIUS,GXP,GYP

    INTERFACE
       SUBROUTINE R2W2B(R,RGB)
         REAL(4),INTENT(IN):: R
         REAL(4),INTENT(OUT):: RGB(3)
       END SUBROUTINE R2W2B
    END INTERFACE

    ALLOCATE(GFM(NDIVM,NDIVM))
    ALLOCATE(RBS(nrmax+1,nthmax+1),ZBS(nrmax+1,nthmax+1))
    ALLOCATE(FBS(nrmax+1,nthmax+1))
    ALLOCATE(RSWF(nsumax+1),ZSWF(nsumax+1))
    ALLOCATE(RSWB(nsumax+1),ZSWB(nsumax+1))
    ALLOCATE(UR(4,4,nrmax+1,nthmax+1),UZ(4,4,nrmax+1,nthmax+1))
    ALLOCATE(UF(4,4,nrmax+1,nthmax+1))
    ALLOCATE(FX(nrmax+1,nthmax+1),FY(nrmax+1,nthmax+1),FXY(nrmax+1,nthmax+1))
    ALLOCATE(U1(4,nsumax+1),F1(nsumax+1))
    ALLOCATE(U2(4,nsumax+1),F2(nsumax+1))
    ALLOCATE(XMP(NDIVM),YMP(NDIVM))

    NHHF=(NHH-1)*factor_nhh+1
    DO NR=1,NRMAX+1
       DO NTH=1,NTHMAX
          NTHF=(NTH-1)*factor_nth+1
          RBS(NR,NTH)=RPST(NTHF,NHHF,NR)
          ZBS(NR,NTH)=ZPST(NTHF,NHHF,NR)
          FBS(NR,NTH)=DBLE(GFL(NR,NTH))
       ENDDO
    ENDDO

    DO NR=1,NRMAX+1
       RBS(NR,NTHMAX+1)=RPST(1,NHHF,NR)
       ZBS(NR,NTHMAX+1)=ZPST(1,NHHF,NR)
       FBS(NR,NTHMAX+1)=DBLE(GFL(NR,1))
    ENDDO

    CALL SPL2D(XR,XTH,RBS,FX,FY,FXY,UR, &
               nrmax+1,NRMAX+1,NTHMAX+1,0,4,IERR)

    CALL SPL2D(XR,XTH,ZBS,FX,FY,FXY,UZ, &
               nrmax+1,NRMAX+1,NTHMAX+1,0,4,IERR)

    CALL SPL2D(XR,XTH,FBS,FX,FY,FXY,UF, &
               nrmax+1,NRMAX+1,NTHMAX+1,0,4,IERR)

    RMIN= 1.D32
    RMAX=-1.D-32
    ZMIN= 1.D32
    ZMAX=-1.D-32

    IF(MODELG.NE.4) THEN
       NHHD=1
    ELSE
       NHHD=NHH
    ENDIF
    DO NSW=1,NSWMAX
       RMIN1=RSW(NSW,NHHD)
       RMAX1=RSW(NSW,NHHD)
       ZMIN1=ZSW(NSW,NHHD)
       ZMAX1=ZSW(NSW,NHHD)
       RMIN=MIN(RMIN,RMIN1)
       RMAX=MAX(RMAX,RMAX1)
       ZMIN=MIN(ZMIN,ZMIN1)
       ZMAX=MAX(ZMAX,ZMAX1)
       IF(RMIN.EQ.RMIN1) THEN
          NRMN=NSW
       ENDIF
    ENDDO

    NDIV=80

    XLNG=RGMAX-RGMIN
    YLNG=ZGMAX-ZGMIN

    DXM =XLNG/NDIV
    DYM =YLNG/NDIV

    DO I=1,NDIV+1
       XMP(I)=RGMIN+DXM*(I-1)
       YMP(I)=ZGMIN+DYM*(I-1)
    ENDDO

    NSWMAX1=NRMN
    NSWMAX2=NSWMAX-NRMN+1

    MN=0
    DO I=NRMN,1,-1
       MN=MN+1
       RSWB(MN)=RSW(I,NHHD)
       ZSWB(MN)=ZSW(I,NHHD)
    ENDDO

    MN=0
    DO I=NRMN,NSWMAX
       MN=MN+1
       RSWF(MN)=RSW(I,NHHD)
       ZSWF(MN)=ZSW(I,NHHD)
    ENDDO

    CALL SPL1D(RSWB,ZSWB,F1,U1,NSWMAX1,0,IERR)
    CALL SPL1D(RSWF,ZSWF,F2,U2,NSWMAX2,0,IERR)

    EPS=1.D-8

    DO NYP=1,NDIV+1
       IXFST=0
       DO NXP=1,NDIV+1
          CALL SPL1DF(XMP(NXP),ZUPL,RSWB,U1,NSWMAX1,IERR)
          IF(IERR.EQ.2) THEN
             ZUPL = 0.D0
          ENDIF
          CALL SPL1DF(XMP(NXP),ZDWL,RSWF,U2,NSWMAX2,IERR)
          IF(IERR.EQ.2) THEN
             ZDWL = 0.D0
          ENDIF

          IF(YMP(NYP).LT.ZDWL.OR.YMP(NYP).GT.ZUPL) THEN
!               GXM(NXP)    =GUCLIP(XMP(NXP))
!               GYM(NYP)    =GUCLIP(YMP(NYP))
             GFM(NXP,NYP)=0.0
             GO TO 100
          ENDIF

!           ###### INITIAL VALUE ######

          IF(IXFST.EQ.0) THEN
!               RHO  =0.5D0
             RHO  =(RB/RA)/2.D0
             THETA=PI
          ENDIF

!           ###### ITERATION ######

          ICNVG=0
10        CONTINUE
          IF(ICNVG.NE.0) THEN
             RHO=RHOP+(RB/RA)/NRMAX*ICNVG
             THETA=PI
          ENDIF
          DO ITER=1,100
!               IF(NYP.EQ.NDIV/2+1) THETA=PI
             CALL SPL2DD(RHO,THETA,B(1),A(1,1),A(2,1),XR,XTH,UR, &
                         nrmax+1,nrmax+1,nthmax+1,IERR)
             IF(IERR.EQ.2.OR.IERR.EQ.4) THEN
                CALL SUBSPL1(RHO,THETA,B(1),A(1,1),A(2,1),UR,IERR)
             ENDIF
             CALL SPL2DD(RHO,THETA,B(2),A(1,2),A(2,2),XR,XTH,UZ, &
                         nrmax+1,nrmax+1,nthmax+1,IERR)
             IF(IERR.EQ.2.OR.IERR.EQ.4) THEN
                CALL SUBSPL1(RHO,THETA,B(2),A(1,2),A(2,2),UZ,IERR)
             ENDIF

             IF(ICNVG.EQ.1) THEN
!                  WRITE(6,*) 'ITER=',ITER,NXP,NYP,IERR
!                  WRITE(6,*) DRHO,DTHETA
!                  WRITE(6,*) RHO,THETA
!                  WRITE(6,*) A(1,1),A(2,1)
!                  WRITE(6,*) A(1,2),A(2,2)
!                  WRITE(6,*) B(1),XMP(NXP)
!                  WRITE(6,*) B(2),YMP(NYP)
             ELSE IF(ICNVG.EQ.2) THEN
                WRITE(6,*) 'NO CONVERGENCE : 2D NEWTON METHOD'
                RETURN
             ENDIF

             B(1)=-B(1)+XMP(NXP)
             B(2)=-B(2)+YMP(NYP)
!     
!           ####### CONVERGENCE CHECK ########
!     
             IF(ITER.GT.1) THEN
                XDIF=DABS(B(1)-XMPP)
                YDIF=DABS(B(2)-YMPP)
                IF(XDIF.LE.EPS.AND.YDIF.LE.EPS) THEN
                   RHOP=RHO
!                     THETAP=THETA
                   GO TO 20
                ENDIF
             ENDIF

             XMPP=B(1)
             YMPP=B(2)
!     
             DRHO=(A(2,2)*B(1)-A(2,1)*B(2)) &
                 /(A(1,1)*A(2,2)-A(2,1)*A(1,2))
             DTHETA=(B(2)-A(1,2)*DRHO)/A(2,2)

             RHO  =RHO  +DRHO
             THETA=THETA+DTHETA

             IF(RHO.LT.0.D0) THEN
                RHO=DABS(RHO)
             ENDIF
             IF(RHO.GT.XR(NRMAX+1)) THEN
!                RHO=0.1D0*XR(NRMAX+1)
             ENDIF
             IF(THETA.GE.2.D0*PI) THEN
                THETA=DMOD(THETA,2.D0*PI)
             ELSE IF(THETA.LT.0.D0) THEN
                THETA=2.D0*PI-DMOD(DABS(THETA),2.D0*PI)
             ENDIF

          ENDDO
          ICNVG=ICNVG+1
          GO TO 10

20        CONTINUE
          CALL SPL2DD(RHO,THETA,FMP,DFX,DFY,XR,XTH,UF, &
                      nrmax+1,nrmax+1,nthmax+1,IERR)
          IF(IERR.EQ.2.OR.IERR.EQ.4) THEN
             FMP = 0.D0
          ENDIF

          GFM(NXP,NYP)=GUCLIP(FMP     )
          IXFST=IXFST+1

100       CONTINUE
       ENDDO
    ENDDO

    GZMIN= 1.E32
    GZMAX=-1.E32
    DO NY=1,NDIV+1
       CALL GMNMX1(GFM(1,NY),1,NDIV+1,1,GZMIN1,GZMAX1)
       GZMIN=MIN(GZMIN,GZMIN1)
       GZMAX=MAX(GZMAX,GZMAX1)
    ENDDO
    IF(GZMAX*GZMIN.LE.0.0) THEN
       GMAX=MAX(ABS(GZMAX),ABS(GZMIN))
       GZMAX= GMAX
       GZMIN=-GMAX
    ENDIF

    GXMAX=GUCLIP(RGMAX)
    GXMIN=GUCLIP(RGMIN)
    GYMAX=GUCLIP(ZGMAX)
    GYMIN=GUCLIP(ZGMIN)

    CALL SETLIN(0,0,7)
!      CALL SETCHS(0.20,0.0)
!      CALL SETFNT(32)

    GXL=7.0
    GYL=GXL*REAL(RKAP)
    GZL= 3.0
    GPHI=-90.0
    GTHETA=30.0
    GRADIUS=100.0

    CALL GDEFIN3D(2.,22.,0.,20.,GXL,GYL,GZL)
    CALL GVIEW3D(GPHI,GTHETA,GRADIUS,0.9,1, &
                 0.5*(GXMIN+GXMAX),0.5*(GYMIN+GYMAX),0.0)
    CALL GDATA3D1(GFM,NDIVM,NDIV+1,NDIV+1, &
                  GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)

    CALL CPLOT3D1(IND,R2W2B)
!      CALL CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,MYR2G2B,2)
!      CALL PERSE3D(3,1)
    CALL PERSE3D(0,30)

!      CALL GAXIS3D(2)
!      CALL GSCALE3DX(0.0,GXSCAL,0.2,2)
!      CALL GSCALE3DY(0.0,GYSCAL,0.2,2)
!      CALL GSCALE3DZ(GZMIN,GZSCAL,0.2,2)
!      CALL GVALUE3DX(0.0,2.*GXSCAL,-6,1)
!      CALL GVALUE3DY(0.0,2.*GYSCAL,-6,1)
!      CALL GVALUE3DZ(0.0,2.*GZSCAL,-2,1)

    CALL SETRGB(0.0,1.0,0.0)
    NSU=1
    GRL=GUCLIP(RSU(NSU,NHH))
    GZL=GUCLIP(ZSU(NSU,NHH))
    CALL GTTTA(GRL,GZL,0.0,GXP,GYP)
    CALL MOVE(GXP,GYP)
    DO NSU=2,NSUMAX
       GRL=GUCLIP(RSU(NSU,NHH))
       GZL=GUCLIP(ZSU(NSU,NHH))
       CALL GTTTA(GRL,GZL,0.0,GXP,GYP)
       CALL DRAW(GXP,GYP)
    ENDDO
    NSU=1
    GRL=GUCLIP(RSU(NSU,NHH))
    GZL=GUCLIP(ZSU(NSU,NHH))
    CALL GTTTA(GRL,GZL,0.0,GXP,GYP)
    CALL DRAW(GXP,GYP)

    IF(MODELG.NE.4) THEN
       NHHD=1
    ELSE
       NHHD=NHH
    ENDIF
    CALL SETRGB(0.0,0.0,0.0)
    NSW=1
    GRL=GUCLIP(RSW(NSW,NHHD))
    GZL=GUCLIP(ZSW(NSW,NHHD))
    CALL GTTTA(GRL,GZL,0.0,GXP,GYP)
    CALL MOVE(GXP,GYP)
    DO NSW=2,NSWMAX
       GRL=GUCLIP(RSW(NSW,NHHD))
       GZL=GUCLIP(ZSW(NSW,NHHD))
       CALL GTTTA(GRL,GZL,0.0,GXP,GYP)
       CALL DRAW(GXP,GYP)
    ENDDO
    RETURN
  END SUBROUTINE wm_g3da

!
!     ****** WRITE PARAMETERS ******

  SUBROUTINE wm_gprm(K1,K3,NTH,NHH,MD,ND)

    USE wmcomm
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K1,K3
    INTEGER,INTENT(IN):: nth,nhh,md,nd
    REAL(4):: XPOS,YPOS,DY
    INTEGER:: NS,MDX,NDX

    CALL SETCHS(0.3,0.0)
    XPOS=22.0
    YPOS=17.5
    DY=0.4

    CALL MOVE(XPOS,YPOS)

    CALL TEXT('MP=',3)
    DO NS=1,NSMAX
       CALL NUMBI(MODELP(NS),'(I3)',3)
    END DO

!     ****** DEVICE PARAMETERS ******

    CALL MOVE(XPOS,YPOS-DY)
    CALL TEXT('BB: ',4)
    CALL NUMBD(BB,'(F8.4)',8)

    CALL MOVE(XPOS,YPOS-2*DY)
    CALL TEXT('RR: ',4)
    CALL NUMBD(RR,'(F8.4)',8)

    CALL MOVE(XPOS,YPOS-3*DY)
    CALL TEXT('Q0: ',4)
    CALL NUMBD(Q0,'(F8.4)',8)

    CALL MOVE(XPOS,YPOS-4*DY)
    CALL TEXT('QA: ',4)
    CALL NUMBD(QA,'(F8.4)',8)
    YPOS=YPOS-5*DY

!     ****** PLASMA PARAMETERS ******

    CALL MOVE(XPOS,YPOS)
    CALL TEXT('PN: ',4)
    DO NS=1,NSMAX
       CALL MOVE(XPOS,YPOS-NS*DY)
       CALL NUMBD(PN(NS),'(F12.4)',12)
    ENDDO

    YPOS=YPOS-7*DY
    CALL MOVE(XPOS,YPOS)
    CALL TEXT('PTPR: ',6)
    DO NS=1,NSMAX
       CALL MOVE(XPOS,YPOS-NS*DY)
       CALL NUMBD(PTPR(NS),'(F12.4)',12)
    ENDDO
    YPOS= YPOS-8*DY

!     ****** WAVE PARAMETERS ******

    IF(MODEEG.EQ.0) THEN
       CALL MOVE(XPOS,YPOS)
       CALL TEXT('RF:',3)
       CALL NUMBD(RF,'(F9.4)',9)
       YPOS= YPOS-2*DY
    ELSE
       CALL MOVE(XPOS,YPOS)
       CALL TEXT('RFR:',4)
       CALL MOVE(XPOS,YPOS-DY)
       CALL NUMBD(RF,'(1PD12.4)',12)
       CALL MOVE(XPOS,YPOS-2*DY)
       CALL TEXT('RFI:',4)
       CALL MOVE(XPOS,YPOS-3*DY)
       CALL NUMBD(RFI,'(1PD12.4)',12)
       CALL MOVE(XPOS,YPOS-4*DY)
       CALL TEXT('AMP:',4)
       CALL MOVE(XPOS,YPOS-5*DY)
       CALL NUMBD(AMP_EIGEN,'(1PD12.4)',12)
       YPOS= YPOS-7*DY
    ENDIF

    CALL MOVE(XPOS,YPOS)
    CALL TEXT('NPH0:',5)
    CALL NUMBI(NPH0,'(I7)',7)

    CALL MOVE(XPOS,YPOS-DY)
    CALL TEXT('NTH0:',5)
    CALL NUMBI(NTH0,'(I7)',7)

    CALL MOVE(XPOS,YPOS-2*DY)
    CALL TEXT('NRMAX:',6)
    CALL NUMBI(NRMAX,'(I6)',6)

    CALL MOVE(XPOS,YPOS-3*DY)
    CALL TEXT('NTHMAX:',7)
    CALL NUMBI(NTHMAX,'(I5)',5)

    CALL MOVE(XPOS,YPOS-4*DY)
    CALL TEXT('NHHMAX:',7)
    CALL NUMBI(NHHMAX,'(I5)',5)
    YPOS=YPOS-6*DY

!     ****** COMPUTATION RESULTS ******

    CALL MOVE(XPOS,YPOS)
    IF(K1.EQ.'R'.AND.K3.EQ.'M') THEN
       MDX=MD-MDMIN+1
       NDX=ND-NDMIN+1
       CALL TEXT('Pabs(MD):',9)
       DO NS=1,NSMAX
          CALL MOVE(XPOS,YPOS-NS*DY)
          CALL NUMBD(PABSKT(MDX,NDX,NS),'(1PD12.4)',12)
       ENDDO
    ELSE
       CALL TEXT('Pabs:',5)
       DO NS=1,NSMAX
          CALL MOVE(XPOS,YPOS-NS*DY)
          CALL NUMBD(PABST(NS),'(1PD12.4)',12)
       ENDDO
    ENDIF
    YPOS=YPOS-7*DY

    IF(MODEEG.EQ.0) THEN
       CALL MOVE(XPOS,YPOS)
       CALL TEXT('CUR:',4)
       CALL MOVE(XPOS,YPOS-DY)
       CALL NUMBD(PCURT,'(1PE12.4)',12)
       YPOS=YPOS-2*DY
       CALL MOVE(XPOS,YPOS)
       IF(K1.EQ.'R'.AND.K3.EQ.'M') THEN
          NDX=ND-NDMIN+1
          MDX=MD-MDMIN+1
          CALL TEXT('Imp(MD):',8)
          CALL MOVE(XPOS,YPOS-DY)
          CALL NUMBD(DBLE(CPRADK(MDX,NDX)),'(1PE12.4)',12)
          CALL MOVE(XPOS,YPOS-2*DY)
          CALL NUMBD(DIMAG(CPRADK(MDX,NDX)),'(1PE12.4)',12)
       ELSE
          CALL TEXT('Imp:',4)
          CALL MOVE(XPOS,YPOS-DY)
          CALL NUMBD(DBLE(CPRAD),'(1PE12.4)',12)
          CALL MOVE(XPOS,YPOS-2*DY)
          CALL NUMBD(DIMAG(CPRAD),'(1PE12.4)',12)
       ENDIF
       YPOS=YPOS-4*DY
    ENDIF

!     ****** THETA/MODE ******

    CALL MOVE(XPOS,YPOS)
    IF(K1.EQ.'R') THEN
       IF(K3.EQ.'M') THEN
          CALL TEXT('MD:',3)
          CALL NUMBI(MD,'(I7)',7)
          CALL MOVE(XPOS,YPOS-DY)
          CALL TEXT('ND:',3)
          CALL NUMBI(ND,'(I7)',7)
       ELSEIF(K3.EQ.'T'.OR.K3.EQ.'A') THEN
          CALL TEXT('Th:',3)
          CALL NUMBD(DBLE(NTH-1)/DBLE(NTHMAX)*360.D0,'(F9.4)',9)
          CALL MOVE(XPOS,YPOS-DY)
          CALL TEXT('Ph:',3)
          CALL NUMBD(DBLE(NHH-1)/DBLE(NHHMAX)*360.D0,'(F9.4)',9)
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE wm_gprm

!     ****** SET DATA COLOR ******

  SUBROUTINE MYR2G2B(R,RGB)

    USE wmcomm
    IMPLICIT NONE
    REAL(4),INTENT(IN):: R
    REAL(4),INTENT(OUT):: RGB(3)
    INTEGER:: IR

    IR=INT(12.0*R)
    IF (IR.LE.4) THEN
       RGB(1)=0.0
       RGB(2)=REAL(IR)/3.0
       RGB(3)=1.0
    ELSE IF (IR.LE.6) THEN
       RGB(1)=0.0
       RGB(2)=1.0
       RGB(3)=1.0-REAL(IR-3)/3.0
    ELSE IF (IR.LE.9) THEN
       RGB(1)=REAL(IR-6)/3.0
       RGB(2)=1.0
       RGB(3)=0.0
    ELSE
       RGB(1)=1.0
       RGB(2)=1.0-REAL(IR-9)/3.0
       RGB(3)=0.0
    ENDIF

    RETURN
  END SUBROUTINE MYR2G2B

!     ****** SET DATA COLOR ******

  SUBROUTINE MYR2W2B(R,RGB)

    USE wmcomm
    IMPLICIT NONE
    REAL(4),INTENT(IN):: R
    REAL(4),INTENT(OUT):: RGB(3)
    REAL(4):: X,F
    INTEGER:: IR

!      0 < R < 1
!     -1 < X < 1

    X=2*R-1.0
    IF(X.LT.0.0) THEN
       F=1.0+(X-0.3*X*X*X)/0.7
       RGB(1)=F
       RGB(2)=F
       RGB(3)=1.0
    ELSE
       F=1.0-(X-0.3*X*X*X)/0.7
       RGB(1)=1.0
       RGB(2)=F
       RGB(3)=F
    ENDIF

    RETURN
  END SUBROUTINE MYR2W2B

!     ****** CALCULATE THE VALUE OF MESH POINTS BY SPLINE ******
!
!  SUBROUTINE MESHVL(RHO,THETA,U,F,DFX,DFY,IERR)
!
!      INCLUDE 'wmcomm.inc'
!
!      COMMON /WMSPL1/ TH(NTHM+1)
!
!      DIMENSION U(4,4,NRM,NTHM+1)
!
!      CALL SPL2DD(RHO,THETA,F,DFX,DFY,XR,TH,U,NRM,NRMAX+1,NTHMAX+1,IERR)
!
!      RETURN
!      END

!     ****** EXCEED THE MAX POINT ******

  SUBROUTINE SUBSPL1(X,Y,Z,DZX,DZY,U,IERR)

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),intent(IN):: X,Y
    REAL(rkind),INTENT(IN):: U(4,4,nrmax+1,nthmax+1)
    REAL(rkind),INTENT(OUT):: Z,DZX,DZY
    INTEGER,INTENT(OUT):: IERR

    IF(IERR.EQ.1) THEN
       CALL SPL2DD(XR(1),Y,Z,DZX,DZY,XR,XTH,U, &
                     nrmax+1,NRMAX+1,NTHMAX+1,IERR)
       IF(IERR.EQ.0) THEN
          Z=Z+DZX*(X-XR(1))
       ELSE IF(IERR.EQ.4) THEN
          CALL SPL2DD(XR(1),XTH(NTHMAX+1),Z,DZX,DZY,XR,XTH,U, &
                        nrmax+1,NRMAX+1,NTHMAX+1,IERR)
          Z=Z+DZX*(X-XR(1))+DZY*(Y-XTH(NTHMAX+1))
       ELSE
          WRITE(6,*) 'XXX : SPLINE IERR=',IERR
       ENDIF

    ELSE IF(IERR.EQ.2) THEN
       CALL SPL2DD(XR(NRMAX+1),Y,Z,DZX,DZY,XR,XTH,U, &
                   nrmax+1,NRMAX+1,NTHMAX+1,IERR)
       IF(IERR.EQ.0) THEN
          Z=Z+DZX*(X-XR(NRMAX+1))
       ELSE IF(IERR.EQ.4) THEN
          CALL SPL2DD(XR(NRMAX+1),XTH(NTHMAX+1),Z,DZX,DZY,XR,XTH,U, &
                      nrmax+1,NRMAX+1,NTHMAX+1,IERR)
          Z=Z+DZX*(X-XR(NRMAX+1))+DZY*(Y-XTH(NTHMAX+1))
       ELSE
          WRITE(6,*) 'XXX : SPLINE IERR=',IERR
       ENDIF

    ELSE IF(IERR.EQ.4) THEN
       CALL SPL2DD(X,XTH(NTHMAX+1),Z,DZX,DZY,XR,XTH,U, &
                   nrmax+1,NRMAX+1,NTHMAX+1,IERR)
       IF(IERR.EQ.0) THEN
          Z=Z+DZY*(Y-XTH(NTHMAX+1))
       ELSE
          WRITE(6,*) 'XXX : SPLINE IERR',IERR
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE SUBSPL1

!     ****** CONTOUR PAINT : XY, VARIABLE STEP ******

  SUBROUTINE CONTF6(Z,X,Y,NXA,NXMAX,NYMAX,ZL,RGB,NSTEP,IPRD)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NXA,NXMAX,NYMAX,NSTEP,IPRD
    REAL(4),INTENT(IN):: Z(NXA,NYMAX),X(NXA,NYMAX),Y(NXA,NYMAX)
    REAL(4),INTENT(IN):: ZL(NSTEP),RGB(3,NSTEP)
    INTERFACE
       SUBROUTINE CONTV1(XA,YA,N,XB,YB,M,NN)
         INTEGER,INTENT(IN):: N,M
         REAL(4),INTENT(IN):: XA(N),YA(N)
         REAL(4),INTENT(OUT):: XB(M),YB(M)
         INTEGER,INTENT(OUT):: NN
       END SUBROUTINE CONTV1
    END INTERFACE

    CALL CONTFD(Z,X,Y,NXA,NXMAX,NYMAX,&
                ZL,RGB,NSTEP,IPRD,3,CONTV1)

    RETURN
  END SUBROUTINE CONTF6

!     ****** CONTOUR PAINT : COMMON SUB *******

  SUBROUTINE CONTFD(Z,X,Y,NXA,NX,NY,ZL,RGB,NSTEP,IPRD,INDX,SUBV)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NXA,NX,NY,NSTEP,IPRD,INDX
    REAL(4),INTENT(IN):: Z(NXA,NY),X(NXA*NY),Y(NXA*NY)
    REAL(4),INTENT(IN):: ZL(NSTEP),RGB(3,0:NSTEP)
    REAL(4):: XA(4),YA(4),ZA(4)
    REAL(4):: XG(101),YG(101)
    INTEGER,PARAMETER:: NH=101
    REAL(4):: XT(8,0:NH),YT(8,0:NH),HT(0:NH+1)
    INTEGER:: JH(0:NH)
    REAL(4):: RGBS(3,0:NH)
    INTERFACE
       SUBROUTINE SUBV(XA,YA,N,XB,YB,M,NN)
         INTEGER,INTENT(IN):: N,M
         REAL(4),INTENT(IN):: XA(N),YA(N)
         REAL(4),INTENT(OUT):: XB(M),YB(M)
         INTEGER,INTENT(OUT):: NN
       END SUBROUTINE SUBV
    END INTERFACE
    INTEGER:: KMAX,J,I,NXMAX,NYMAX,JJ,II,IH,NN
    REAL(4):: RS,GS,BS

    CALL INQRGB(RS,GS,BS)

    KMAX=NSTEP
    IF(KMAX.GT.NH) KMAX=NH

    IF(ZL(1).LE.ZL(KMAX)) THEN
       RGBS(1,0)=RGB(1,0)
       RGBS(2,0)=RGB(2,0)
       RGBS(3,0)=RGB(3,0)
       DO J=1,KMAX
          HT(J)=ZL(J)
          RGBS(1,J)=RGB(1,J)
          RGBS(2,J)=RGB(2,J)
          RGBS(3,J)=RGB(3,J)
       ENDDO
    ELSE
       RGBS(1,0)=RGB(1,KMAX)
       RGBS(2,0)=RGB(2,KMAX)
       RGBS(3,0)=RGB(3,KMAX)
       DO J=1,KMAX
          HT(J)=ZL(KMAX-J+1)
          RGBS(1,J)=RGB(1,KMAX-J+1)
          RGBS(2,J)=RGB(2,KMAX-J+1)
          RGBS(3,J)=RGB(3,KMAX-J+1)
       ENDDO
    ENDIF

    HT(0)=-1.E35
    HT(KMAX+1)=1.E35

    NXMAX=NX-1
    NYMAX=NY-1
    IF(IPRD.EQ.1.OR.IPRD.EQ.3) NXMAX=NX
    IF(IPRD.EQ.2.OR.IPRD.EQ.3) NYMAX=NY

    DO J=1,NYMAX
       DO I=1,NXMAX
          IF(INDX.EQ.0.OR.INDX.EQ.2) THEN
             XA(1)=X(1)*(NXA*(J-1)+I-1)+X(2)
             XA(2)=X(1)*(NXA*(J-1)+I  )+X(2)
          ELSE
             XA(1)=X(NXA*(J-1)+I  )
             XA(2)=X(NXA*(J-1)+I+1)
             IF(J.EQ.NYMAX.AND.(IPRD.EQ.2.OR.IPRD.EQ.3)) THEN
                XA(3)=X(I+1)
                XA(4)=X(I)
             ELSE
                XA(3)=X(NXA*(J-1)+I+1+NXA)
                XA(4)=X(NXA*(J-1)+I  +NXA)
             ENDIF
          ENDIF
          IF(INDX.EQ.0.OR.INDX.EQ.1) THEN
             YA(1)=Y(1)*(NXA*(J-1)+I-1)+Y(2)
             YA(3)=Y(1)*(NXA*(J-1)+I  )+Y(2)
          ELSE
             YA(1)=Y(NXA*(J-1)+I  )
             YA(2)=Y(NXA*(J-1)+I+1)
             IF(J.EQ.NYMAX.AND.(IPRD.EQ.2.OR.IPRD.EQ.3)) THEN
                YA(3)=Y(I+1)
                YA(4)=Y(I)
             ELSE
                YA(3)=Y(NXA*(J-1)+I+1+NXA)
                YA(4)=Y(NXA*(J-1)+I  +NXA)
             ENDIF
          ENDIF

          JJ=J+1
          IF(J.EQ.NY) JJ=1
          II=I+1
          IF(I.EQ.NX) II=1
          ZA(1)=Z(I, J )
          ZA(2)=Z(II,J )
          ZA(3)=Z(II,JJ)
          ZA(4)=Z(I, JJ)

          CALL CONTFX(XA,YA,ZA,4,XT,YT,JH,HT,KMAX)

          DO IH=0,KMAX
             IF(JH(IH).GE.3) THEN
                CALL SETRGB(RGBS(1,IH),RGBS(2,IH),RGBS(3,IH))
                CALL SUBV(XT(1,IH),YT(1,IH),JH(IH),XG,YG,100,NN)
                XG(NN+1)=XG(1)
                YG(NN+1)=YG(1)
                CALL POLY(XG,YG,NN+1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    CALL SETRGB(RS,GS,BS)
    RETURN
  END SUBROUTINE CONTFD

!     ****** DRAW METRIC GRAPH ******

  SUBROUTINE wm_ggr(GY,NGMAX,RTITL,GP)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NGMAX
    REAL(4),INTENT(IN):: GY(nrmax+1,nthmax),GP(4)
    CHARACTER(LEN=6),INTENT(IN):: RTITL
    INTEGER:: ICL(6)=[7,6,5,4,3,2]
    INTEGER:: IPAT(6)=[0,2,4,6,3,1]
    REAL(4):: GX(nrmax+1)
    INTEGER:: NR,I,ILD
    REAL(4):: GXMIN,GXMAX,GYMIN,GYMAX,GSXMIN,GSXMAX,GSX,GSYMIN,GSYMAX,GSY

    GXMIN=GUCLIP(XRHO(1))
    GXMAX=GUCLIP(XRHO(NRMAX+1))
    DO NR=1,NRMAX+1
       GX(NR)=GUCLIP(XRHO(NR))
    ENDDO
    CALL GMNMX2(GY,nrmax+1,1,NRMAX+1,1,1,NGMAX,1,GYMIN,GYMAX)

    CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)

    IF(GYMIN*GYMAX.GT.0.0) THEN
       IF(GYMIN.GT.0.0) THEN
          GYMIN=0.0
       ELSE
          GYMAX=0.0
       ENDIF
    ENDIF
    CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)

    IF(ABS(GSYMAX-GSYMIN).LT.1.E-32) GOTO 9000

    CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),GXMIN,GXMAX,GSYMIN,GSYMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
    CALL GVALUE(GXMIN,GSX*5,0.0,0.0,NGULEN(GSX*5))
    CALL GSCALE(0.0,0.0,0.0,GSY,0.1,9)
    CALL GVALUE(0.0,0.0,0.0,GSY*2,NGULEN(GSY*2))

    DO I=1,NGMAX
       ILD=MOD(I-1,6)+1
       CALL SETLIN(0,2,ICL(ILD))
       CALL GPLOTP(GX,GY(1,I),1,NRMAX+1,1,0,0,IPAT(ILD))
    ENDDO
    CALL SETLIN(0,2,7)

9000 CONTINUE
    CALL MOVE(GP(1),GP(4)+0.2)
    CALL TEXT(RTITL,6)

    RETURN
  END SUBROUTINE wm_ggr

  !     ****** DRAW CIRCLE ******

  SUBROUTINE wm_circle(RMAX)

    IMPLICIT NONE
    REAL(4),INTENT(IN):: RMAX
    REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX
    REAL(4):: DX,DY,DT1,X1,Y1
    INTEGER:: I
    
    CALL INQGDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
    DX=(PXMAX-PXMIN)/(GXMAX-GXMIN)
    DY=(PYMAX-PYMIN)/(GYMAX-GYMIN)

    DT1=2*3.1415926/100
    X1=DX*(RMAX-GXMIN)+PXMIN
    Y1=DY*(    -GYMIN)+PYMIN
    CALL MOVE(X1,Y1)
    DO I=1,101
       X1=DX*(RMAX*COS(DT1*(I-1))-GXMIN)+PXMIN
       Y1=DY*(RMAX*SIN(DT1*(I-1))-GYMIN)+PYMIN
       CALL DRAW(X1,Y1)
    ENDDO

    RETURN
  END SUBROUTINE wm_circle
END MODULE wmgsub
