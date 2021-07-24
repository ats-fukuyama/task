MODULE w1grf1

CONTAINS

!     ****** 1-D GRAPHIC ROUTINE ******

  SUBROUTINE W1GRUD
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZ,NX,NN,NS,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL
    REAL:: GAMIN1,GAMIN2,GAMIN3,GAMAX1,GAMAX2,GAMAX3

    IF(NGRAPH.EQ.0) RETURN

!    IF(PANT.EQ.0.D0) THEN
       RNORM=1.D0
!    ELSE
!       RNORM=1.D0/SQRT(PANT)
!    END IF
    RNORM2=RNORM*RNORM

1   CONTINUE
    IF(NZMAX.EQ.1) THEN
       NZ=1
    ELSE
100    CONTINUE
       WRITE(6,*) '## GRAPH INPUT : NZ  (0 : QUIT GRAPHIC)'
       READ(5,*,END=9990) NZ
       IF(NZ.GT.NZMAX) GOTO 100
       IF(NZ.LE.0) RETURN
    ENDIF

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,2,4)
    CALL MOVE( 2.5,17.7)
    CALL TEXT('EX,EY,EZ,PABS',13)

    GA=RA
    GB=RB
    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GX    (NX)=XAM(NN)
       GXM   (NX)=XAM(NX)
       GDATA1(NX)= ABS(CE2DA(NZ,NN,1))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,1))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,1))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,13.5,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,7)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,4)

    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GDATA1(NX)= ABS(CE2DA(NZ,NN,2))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,2))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,2))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,9.5,13.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,7)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,4)

    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GDATA1(NX)= ABS(CE2DA(NZ,NN,3))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,3))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,3))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,5.5,9.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,7)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,4)

    GMIN=0.0
    GMAX=0.0
    DO NS=1,NSMAX
       SELECT CASE(NMODEL)
       CASE(0:5,12)
          DO NX=1,NXPMAX
             GDATA1(NX)=PABSX(NX,NS)*RNORM2
          END DO
          CALL GMNMX1(GDATA1,2,NXPMAX-1,1,GAMIN,GAMAX)
       CASE(6:11)
          DO NX=1,NXMAX
             GDATA1(NX)=PABSX(NX,NS)*RNORM2
          END DO
          CALL GMNMX1(GDATA1,2,NXMAX-1,1,GAMIN,GAMAX)
       END SELECT
       GMIN=MIN(GMIN,GAMIN)
       GMAX=MAX(GMAX,GAMAX)
    END DO
    GAMAX=MAX(ABS(GMAX),ABS(GMIN))
    CALL GQSCAL(0.,GAMAX,GMIN,GMAX,SCAL)

    CALL GDEFIN(2.5,12.5,1.5,5.5,-GB,GB,-0.1*GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,0.,0.1,9)
    CALL GSCALE(0.,0.,0.,SCAL,0.1,1)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,20.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       CALL SETLIN(0,2,7-NS)
       SELECT CASE(NMODEL)
       CASE(0:5,12)
          DO NX=1,NXPMAX
             GDATA1(NX)=PABSX(NX,NS)*RNORM2
          END DO
          CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
       CASE(6:11)
          DO NX=1,NXMAX
             GDATA1(NX)=PABSX(NX,NS)*RNORM2
          END DO
          CALL GPLOTP(GXM,GDATA1,1,NXMAX,1,-NS,10,0)
       END SELECT
    END DO
    CALL SETLIN(0,2,4)

    IF(NMODEL.EQ.5.OR.NMODEL.EQ.3.OR.NMODEL.EQ.1) THEN
       DO NX=1,NXMAX
          NN        =NXPRNT(NX)
          GDATA1(NX)=FLUXX(NN)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
       IF(ABS(GAMIN).GE.ABS(GAMAX)) THEN
          DO NX=1,NXMAX
             GDATA1(NX)=-GDATA1(NX)
          END DO
          GAMAX=-GAMIN
       ENDIF
       CALL GQSCAL(0.,GAMAX*1.02,GMIN,GMAX,SCAL)
       CALL GDEFIN(2.5,12.5,1.5,5.5,-GB,GB,-0.1*GMAX,GMAX)
       CALL GSCALE(0.,0.,0.,SCAL,0.1,5)
       CALL SETLIN(0,2,7)
       CALL GPLOTP(GX, GDATA1,1,NXMAX,1,  0, 0,0)
    END IF

    CALL W1GPRM(NZ)
    CALL W1GPRX

    SELECT CASE(NMODEL)
    CASE(0:5,12)
       DO NX=1,NXPMAX
          GDATA1(NX)=(AHL(NX,1,3)+AHL(NX,1,4))*RNORM2
          GDATA2(NX)=AHL(NX,1,3)*RNORM2
          GDATA3(NX)=AHL(NX,1,4)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       CALL GMNMX1(GDATA2,1,NXPMAX,1,GAMIN2,GAMAX2)
       CALL GMNMX1(GDATA3,1,NXPMAX,1,GAMIN3,GAMAX3)
       GAMIN=MIN(GAMIN1,GAMIN2,GAMIN3)
       GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3)
       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
       CALL GDEFIN(15.0,25.0,1.5,5.5,-GB,GB,GMIN,GMAX)
       CALL GFRAME
       CALL GSCALE(0.,0.2,0.,2*SCAL,0.1,9)
       CALL GSCALE(0.,100.,0.,0.,0.2,9)
       CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
       CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
       CALL GVALUE(0.,REAL(RA),0.,0.,2)
       CALL SETLIN(0,2,7)
       !     CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
       CALL SETLIN(0,2,6)
       CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
       CALL SETLIN(0,2,5)
       CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1, 0, 0,0)
       CALL SETLIN(0,2,4)
    CASE(6:11)
       DO NX=1,NXMAX
          GDATA1(NX)=(AHL(NX,1,3)+AHL(NX,1,4))*RNORM2
          GDATA2(NX)=AHL(NX,1,3)*RNORM2
          GDATA3(NX)=AHL(NX,1,4)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN1,GAMAX1)
       CALL GMNMX1(GDATA2,1,NXMAX,1,GAMIN2,GAMAX2)
       CALL GMNMX1(GDATA3,1,NXMAX,1,GAMIN3,GAMAX3)
       GAMIN=MIN(GAMIN1,GAMIN2,GAMIN3)
       GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3)
       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
       CALL GDEFIN(15.0,25.0,1.5,5.5,-GB,GB,GMIN,GMAX)
       CALL GFRAME
       CALL GSCALE(0.,0.2,0.,2*SCAL,0.1,9)
       CALL GSCALE(0.,100.,0.,0.,0.2,9)
       CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
       CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
       CALL GVALUE(0.,REAL(RA),0.,0.,2)
       CALL SETLIN(0,2,7)
       !     CALL GPLOTP(GXM,GDATA1,1,NXMAX,1, 0, 0,0)
       CALL SETLIN(0,2,6)
       CALL GPLOTP(GXM,GDATA2,1,NXMAX,1,-7,10,0)
       CALL SETLIN(0,2,5)
       CALL GPLOTP(GXM,GDATA3,1,NXMAX,1, 0, 0,0)
       CALL SETLIN(0,2,4)
    END SELECT

    CALL PAGEE
    IF(NZMAX.GT.1) GOTO 1
    RETURN

9990 CONTINUE
    RETURN
  END SUBROUTINE W1GRUD

!     ****** 1-D GRAPHIC PARAMETER ******

  SUBROUTINE W1GPRM(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NS
    REAL(rkind):: RNORM2

    IF(ABS(PANT).LE.1.D-30) THEN
       RNORM2=1.D0
    ELSE
       RNORM2=1.D0/PANT
    ENDIF

    CALL SETLIN(0,0,7)
    CALL MOVE(14.,17.5)
    CALL TEXT('FREQ =',6)
    CALL NUMBD(RF,'(F9.4)',9)
    IF(NZMAX.EQ.1) THEN
       CALL TEXT('   K-PR =',9)
       CALL NUMBD(RKZ,'(F9.4)',9)
    ELSE
       CALL TEXT('   Z-POS=',9)
       CALL NUMBD(ZA(NZ),'(F9.4)',9)
    ENDIF
    CALL MOVE(14.,17.0)
    CALL TEXT('B    =',6)
    CALL NUMBD(BB,'(F9.4)',9)
    CALL TEXT('   WALLR=',9)
    CALL NUMBD(WALLR,'(1PE9.1)',9)

    CALL MOVE(14.,16.4)
    CALL TEXT('J-HI =',6)
    CALL NUMBD(AJYH(1),'(F9.4)',9)
    CALL TEXT('/',1)
    CALL NUMBD(RANT1(1),'(F8.4)',8)
    CALL NUMBD(XANT1(1),'(F9.4)',9)
    CALL MOVE(14.,15.9)
    CALL TEXT('J-LO =',6)
    CALL NUMBD(AJYL(1),'(F9.4)',9)
    CALL TEXT('/',1)
    CALL NUMBD(RANT2(1),'(F8.4)',8)
    CALL NUMBD(XANT2(1),'(F9.4)',9)

    CALL MOVE(14.,15.3)
    CALL TEXT('      ',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL TEXT('     NS',7)
       CALL NUMBI(NS,'(I1)',1)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14.,14.8)
    CALL TEXT('N-20 =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PN(NS),'(F8.3)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14.,14.3)
    CALL TEXT('T-PP =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PTPP(NS),'(F8.3)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14.,13.8)
    CALL TEXT('T-PR =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PTPR(NS),'(F8.3)',8)
    END DO

    IF(NZMAX.EQ.1.AND.NZ.EQ.0) RETURN

    CALL SETLIN(0,0,7)
    CALL MOVE(14.,13.3)
    CALL TEXT('P/TOT=',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PABSXZ(NS)*RNORM2,'(F8.3)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(14.,12.7)
    CALL TEXT('PTOT =',6)
    CALL NUMBD(PABSTT*RNORM2,'(F9.4)',9)
    CALL TEXT('   JTOT =',9)
    CALL NUMBD(AJCDT*RNORM2,'(F9.4)',9)

    RETURN
  END SUBROUTINE W1GPRM

!     ****** 2-D GRAPHIC ROUTINE ******

  SUBROUTINE W1GR2DW(MODE_2D)
    USE w1comm,ONLY: rkind,XAM,NXMAX,NZMAX,NXPRNT,CE2DA,PANT,RB,ZA,RZ
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE_2D
    INTEGER:: NX,NZ
    REAL(rkind):: RNORM
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: GX,GZ
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: GFR,GFI
    REAL(rkind),DIMENSION(1):: THICK,THIN

    THIN(1)=0.02D0
    THICK(1)=0.035D0

!    IF(PANT.EQ.0.D0) THEN
       RNORM=1.D0
!    ELSE
!       RNORM=1.D0/SQRT(PANT)
!    END IF
    ALLOCATE(GX(NXMAX),GZ(NZMAX))
    ALLOCATE(GFR(NXMAX,NZMAX),GFI(NXMAX,NZMAX))

    DO NX=1,NXMAX
       GX(NX)=XAM(NX)
    END DO

    DO NZ=1,NZMAX
       GZ(NZ)=ZA(NZ)
    END DO

    CALL PAGES

    DO NX=1,NXMAX
       DO NZ=1,NZMAX
          GFR(NX,NZ)= DBLE(CE2DA(NZ,NX,1))*RNORM
          GFI(NX,NZ)=AIMAG(CE2DA(NZ,NX,1))*RNORM
       END DO
    END DO

    CALL GRD2D(1,GX,GZ,GFR,NXMAX,NXMAX,NZMAX,'@Re(Ex)@', &
               XMIN=-RB/3.d0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=2.D0,GPXMAX=12.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               NLMAX=30,LINE_THICKNESS=THIN,MODE_2D=MODE_2D)
    CALL GRD2D(2,GX,GZ,GFI,NXMAX,NXMAX,NZMAX,'@Im(Ex)@', &
               XMIN=-RB/3.D0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=14.D0,GPXMAX=24.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               NLMAX=30,LINE_THICKNESS=THIN,MODE_2D=MODE_2D)

    CALL PAGEE

    CALL PAGES

    DO NX=1,NXMAX
       DO NZ=1,NZMAX
          GFR(NX,NZ)= DBLE(CE2DA(NZ,NX,2))*RNORM
          GFI(NX,NZ)=AIMAG(CE2DA(NZ,NX,2))*RNORM
       END DO
    END DO

    CALL GRD2D(3,GX,GZ,GFR,NXMAX,NXMAX,NZMAX,'@Re(Ey)@', &
               XMIN=-RB/3.D0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=2.D0,GPXMAX=12.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               LINE_THICKNESS=THICK,MODE_2D=MODE_2D)
    CALL GRD2D(4,GX,GZ,GFI,NXMAX,NXMAX,NZMAX,'@Im(Ey)@', &
               XMIN=-RB/3.D0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=14.D0,GPXMAX=24.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               LINE_THICKNESS=THICK,MODE_2D=MODE_2D)

    CALL PAGEE

    CALL PAGES

    DO NX=1,NXMAX
       DO NZ=1,NZMAX
          GFR(NX,NZ)= DBLE(CE2DA(NZ,NX,3))*RNORM
          GFI(NX,NZ)=AIMAG(CE2DA(NZ,NX,3))*RNORM
       END DO
    END DO

    CALL GRD2D(1,GX,GZ,GFR,NXMAX,NXMAX,NZMAX,'@Re(Ez)@', &
               XMIN=-RB/3.D0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=2.D0,GPXMAX=12.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               LINE_THICKNESS=THICK,MODE_2D=MODE_2D)
    CALL GRD2D(2,GX,GZ,GFI,NXMAX,NXMAX,NZMAX,'@Im(Ez)@', &
               XMIN=-RB/3.D0,XMAX=RB,YMIN=0.D0,YMAX=RZ, &
               GPXMIN=14.D0,GPXMAX=24.D0,GPYMIN=2.D0,GPYMAX=17.D0, &
               LINE_THICKNESS=THICK,MODE_2D=MODE_2D)

    CALL PAGEE

    DEALLOCATE(GX,GZ,GFR,GFI)
    RETURN
  END SUBROUTINE W1GR2DW

  SUBROUTINE W1GR2DR
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NS,NX,NZ,NZZ,NN,NGULEN
    REAL:: RNORM2,PGAMIN,PGAMAX
    REAL:: GA,GB,GZMIN,GZMAX,GZSCL,GAMIN,GAMAX,GMIN,GMAX,GSCAL
    REAL:: PGMIN,PGMAX,PSCAL,SGAMAX

    GA=RA
    GB=RB
    CALL PAGES
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,0,7)

    CALL MOVE(1.1,17.0)
    CALL TEXT(' FREQ=',6)
    CALL NUMBD(RF,'(F8.3)',8)
    CALL TEXT('     BB=',8)
    CALL NUMBD(BB,'(F8.3)',8)
    CALL TEXT('     RR=',8)
    CALL NUMBD(RR,'(F8.3)',8)
    CALL TEXT('     RZ=',8)
    CALL NUMBD(RZ,'(F8.3)',8)
    CALL TEXT('  WALLR=',8)
    CALL NUMBD(WALLR,'(1PE8.1)',8)

    CALL MOVE(1.1,16.5)
    CALL TEXT('  NXP=',6)
    CALL NUMBI(NXPMAX,'(I8)',8)
    CALL TEXT('    NZP=',8)
    CALL NUMBI(NZMAX,'(I8)',8)
    CALL TEXT('     RA=',8)
    CALL NUMBD(RA,'(F8.3)',8)
    CALL TEXT('     RD=',8)
    CALL NUMBD(RD,'(F8.3)',8)
    CALL TEXT('     RB=',8)
    CALL NUMBD(RB,'(F8.3)',8)

    CALL MOVE(1.1,16.0)
    CALL TEXT('              ',14)
    CALL TEXT('                ',16)
    CALL TEXT('        ',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL TEXT('    NS-',7)
       CALL NUMBI(NS,'(I1)',1)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,15.5)
    CALL TEXT('  J-H=',6)
    CALL NUMBD(AJYH(1),'(F8.3)',8)
    CALL TEXT('    J-L=',8)
    CALL NUMBD(AJYL(1),'(F8.3)',8)
    CALL TEXT('  N-20 :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PN(NS),'(F8.3)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,15.0)
    CALL TEXT('  R-H=',6)
    CALL NUMBD(RANT1(1),'(F8.3)',8)
    CALL TEXT('    R-L=',8)
    CALL NUMBD(RANT2(1),'(F8.3)',8)
    CALL TEXT('  T-PR :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PTPR(NS),'(F8.3)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,14.5)
    CALL TEXT('  X-H=',6)
    CALL NUMBD(XANT1(1),'(F8.3)',8)
    CALL TEXT('    X-L=',8)
    CALL NUMBD(XANT2(1),'(F8.3)',8)
    CALL TEXT('  T-PP :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PTPP(NS),'(F8.3)',8)
    END DO

    IF(PANT.EQ.0.D0) THEN
       RNORM2=1.D0
    ELSE
       RNORM2=1.D0/PANT
    END IF

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,14.0)
    CALL TEXT('PANTH=',6)
    CALL NUMBD(PANT1*RNORM2,'(F8.3)',8)
    CALL TEXT('  PANTL=',8)
    CALL NUMBD(PANT2*RNORM2,'(F8.3)',8)
    CALL TEXT('  N-S  :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PNS(NS),'(F8.3)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,13.5)
    CALL TEXT('PWALH=',6)
    CALL NUMBD(PWALL1*RNORM2,'(F8.3)',8)
    CALL TEXT('  PWALL=',8)
    CALL NUMBD(PWALL2*RNORM2,'(F8.3)',8)
    CALL TEXT('  T-S  :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PTS(NS),'(F8.3)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,13.0)
    CALL TEXT('PABST=',6)
    CALL NUMBD(PABSTT*RNORM2,'(F8.3)',8)
    CALL TEXT('  AJCDT=',8)
    CALL NUMBD(AJCDT*RNORM2,'(F8.3)',8)
    CALL TEXT('  PABS :',8)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(PABSXZ(NS)*RNORM2,'(F8.3)',8)
    END DO

    DO NZ=1,NZMAX+1
       IF(NZ.LE.NZMAX/2) THEN
          NZZ = NZ + NZMAX/2
       ELSE
          NZZ = NZ - NZMAX/2
       ENDIF
       GZ(NZ)=AKZ(NZZ)
    END DO
    GZ(1    )=-AKZ(NZMAX/2+1)
    CALL GQSCAL(GZ(1),GZ(NZMAX+1),GZMIN,GZMAX,GZSCL)

    DO NZ=1,NZMAX
       GDATAZ(NZ)=REAL(CPANTK(NZ))*RNORM2
    END DO
    CALL GMNMX1(GDATAZ,1,NZMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,GSCAL)

    CALL SETLIN(0,0,4)
    IF(NGRAPH.EQ.3) THEN
       CALL GDEFIN(2.5,12.5,1.0,12.0,GZ(1),GZ(NZMAX+1),-0.1*GMAX,GMAX)
    ELSE
       CALL GDEFIN(2.5,12.5,1.0,6.5,GZ(1),GZ(NZMAX+1),-0.1*GMAX,GMAX)
    ENDIF
    CALL GFRAME
    CALL GSCALE(0.,GZSCL,0.,   GSCAL,0.1,9)
    CALL GSCALE(0.,0.,0.,2.*GMAX,0.0,0)
    CALL GVALUE(0.,0.,0.,4*GSCAL,NGULEN(4*GSCAL))
    CALL GVALUE(0.,2*GZSCL,0.,0.,     0)
    CALL SETLIN(0,0,7)
    DO NZ=1,NZMAX+1
       IF(NZ.LE.NZMAX/2) THEN
          NZZ = NZ + NZMAX/2
       ELSE
          NZZ = NZ - NZMAX/2
       ENDIF
       GDATAZ(NZ)=REAL(CPANTK(NZZ)*RNORM2)
    END DO
    CALL GPLOTP(GZ,GDATAZ,1,NZMAX+1,1,0,0,0)
    DO NZ=1,NZMAX+1
       GDATAZ(NZ)=0.0
    END DO
    DO NS=1,NSMAX
       DO NZ=1,NZMAX+1
          IF(NZ.LE.NZMAX/2) THEN
             NZZ = NZ + NZMAX/2
          ELSE
             NZZ = NZ - NZMAX/2
          ENDIF
          GDATAZ(NZ)=GDATAZ(NZ)+PAK(NZZ,NS)*RNORM2
       END DO
       CALL SETLIN(0,0,7-NS)
       CALL GPLOTP(GZ,GDATAZ,1,NZMAX+1,1,-NS,10,0)
    END DO

    IF(NGRAPH.NE.2.AND.NGRAPH.NE.3) THEN
       DO NZ=1,NZMAX+1
          IF(NZ.LE.NZMAX/2) THEN
             NZZ = NZ + NZMAX/2
          ELSE
             NZZ = NZ - NZMAX/2
          ENDIF
          GDATAZ(NZ)=AJCDK(NZZ)*RNORM2
       END DO
       CALL GMNMX1(GDATAZ,1,NZMAX,1,GAMIN,GAMAX)
       CALL SETLIN(0,0,4)
       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,GSCAL)
       CALL GDEFIN(2.5,12.5,6.5,12.0,GZ(1),GZ(NZMAX+1),GMIN,GMAX)
       CALL GFRAME
       CALL GSCALE(0.,  GZSCL,0.,  GSCAL,0.1,9)
       CALL GSCALE(0.,2*GZSCL,0.,2*GSCAL,0.2,9)
       CALL GSCALE(0.,0.,0.,2.*MAX(ABS(GMAX),ABS(GMIN)),0.0,0)
       CALL GVALUE(0.,0.,0.,4*GSCAL,NGULEN(4*GSCAL))
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GZ,GDATAZ,1,NZMAX+1,1,0,0,0)
    ENDIF

    DO NX=1,NXMAX
       NN=NXPRNT(NX)
       GX(NX)=XAM(NN)
       GXM(NX)=XAM(NX)
    END DO

    PGAMIN=0.0
    PGAMAX=0.0
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=PABSX(NX,NS)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
       PGAMIN=MIN(PGAMIN,GAMIN)
       PGAMAX=MAX(PGAMAX,GAMAX)
    END DO
    PGAMAX=MAX(ABS(PGAMAX),ABS(PGAMIN))
    CALL GQSCAL(0.,PGAMAX,PGMIN,PGMAX,PSCAL)

    IF(NGRAPH.NE.2.AND.NGRAPH.NE.3) THEN
       DO NX=1,NXPMAX
          GDATA1(NX)=AJCDX(NX)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
       IF(ABS(GAMAX).GE.ABS(GAMIN)) THEN
          SGAMAX=GAMAX
       ELSE
          SGAMAX=GAMIN
       ENDIF

       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,GSCAL)
       CALL GDEFIN(15.0,25.0,6.5,12.0,-GB,GB,GMIN,GMAX)
       CALL SETLIN(0,0,4)
       CALL GFRAME
       CALL GSCALE(0.,0.1,0.,0.,0.1,9)
       CALL GSCALE(0.,0.,0.,GSCAL,0.1,9)
       CALL GSCALE(0.,100.,0.,0.,0.2,9)
       CALL GSCALE(-GA,2.*GA,0.,20.*GSCAL,0.1,0)
       CALL GVALUE(0.,0.,0.,4*GSCAL,NGULEN(4*GSCAL))
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,0,0,0)
    ENDIF
    CALL SETLIN(0,0,4)

    IF(NGRAPH.EQ.3) THEN
       CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*PGMAX,PGMAX)
    ELSE
       CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*PGMAX,PGMAX)
    ENDIF
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,0.,0.1,9)
    CALL GSCALE(0.,0.,0.,PSCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,20.*PSCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*PSCAL,NGULEN(4*PSCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       DO NX=1,NXPMAX
          GDATA1(NX)=PABSX(NX,NS)*RNORM2
       END DO
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO

    IF(NMODEL.EQ.5.OR.NMODEL.EQ.3.OR.NMODEL.EQ.1) THEN
       DO NX=1,NXMAX
          NN        =NXPRNT(NX)
          GDATA1(NX)=FLUXX(NN)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
       IF(ABS(GAMIN).GE.ABS(GAMAX)) THEN
          DO NX=1,NXMAX
             GDATA1(NX)=-GDATA1(NX)
          END DO
          GAMAX=-GAMIN
       ENDIF
       CALL GQSCAL(0.,GAMAX*1.02,GMIN,GMAX,GSCAL)
       IF(NGRAPH.EQ.3) THEN
          CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*GMAX,GMAX)
       ELSE
          CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*GMAX,GMAX)
       ENDIF
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GX,GDATA1,1,NXMAX,1,0,0,0)
    ENDIF

    CALL SETLIN(0,0,7)
    CALL MOVE(1.1,12.5)
    CALL TEXT('PGAMAX=',6)
    CALL NUMBR(PGAMAX,'(F8.3)',8)
    CALL TEXT('  SGAMAX=',8)
    CALL NUMBR(SGAMAX,'(F8.3)',8)

    CALL PAGEE
    RETURN
  END SUBROUTINE W1GR2DR

  SUBROUTINE W1GPRX
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NS
    REAL(rkind):: RNORM2

    IF(PANT.EQ.0.D0) THEN
       RNORM2=1.D0
    ELSE
       RNORM2=1.D0/PANT
    END IF

    CALL SETLIN(0,0,7)
    CALL MOVE(14.,12.0)
    CALL TEXT('IR   =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(AHLT(NS,1)*RNORM2,'(F8.5)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14.,11.5)
    CALL TEXT('INR  =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(AHLT(NS,2)*RNORM2,'(F8.5)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14.,11.0)
    CALL TEXT('ITOT =',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD((AHLT(NS,1)+AHLT(NS,2))*RNORM2,'(F8.5)',8)
    END DO

    CALL SETLIN(0,0,7)
    CALL MOVE(14.,10.3)
    CALL TEXT('IR  T=',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(AHLT(NS,3)*RNORM2,'(F8.5)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14., 9.8)
    CALL TEXT('INR T=',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD(AHLT(NS,4)*RNORM2,'(F8.5)',8)
    END DO
    CALL SETLIN(0,0,7)
    CALL MOVE(14., 9.3)
    CALL TEXT('ITOTT=',6)
    DO NS=1,MIN(NSMAX,4)
       CALL SETLIN(0,0,7-NS)
       CALL NUMBD((AHLT(NS,3)+AHLT(NS,4))*RNORM2,'(F8.5)',8)
    END DO
    CALL SETLIN(0,2,4)

    RETURN
  END SUBROUTINE W1GPRX

  SUBROUTINE W1GR1DJ
    USE w1comm
    USE w1sub,ONLY: w1fftl
    USE libgrf
    IMPLICIT NONE
    INTEGER:: NZ
    INTEGER,DIMENSION(:),ALLOCATABLE:: NSHIFT
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA

    ALLOCATE(XDATA(NZMAX),FDATA(NZMAX,2),NSHIFT(NZMAX))
    DO NZ=1,NZMAX
       XDATA(NZ)=ZA(NZ)
    END DO

    SELECT CASE(MDLWG)
    CASE(1)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG1(NZ))
          FDATA(NZ,2)=AIMAG(CWG1(NZ))
       END DO
    CASE(2)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG2(NZ))
          FDATA(NZ,2)=AIMAG(CWG2(NZ))
       END DO
    CASE(3)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG3(NZ))
          FDATA(NZ,2)=AIMAG(CWG3(NZ))
       END DO
    CASE(4)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG4(NZ))
          FDATA(NZ,2)=AIMAG(CWG4(NZ))
       END DO
    END SELECT

    CALL W1FFTL(CWG1,NZMAX,0)
    CALL W1FFTL(CWG2,NZMAX,0)
    CALL W1FFTL(CWG3,NZMAX,0)
    CALL W1FFTL(CWG4,NZMAX,0)

    CALL PAGES
    CALL GRD1D(1,XDATA,FDATA,NZMAX,NZMAX,2,'@CWG1(Z)@')

    DO NZ=1,NZMAX/2
       NSHIFT(NZ        )=NZ+NZMAX/2
       NSHIFT(NZ+NZMAX/2)=NZ
    END DO

    DO NZ=1,NZMAX
       XDATA(NZ)=AKZ(NSHIFT(NZ))
    END DO

    SELECT CASE(MDLWG)
    CASE(1)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG1(NSHIFT(NZ)))
          FDATA(NZ,2)=AIMAG(CWG1(NSHIFT(NZ)))
       END DO
    CASE(2)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG2(NSHIFT(NZ)))
          FDATA(NZ,2)=AIMAG(CWG2(NSHIFT(NZ)))
       END DO
    CASE(3)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG3(NSHIFT(NZ)))
          FDATA(NZ,2)=AIMAG(CWG3(NSHIFT(NZ)))
       END DO
    CASE(4)
       DO NZ=1,NZMAX
          FDATA(NZ,1)= REAL(CWG4(NSHIFT(NZ)))
          FDATA(NZ,2)=AIMAG(CWG4(NSHIFT(NZ)))
       END DO
    END SELECT

    CALL GRD1D(2,XDATA,FDATA,NZMAX,NZMAX,2,'@CWG1(KZ)@')

    CALL PAGEE
    DEALLOCATE(XDATA,FDATA,NSHIFT)

    CALL W1FFTL(CWG1,NZMAX,1)
    CALL W1FFTL(CWG2,NZMAX,1)
    CALL W1FFTL(CWG3,NZMAX,1)
    CALL W1FFTL(CWG4,NZMAX,1)

  END SUBROUTINE W1GR1DJ

END MODULE w1grf1
