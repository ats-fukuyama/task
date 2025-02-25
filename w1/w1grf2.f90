MODULE w1grf2

CONTAINS

!     ****** 1-D GRAPHIC ROUTINE (B AND A) ******

  SUBROUTINE W1GR1B
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL

    IF(NGRAPH.EQ.0) RETURN

!    RNORM=1.D0/SQRT(PANT)
    RNORM=1.D0
    RNORM2=RNORM*RNORM

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,0,4)
    GA=RA
    GB=RB
    CALL MOVE( 2.5,17.7)
    CALL TEXT('BX,BY,BZ',8)
    CALL MOVE(15.0,17.7)
    CALL TEXT('AX,AZ,PHI',9)

    DO NX=1,NXPMAX
       GXM   (NX)=XAM(NX)
       GDATA1(NX)=  ABS(CBF(NX,1))*RNORM
       GDATA2(NX)= REAL(CBF(NX,1))*RNORM
       GDATA3(NX)=AIMAG(CBF(NX,1))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,12.,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=  ABS(CBF(NX,2))*RNORM
       GDATA2(NX)= REAL(CBF(NX,2))*RNORM
       GDATA3(NX)=AIMAG(CBF(NX,2))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,6.5,12.,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=  ABS(CBF(NX,3))*RNORM
       GDATA2(NX)= REAL(CBF(NX,3))*RNORM
       GDATA3(NX)=AIMAG(CBF(NX,3))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,1.,6.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=  ABS(CAF(NX,1))*RNORM
       GDATA2(NX)= REAL(CAF(NX,1))*RNORM
       GDATA3(NX)=AIMAG(CAF(NX,1))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,12.,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=  ABS(CAF(NX,3))*RNORM
       GDATA2(NX)= REAL(CAF(NX,3))*RNORM
       GDATA3(NX)=AIMAG(CAF(NX,3))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,6.5,12.,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=  ABS(CAF(NX,4))*RNORM
       GDATA2(NX)= REAL(CAF(NX,4))*RNORM
       GDATA3(NX)=AIMAG(CAF(NX,4))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,1.,6.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    CALL PAGEE
    RETURN
  END SUBROUTINE W1GR1B

  SUBROUTINE W1GRUF
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL
    REAL:: GAMIN1,GAMIN2,GAMIN3,GAMAX1,GAMAX2,GAMAX3

    IF(NGRAPH.EQ.0) RETURN

!    RNORM=1.D0/SQRT(PANT)
    RNORM=1.D0
    RNORM2=RNORM*RNORM

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,2,4)
    GA=RA
    GB=RB
    CALL MOVE( 2.5,17.7)
    CALL TEXT('QH,E.B,EPARA,LAMBDA',19)
    CALL MOVE(15.0,17.7)
    CALL TEXT('FR,FNR,JE,JET',13)

    DO NX=1,NXPMAX
       GXM   (NX)=XAM(NX)
       GDATA1(NX)=RHL(NX,2)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    GAMAX=MAX(GAMAX,ABS(GAMIN))
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,13.5,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,2,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=RHL(NX,3)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,5,NXPMAX-4,1,GAMIN,GAMAX)
    GAMAX=MAX(GAMAX,ABS(GAMIN))
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,9.5,13.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
!     CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,2,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)

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
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,2,7)
!     CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)
    CALL GVALUE(0.,REAL(RA),0.,0.,2)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,1)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,13.5,17.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,1)*RNORM2
       END DO
       CALL SETLIN(0,2,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,2,4)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,2)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,5,NXPMAX-4,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,9.5,13.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
!     CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,2)*RNORM2
       END DO
       CALL SETLIN(0,2,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,2,4)


    DO NX=1,NXPMAX
       GDATA1(NX)=RHL(NX,6)*RNORM2
       GDATA2(NX)=RHL(NX,7)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,5,NXPMAX-4,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,5,NXPMAX-4,1,GAMIN2,GAMAX2)
    GAMAX=MAX(GAMAX1,GAMAX2)
    GAMIN=MIN(GAMIN1,GAMIN2)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,5.5,9.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,2,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)

    DO NX=1,NXPMAX
       IF(RHL(NX,6).EQ.0.D0) THEN
          GDATA1(NX)=0.0
       ELSE
          GDATA1(NX)=-RHL(NX,7)/RHL(NX,6)
       END IF
!        IF(F.GT.1.D0) THEN  
!            GDATA1(NX)= 2.D0-1.D0/F
!         ELSEIF(F.LT.-1.D0) THEN
!           GDATA1(NX)=-2.D0-1.D0/F
!        ELSE
!           GDATA1(NX)= F
!        ENDIF
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    SCAL=0.25
    CALL GDEFIN(2.5,12.5,1.5,5.5,-GB,GB,-0.5,1.5)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL SETLIN(1,2,4)
    CALL GSCALE(0.,0.,1.,13.*SCAL,0.1,0)
    CALL SETLIN(0,2,4)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,2,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=(AHL(NX,1,1)+AHL(NX,1,2))*RNORM2
       GDATA2(NX)=AHL(NX,1,1)*RNORM2
       GDATA3(NX)=AHL(NX,1,2)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,1,NXPMAX,1,GAMIN2,GAMAX2)
    CALL GMNMX1(GDATA3,1,NXPMAX,1,GAMIN3,GAMAX3)
    GAMIN=MIN(GAMIN1,GAMIN2,GAMIN3)
    GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,5.5,9.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.2,0.,2*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
!     CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,2,7)
!     CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,2,4)


    CALL PAGEE
    RETURN
  END SUBROUTINE W1GRUF

!     ****** 1-D GRAPHIC ROUTINE (HELICITY) ******

  SUBROUTINE W1GR1H
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL
    REAL:: GAMIN1,GAMIN2,GAMIN3,GAMAX1,GAMAX2,GAMAX3

    IF(NGRAPH.EQ.0) RETURN

!    RNORM=1.D0/SQRT(PANT)
    RNORM=1.D0
    RNORM2=RNORM*RNORM

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,0,4)
    GA=RA
    GB=RB
    CALL MOVE( 2.5,17.7)
    CALL TEXT('FTOT,HEL,RAT',11)
    CALL MOVE(15.0,17.7)
    CALL TEXT('JE,JR,JNR',9)

    DO NX=1,NXPMAX
       GXM   (NX)=XAM(NX)
       GDATA1(NX)=RHL(NX,4)*RNORM2
       GDATA2(NX)=RHL(NX,5)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,1,NXPMAX,1,GAMIN2,GAMAX2)
    GAMAX=MAX(GAMAX1,GAMAX2)
    GAMIN=MIN(GAMIN1,GAMIN2)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,12.,17.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1, -3, 10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=RHL(NX,6)*RNORM2
       GDATA2(NX)=RHL(NX,7)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,2,NXPMAX-1,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,2,NXPMAX-1,1,GAMIN2,GAMAX2)
    GAMAX=MAX(GAMAX1,GAMAX2)
    GAMIN=MIN(GAMIN1,GAMIN2)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,6.5,12.,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1, -3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       IF(RHL(NX,6).EQ.0.D0) THEN
          GDATA1(NX)=0.0
       ELSE
          GDATA1(NX)=-RHL(NX,7)/RHL(NX,6)
       END IF
!        IF(F.GT.1.D0) THEN  
!            GDATA1(NX)= 2.D0-1.D0/F
!         ELSEIF(F.LT.-1.D0) THEN
!           GDATA1(NX)=-2.D0-1.D0/F
!        ELSE
!           GDATA1(NX)= F
!        ENDIF
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    SCAL=0.25
    CALL GDEFIN(2.5,12.5,1.,6.5,-GB,GB,-1.999,1.999)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=(AHL(NX,1,1)+AHL(NX,1,2))*RNORM2
       GDATA2(NX)=AHL(NX,1,1)*RNORM2
       GDATA3(NX)=AHL(NX,1,2)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,1,NXPMAX,1,GAMIN2,GAMAX2)
    CALL GMNMX1(GDATA3,1,NXPMAX,1,GAMIN3,GAMAX3)
    GAMIN=MIN(GAMIN1,GAMIN2,GAMIN3)
    GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,12.,17.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=AHL(NX,NS,1)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,6.5,12.,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=AHL(NX,NS,1)*RNORM2
       END DO
       CALL SETLIN(0,0,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,0,4)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=AHL(NX,NS,2)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,1.,6.5,-GB,GB,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=AHL(NX,NS,2)*RNORM2
       END DO
       CALL SETLIN(0,0,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,0,4)

    CALL PAGEE
    RETURN
  END SUBROUTINE W1GR1H

!  

!     ****** 1-D GRAPHIC ROUTINE ******

  SUBROUTINE W1GR1D
    USE w1comm
    USE w1grf1,ONLY: w1gprm
    IMPLICIT NONE
    INTEGER:: NZ,NX,NN,NS,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL
    REAL:: PGAMIN,PGAMAX,PGMIN,PGMAX,SGAMAX,PSCAL

    IF(NGRAPH.EQ.0) RETURN

    IF(ABS(PANT).LE.0.D0) THEN
       RNORM=1.D0
    ELSE
       RNORM=1.D0/SQRT(ABS(PANT))
    ENDIF
    RNORM2=RNORM*RNORM

1   IF(NZMAX.EQ.1) THEN
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
    CALL SETLIN(0,0,4)
    CALL MOVE( 2.5,17.7)
    CALL TEXT('EX,EY,EZ',8)

    GA=RA
    GB=RB
    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GX    (NX)=XAM(NN)
       GXM   (NX)=XAM(NX)
       GDATA1(NX)=  ABS(CE2DA(NZ,NN,1))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,1))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,1))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,12.,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GX,GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GX,GDATA2,1,NXMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GX,GDATA3,1,NXMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GDATA1(NX)=  ABS(CE2DA(NZ,NN,2))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,2))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,2))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,6.5,12.,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GX,GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GX,GDATA2,1,NXMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GX,GDATA3,1,NXMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXMAX
       NN        =NXPRNT(NX)
       GDATA1(NX)=  ABS(CE2DA(NZ,NN,3))*RNORM
       GDATA2(NX)= REAL(CE2DA(NZ,NN,3))*RNORM
       GDATA3(NX)=AIMAG(CE2DA(NZ,NN,3))*RNORM
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,1.,6.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GX,GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GX,GDATA2,1,NXMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GX,GDATA3,1,NXMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    SELECT CASE(NMODEL)
    CASE(0:5)

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

    DO NX=1,NXPMAX
       GDATA1(NX)=AJCDX(NX)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    IF(ABS(GAMAX).GE.ABS(GAMIN)) THEN
       SGAMAX=GAMAX
    ELSE
       SGAMAX=GAMIN
    ENDIF
      
    IF(MOD(NGRAPH,4).NE.2.AND.MOD(NGRAPH,4).NE.3) THEN
       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
       CALL GDEFIN(15.0,25.0,6.5,12.0,-GB,GB,GMIN,GMAX)
       CALL GFRAME
       CALL GSCALE(0.,0.1,0.,0.,0.1,9)
       CALL GSCALE(0.,0.,0.,SCAL,0.1,9)
       CALL GSCALE(0.,100.,0.,0.,0.2,9)
       CALL GSCALE(-GA,2.*GA,0.,20.*SCAL,0.1,0)
       CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,0,0,0)
       CALL SETLIN(0,0,4)
    ENDIF

    IF(MOD(NGRAPH,4).EQ.3) THEN
       CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*PGMAX,PGMAX)
    ELSE
       CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*PGMAX,PGMAX)
    ENDIF
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,0.,0.1,9)
    CALL GSCALE(0.,0.,0.,PSCAL,0.1,1)
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
    CALL SETLIN(0,0,4)

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
       IF(MOD(NGRAPH,4).EQ.3) THEN
          CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*GMAX,GMAX)
       ELSE
          CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*GMAX,GMAX)
       ENDIF
       CALL GSCALE(0.,0.,0.,SCAL,0.1,5)
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GX,GDATA1,1,NXMAX,1,0,0,0)
    ENDIF

    CASE(6:11)

    PGAMIN=0.0
    PGAMAX=0.0
    DO NS=1,NSMAX
       DO NX=1,NXMAX
          GDATA1(NX)=PABSX(NX,NS)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
       PGAMIN=MIN(PGAMIN,GAMIN)
       PGAMAX=MAX(PGAMAX,GAMAX)
    END DO
    PGAMAX=MAX(ABS(PGAMAX),ABS(PGAMIN))
    CALL GQSCAL(0.,PGAMAX,PGMIN,PGMAX,PSCAL)

    DO NX=1,NXMAX
       GDATA1(NX)=AJCDX(NX)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    IF(ABS(GAMAX).GE.ABS(GAMIN)) THEN
       SGAMAX=GAMAX
    ELSE
       SGAMAX=GAMIN
    ENDIF
      
    IF(MOD(NGRAPH,4).NE.2.AND.MOD(NGRAPH,4).NE.3) THEN
       CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
       CALL GDEFIN(15.0,25.0,6.5,12.0,-GB,GB,GMIN,GMAX)
       CALL GFRAME
       CALL GSCALE(0.,0.1,0.,0.,0.1,9)
       CALL GSCALE(0.,0.,0.,SCAL,0.1,9)
       CALL GSCALE(0.,100.,0.,0.,0.2,9)
       CALL GSCALE(-GA,2.*GA,0.,20.*SCAL,0.1,0)
       CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GXM,GDATA1,1,NXMAX,1,0,0,0)
       CALL SETLIN(0,0,4)
    ENDIF

    IF(MOD(NGRAPH,4).EQ.3) THEN
       CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*PGMAX,PGMAX)
    ELSE
       CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*PGMAX,PGMAX)
    ENDIF
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,0.,0.1,9)
    CALL GSCALE(0.,0.,0.,PSCAL,0.1,1)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,20.*PSCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*PSCAL,NGULEN(4*PSCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       CALL SETLIN(0,0,7-NS)
       DO NX=1,NXMAX
          GDATA1(NX)=PABSX(NX,NS)*RNORM2
       END DO
       CALL GPLOTP(GXM,GDATA1,1,NXMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,0,4)

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
       IF(MOD(NGRAPH,4).EQ.3) THEN
          CALL GDEFIN(15.0,25.0,1.0,12.0,-GB,GB,-0.1*GMAX,GMAX)
       ELSE
          CALL GDEFIN(15.0,25.0,1.0,6.5,-GB,GB,-0.1*GMAX,GMAX)
       ENDIF
       CALL GSCALE(0.,0.,0.,SCAL,0.1,5)
       CALL SETLIN(0,0,7)
       CALL GPLOTP(GX,GDATA1,1,NXMAX,1,0,0,0)
    ENDIF
    END SELECT

    CALL W1GPRM(NZ)

    CALL PAGEE
    IF(NZMAX.GT.1) GOTO 1
    RETURN

9990 CONTINUE
    RETURN
  END SUBROUTINE W1GR1D

!     ****** 1-D GRAPHIC ROUTINE (HELICITY) ******

  SUBROUTINE W1GR1F
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS,NGULEN
    REAL(rkind):: RNORM,RNORM2
    REAL:: GA,GB,GAMIN,GAMAX,GMIN,GMAX,SCAL
    REAL:: GAMIN1,GAMIN2,GAMIN3,GAMAX1,GAMAX2,GAMAX3

    IF(NGRAPH.EQ.0) RETURN

!    RNORM=1.D0/SQRT(PANT)
    RNORM=1.D0
    RNORM2=RNORM*RNORM

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,0,4)
    GA=RA
    GB=RB
    CALL MOVE( 2.5,17.7)
    CALL TEXT('A.B,QH,E.B',10)
    CALL MOVE(15.0,17.7)
    CALL TEXT('FE,FR,FNR',9)

    DO NX=1,NXPMAX
       GXM   (NX)=XAM(NX)
       GDATA1(NX)=RHL(NX,1)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    GAMAX=MAX(GAMAX,ABS(GAMIN))
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,12.,17.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=RHL(NX,2)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN,GAMAX)
    GAMAX=MAX(GAMAX,ABS(GAMIN))
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,6.5,12.,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=RHL(NX,3)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,2,NXPMAX-1,1,GAMIN,GAMAX)
    GAMAX=MAX(GAMAX,ABS(GAMIN))
    CALL GQSCAL(-GAMAX,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(2.5,12.5,1.,6.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,4)

    DO NX=1,NXPMAX
       GDATA1(NX)=(FHL(NX,1,1)+FHL(NX,1,2))*RNORM2
       GDATA2(NX)=FHL(NX,1,1)*RNORM2
       GDATA3(NX)=FHL(NX,1,2)*RNORM2
    END DO
    CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,1,NXPMAX,1,GAMIN2,GAMAX2)
    CALL GMNMX1(GDATA3,1,NXPMAX,1,GAMIN3,GAMAX3)
    GAMIN=MIN(GAMIN1,GAMIN2,GAMIN3)
    GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3,ABS(GAMIN))
    CALL GQSCAL(0.0,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,12.,17.5,-GB,GB,-0.1*GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,0,7)
    CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1, 0, 0,0)
    CALL SETLIN(0,0,6)
    CALL GPLOTP(GXM,GDATA2,1,NXPMAX,1,-7,10,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(GXM,GDATA3,1,NXPMAX,1,-3,10,0)
    CALL SETLIN(0,0,4)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,1)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(0.0,MAX(GAMAX,ABS(GAMIN)),GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,6.5,12.,-GB,GB,-0.1*GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,1)*RNORM2
       END DO
       CALL SETLIN(0,0,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,0,4)

    GAMIN= 1.E32
    GAMAX=-1.E32
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,2)*RNORM2
       END DO
       CALL GMNMX1(GDATA1,1,NXPMAX,1,GAMIN1,GAMAX1)
       GAMIN=MIN(GAMIN,GAMIN1)
       GAMAX=MAX(GAMAX,GAMAX1)
    END DO
    CALL GQSCAL(0.0,MAX(GAMAX,ABS(GAMIN)),GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,1.,6.5,-GB,GB,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.1,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,100.,0.,0.,0.2,9)
    CALL GSCALE(-GA,2.*GA,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,REAL(RA),0.,0.,2)
    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          GDATA1(NX)=FHL(NX,NS,2)*RNORM2
       END DO
       CALL SETLIN(0,0,7-NS)
       CALL GPLOTP(GXM,GDATA1,1,NXPMAX,1,-NS,10,0)
    END DO
    CALL SETLIN(0,0,4)

    CALL PAGEE
    RETURN
  END SUBROUTINE W1GR1F
END MODULE w1grf2
