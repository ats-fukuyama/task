! wfhout.f90

MODULE wfhout

  INTEGER,PARAMETER:: dp = KIND(1.0D0)
  
  PRIVATE
  PUBLIC wf_hout

CONTAINS

! New graphic module

  SUBROUTINE wf_hout

    INCLUDE 'wfcomm.inc'
    INTEGER,PARAMETER:: NCHMAX=80
    CHARACTER(LEN=NCHMAX):: KLINE,KWORD(12),KWD,KWTEMP
    CHARACTER(LEN=1):: KID,KG1,KG2,KG3,KG4
    INTEGER:: NL,NWD,NCH
    REAL(dp),save:: xdmin=0.D0,xdmax=0.D0,ydmin=0.D0,ydmax=0.D0
    INTEGER,SAVE:: nconf=1,n1max=100,n2max=100
    REAL(dp),ALLOCATABLE:: xg(:),yg(:)
    INTEGER,ALLOCATABLE:: nelmxy(:,:)

! --- initializaion ---

    INIT=0
    IREGION=0
    IPAGE=0

! --- set boundary of x and y ---

    IF(xdmin.EQ.xdmax) THEN
       xdmin=MIN(xd(1:nnod))+1.D-6
       xdmax=MAX(xd(1:nnod))-1.D-6
    END IF
    IF(ydmin.EQ.ydmax) THEN
       ydmin=MIN(yd(1:nnod))+1.D-6
       ydmax=MAX(yd(1:nnod))-1.D-6
    END IF

! --- input graphic commands ---

1   WRITE(6,'(A)') '# INPUT : A,X/Y/Z/F,F/R,R/I'
    WRITE(6,'(A)') '          R=region X=exit'
    READ(5,'(A80)',ERR=1,END=9000) KLINE

! --- separate KLINE int KWORD ---

9   NL=0
    NWD=0
    NCH=0
10  IF(NL.GE.80) GOTO 20
    NL=NL+1
    KID=KLINE(NL:NL)
    CALL GUCPTL(KID)
    IF(KID.NE.','.AND.KID.NE.' ') THEN
       IF(NCH.LT.NCHMAX) NCH=NCH+1
       KWD(NCH:NCH)=KID
    ELSE
       IF(NCH.NE.0) THEN
          IF(NWD.LT.NWDM) NWD=NWD+1
          KWORD(NWD)=KWD(1:NCH)
          NCH=0
       ENDIF
    ENDIF
    GOTO 10

20  NWMAX=NWD

! --- set first KWORD ---

    NWD=0

! --- set first KWORD and extract KG1..4  ---

21  NWD=NWD+1
    IF(NWD.GT.NWMAX) GOTO 30
    KWD=KWORD(NWD)
    KG1=KWD(1:1)
    KG2=KWD(2:2)
    KG3=KWD(3:3)
    KG4=KWD(4:4)

! --- KG1='X' for exit  ---

    IF(KG1.EQ.'X') GO TO 9000

! --- KG1='R' for setting plottin region  ---

    IF(KG1.EQ.'R'.OR.INIT.EQ.0) THEN
101    WRITE(6,'(A)') '# INPUT: nconf,xdmin,xdmax,ydmin,ydmax,n1max,n2max'
       WRITE(6,'(A)') '           nconf=1,4:XY, 2,5:XZ, 3,6:YZ, 0:end'
       WRITE(6,'(A,I3,1P4E12.4,2I5)') &
            '       ',nconf,xdmin,xdmax,ydmin,ydmax,n1max,n2max
       nconf_temp=nconf
       READ(5,*,ERR=102,END=109) nconf_temp,xdmin,xdmax,ydmin,ydmax,nxmax,nxmax
       IF(nconf_temp.EQ.0) GOTO 109
       nconf=nconf_temp

! --- setup nelmxy table for quick local data access ---

       SELECT CASE(nconf)
       CASE(1)
          IF(ALLOCATED(xg)) DEALLOCATE(xg)
          ALLOCATE(xg(n1max))
          delx=(xdmax-xdmin)/(nxmax-1)
          DO nx=1,nxmax
             xg(nx)=xdmin+delx*(nx-1)
          END DO
          IF(ALLOCATED(yg)) DEALLOCATE(yg)
          ALLOCATE(yg(n2max))
          dely=(ydmax-ydmin)/(n2max-1)
          DO ny=1,n2max
             yg(nx)=ydmin+delx*(ny-1)
          END DO
          DEALLOCATE(nelmxy)
          ALLOCATE(nelmxy(n1max,n2max))
          nelm=1
          DO ny=1,n2max
             DO nx=1,n1max
                CALL WFFEP(xg(nx),yg(ny),nelmxy(nx,ny))
             END DO
          END DO
       END SELECT
       INIT=1
       IREGION=1
       GO TO 1
    END IF

! --- setup field data CA and its derivatives for nx and ny' ---

    IF(IREGION.EQ.1) THEN
       DO ny=1,n2max
          DO nx=1,n1max
             CALL FIELDA(xg(nx),yg(ny),nelmxy(nx,ny),NZ,0, &
                         CAFA(1:4,nx,ny),CAFXA(1:4,nx,ny), &
                         CAFYA(1:4,nx,ny),CAFZA(1:4,nx,ny))
          END DO
       END DO
       DO ny=1,n2max
          DO nx=1,n1max
             CALL FIELDA(nelmxy(nx,ny),NZ,1, &
                  CARA(1:4,nx,ny),CARXA(1:4,nx,ny), &
                  CARYA(1:4,nx,ny),CARZA(1:4,nx,ny))
          END DO
       END DO
       IREGION=2
       GO TO 1
    END IF

! --- Extract field data from CA for (x,y)' ---

    SELECT CASE(KG1)
    CASE('A')
       SELECT CASE(KG2)
       CASE('X')
          ID2=1
       CASE('Y')
          ID2=2
       CASE('Z')
          ID2=3
       CASE('F')
          ID2=4
       CASE DEFAULT
          WRITE(6,*) 'XX Undefined KG2:',KG2
          GO TO 1
       END SELECT ! KG2

       SELECT CASE(KG3)
       CASE('F')
          SELECT CASE(KG4)
          CASE('R')
             DO ny=1,n2max
                DO nx=1,n1max
                   FG(nx,ny)=REAL(CAFA(ID2,nx,ny))
                END DO
             END DO
          CASE('I')
             DO ny=1,n2max
                DO nx=1,n1max
                   FG(nx,ny)=AIMAG(CAFA(ID2,nx,ny))
                END DO
             END DO
          END SELECT
       CASE('R')
          SELECT CASE(KG4)
          CASE('R')
             DO ny=1,n2max
                DO nx=1,n1max
                   FG(nx,ny)=REAL(CARA(ID2,nx,ny))
                END DO
             END DO
          CASE('I')
             DO ny=1,n2max
                DO nx=1,n1max
                   FG(nx,ny)=AIMAG(CARA(ID2,nx,ny))
                END DO
             END DO
          END SELECT ! KG4
       END SELECT ! KG3
    END SELECT ! KG1

    
    
