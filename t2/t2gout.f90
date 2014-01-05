!
!
!
Module T2GOUT

  USE T2COMM, ONLY: &
         i0ikind,i0rkind

  PRIVATE
  PUBLIC T2_GOUT
 
  REAL(i0rkind),DIMENSION(:),ALLOCATABLE:: rhonrho,chinchi
  INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE:: nlnrho,nnnrho
  INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE:: nrhonl,nnnl
  REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: chinl
  INTEGER(i0ikind):: nlmax,nnmax

CONTAINS

  SUBROUTINE T2_GOUT

    USE T2CNST,ONLY: i0ikind,i0rkind
    USE T2PARM,ONLY: T2_PARM
    USE T2COMM
    USE libgrf
    IMPLICIT NONE
    INTEGER(i0ikind)    :: ierr,mode,ind
    CHARACTER(LEN=80) :: line,kw
    CHARACTER(LEN=1) :: kid,kch
    INTEGER(i0ikind) :: nwmax,iloc0,nw,i,ich0,nch,ich
    CHARACTER(LEN=80),DIMENSION(40):: kword,kwid,knum
    INTEGER(i0ikind),DIMENSION(40):: inum

    CALL T2_GSETUP

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '#### T2 GOUT: Rn/radial An/average Cn,Pn/contour Bn/birdseye X/exit'
    CALL TASK_KLIN(line,kid,mode,T2_PARM)
    IF(mode == 2 .OR. mode == 3) GOTO 1

!   --- separate words in a line ---
    ILOC0=1
    NWMAX=0
    DO I=1,80
       IF(ILOC0==0) THEN
          IF(LINE(I:I)/=' ' .AND. LINE(I:I)/=',') THEN
             IF(NWMAX/=0) THEN
                ILOC0=I
             END IF
          END IF
       ELSE
          IF(LINE(I:I)==' ' .OR. LINE(I:I)==',') THEN
             NWMAX=NWMAX+1
             KWORD(NWMAX)=LINE(ILOC0:I)
             ILOC0=0
          END IF
       END IF
    END DO
    IF(ILOC0 /= 0) THEN
       NWMAX=NWMAX+1
       KWORD(NWMAX)=LINE(ILOC0:I)
    END IF

!    DO NW=1,NWMAX
!       WRITE(6,'(A,I5,4X,A)') 'NW,KWORD=',NW,TRIM(KWORD(NW))
!    END DO

!   --- separate id and number in a word ---

    DO NW=1,NWMAX
       KW=KWORD(NW)
       ICH0=0
       DO NCH=1,LEN(KW)
!         --- lower-case char to upper-case char ---
          ICH=ICHAR(KW(NCH:NCH))
          IF(ICH.GE.97.AND.ICH.LE.122) ICH=ICH-32
          KCH=CHAR(ICH)
          KW(NCH:NCH)=KCH

          IF(KCH.GE.'A'.AND.KCH.LE.'Z') THEN
             ICH0=ICH0+1
          ELSE
             EXIT
          ENDIF
       END DO
       IF(ICH0.EQ.0) THEN
          KWID(NW)=' '
          KNUM(NW)=KW
       ELSE
          KWID(NW)=KW(1:ICH0)
          KNUM(NW)=KW(ICH0+1:LEN(KW))
       END IF
       INUM(NW)=0
       READ(KNUM(NW),'(I)',ERR=901,END=901) INUM(NW)
901    CONTINUE

    END DO

    DO NW=1,NWMAX
!       WRITE(6,'(A,A,A,A,A,I5)') &
!            'KID,KNUM=:',TRIM(KWID(NW)),':',TRIM(KNUM(NW)),':',INUM(NW)
       SELECT CASE(LEN_TRIM(KWID(NW)))
       CASE(1)
          SELECT CASE(KWID(NW))
          CASE('R')
             CALL T2_GR(INUM(NW))
          CASE('C')
             CALL T2_GC(INUM(NW),1)
          CASE('D')
             CALL T2_GC(INUM(NW),2)
          CASE('P')
             CALL T2_GC(INUM(NW),3)
          CASE('B')
             CALL T2_GC(INUM(NW),14)
          CASE('Y')
             CALL T2_GC(INUM(NW),12)
          CASE('X')
             GO TO 9000
          END SELECT
       END SELECT
    END DO

    GO TO 1

9000 CONTINUE
    CALL T2_GRELEASE

    RETURN

  END SUBROUTINE T2_GOUT

  SUBROUTINE T2_GSETUP
    USE T2COMM, ONLY: &
         i0ikind,i0rkind,twopi,i0xmax,d1guv,i0vmax, &
         i0lmax,i0pdiv_number,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax
    IMPLICIT NONE
    INTEGER(i0ikind):: nchi,nl,nrho,nchimaxl,nr,ierr
    REAL(i0rkind):: dchi,drho

    nlmax=i0lmax
    nchimax=i0pdiv_number*2**(i1mlvl(nlmax)-1)
    ALLOCATE(nrhonl(0:nlmax),nnnl(0:nlmax))
    ALLOCATE(chinchi(nchimax+1),chinl(nchimax+1,nlmax))

    dchi=twopi/nchimax
    DO nchi=1,nchimax+1
       chinchi(nchi)=dchi*(nchi-1)
    END DO

    nrhonl(0)=1
    nnnl(0)=1
    nrhomax=2
    nnmax=2
    DO nl=1,nlmax
       nrhonl(nl)=nrhomax
       nnnl(nl)=nnmax
       nchimaxl=i0pdiv_number*2**(i1mlvl(nl)-1)
       nrhomax=nrhomax+i1rdn2(nl)
       nnmax=nnmax+nchimaxl*i1rdn2(nl)
    END DO
    nrhomax=nrhomax-1
    nnmax=nnmax-1
    AlLOCATE(nlnrho(nrhomax),rhonrho(nrhomax),nnnrho(nrhomax))

    nrho=1
    nlnrho(nrho)=0
    rhonrho(nrho)=0.D0
    nnnrho(nrho)=1
    DO nl=1,nlmax
       drho=(d1rec(nl)-d1rec(nl-1))/i1rdn2(nl)
       nchimaxl=i0pdiv_number*2**(i1mlvl(nl)-1)
       dchi=twopi/nchimaxl
       DO nchi=1,nchimaxl+1
          chinl(nchi,nl)=dchi*(nchi-1)
       END DO
       DO nr=1,i1rdn2(nl)
          nrho=nrhonl(nl)+(nr-1)
          nlnrho(nrho)=nl
          rhonrho(nrho)=d1rec(nl-1)+drho*nr
          nnnrho(nrho)=nnnl(nl)+nchimaxl*(nr-1)
!          write(6,'(A,4I5,1PE12.4)') &
!               'nl,nr,nrho,nnnrho(nrho),rhonrho(nrho)=', &
!               nl,nr,nrho,nnnrho(nrho),rhonrho(nrho)
       END DO
    END DO

    RETURN
  END SUBROUTINE T2_GSETUP

  SUBROUTINE T2_GRELEASE

    DEALLOCATE(nrhonl,nnnl)
    DEALLOCATE(chinchi,chinl)
    DEAlLOCATE(nlnrho,rhonrho,nnnrho)

    RETURN
  END SUBROUTINE T2_GRELEASE

  SUBROUTINE T2_GC(INUM,ID)
    USE libgrf,ONLY: GRD2D
    USE T2COMM, ONLY: &
         i0ikind,i0rkind,twopi,i0xmax,d1guv,i0vmax, &
         i0lmax,i0pdiv_number,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax,d0rw
    IMPLICIT NONE
    INTEGER,PARAMETER:: nxmax=41,nymax=41
    INTEGER(i0ikind),INTENT(IN):: inum,id
    REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: gz
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE:: gzl,dgzl,chig
    REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: ugzl
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE:: gx,gy
    REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: ddx,ddy,ddxy,gxy
    REAL(i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: ugz
    INTEGER(i0ikind):: nchi,nl,nrho,nchimaxl,nr,ierr,nchig,nx,ny
    REAL(i0rkind):: dchig,xmin,xmax,ymin,ymax,dx,dy,x,y,r,th
    CHARACTER(LEN=80):: LINE

    nchig=MAX(nchimax,72)
    ALLOCATE(gz(nrhomax,nchig+1),chig(nchig+1))
    ALLOCATE(gzl(nchimax+1),dgzl(nchimax+1),ugzl(4,nchimax+1))
    dchig=TWOPI/nchig
    DO nchi=1,nchig+1
       chig(nchi)=dchig*(nchi-1)
       gz(1,nchi)=d1guv(inum)
    END DO
    DO nrho=2,nrhomax
       nl=nlnrho(nrho)
       nchimaxl=i0pdiv_number*2**(i1mlvl(nl)-1)
       DO nchi=1,nchimaxl
          gzl(nchi)=d1guv(i0vmax*(nnnrho(nrho)+nchi-2)+inum)
       END DO
       gzl(nchimaxl+1)=d1guv(i0vmax*(nnnrho(nrho)-1)+inum)
       CALL SPL1D(chinl(1:nchimaxl+1,nl),gzl,dgzl,ugzl,nchimaxl+1,4,ierr)
       DO nchi=1,nchig+1
          CALL SPL1DF(chig(nchi),gz(nrho,nchi), &
                      chinl(1:nchimaxl+1,nl),ugzl,nchimaxl+1,ierr)
       END DO
    END DO

    IF(ID.GT.10) THEN
       ALLOCATE(ddx(nrhomax,nchig),ddy(nrhomax,nchig),ddxy(nrhomax,nchig))
       ALLOCATE(ugz(4,4,nrhomax,nchig))
       ALLOCATE(gx(nxmax),gy(nymax),gxy(nxmax,nymax))

       CALL SPL2D(rhonrho,chig,gz,ddx,ddy,ddxy,ugz, &
                  nrhomax,nrhomax,nchig+1,0,0,ierr)
       xmin=-rhonrho(nrhomax)
       xmax= rhonrho(nrhomax)
       ymin=-rhonrho(nrhomax)
       ymax= rhonrho(nrhomax)
       dx=(xmax-xmin)/(nxmax-1)
       dy=(ymax-ymin)/(nymax-1)
       DO nx=1,nxmax
          gx(nx)=xmin+dx*(nx-1)
       END DO
       DO ny=1,nymax
          gy(ny)=ymin+dy*(ny-1)
       END DO
       DO nx=1,nxmax
          x=gx(nx)
          DO ny=1,nymax
             y=gy(ny)
             r=SQRT(x*x+y*y)
             th=ATAN2(x,y)
             IF(th.LT.0.D0) th=th+TWOPI
             IF(r >= d0rw) THEN
                gxy(nx,ny)=0.D0
             ELSE
                CALL SPL2DF(r,th,gxy(nx,ny),rhonrho,chig,ugz, &
                            nrhomax,nrhomax,nchig,ierr)
             END IF
!             WRITE(6,'(2I5,1P5E12.4)') nx,ny,x,y,r,th,gxy(nx,ny)
          END DO
       END DO
    END IF
             
    WRITE(LINE,'(A,I3,A)') '@diguv(',inum,')@'

    CALL PAGES
    SELECT CASE(ID)
    CASE(1)
       CALL GRD2D(0,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=4,TITLE_SIZE=0.4D0)
    CASE(2)
       CALL GRD2D(0,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=1,TITLE_SIZE=0.4D0)
    CASE(3)
       CALL GRD2D(0,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=2,TITLE_SIZE=0.4D0)
    CASE(11)
       CALL GRD2D(0,gx,gy,gxy,nxmax,nxmax,nymax, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=4,TITLE_SIZE=0.4D0)
    CASE(12)
       CALL GRD2D(0,gx,gy,gxy,nxmax,nxmax,nymax, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=1,TITLE_SIZE=0.4D0,&
                  XMIN=-d0rw,XMAX=d0rw,YMIN=-d0rw,YMAX=d0rw)
    CASE(13)
       CALL GRD2D(0,gx,gy,gxy,nxmax,nxmax,nymax, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=2,TITLE_SIZE=0.4D0)
    CASE(14)
       CALL GRD2D(0,gx,gy,gxy,nxmax,nxmax,nymax, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=11,TITLE_SIZE=0.4D0)
    CASE(15)
       CALL GRD2D(0,gx,gy,gxy,nxmax,nxmax,nymax, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=12,TITLE_SIZE=0.4D0)
    END SELECT
    CALL PAGEE

    IF(ID.GT.10) THEN
       DEALLOCATE(ddx,ddy,ddxy)
       DEALLOCATE(ugz)
       DEALLOCATE(gx,gy,gxy)
    END IF
    RETURN
  END SUBROUTINE T2_GC

  SUBROUTINE T2_GR(INUM)
    USE libgrf,ONLY: GRD1D
    USE T2COMM, ONLY: &
         i0ikind,i0rkind,twopi,i0xmax,d1guv,i0vmax, &
         i0lmax,i0pdiv_number,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax
    IMPLICIT NONE
    INTEGER(i0ikind),INTENT(IN):: inum
    REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: gz
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE:: gzl,dgzl
    REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: ugzl
    INTEGER(i0ikind):: nchi,nl,nrho,nchimaxl,nr,ierr
    CHARACTER(LEN=80):: LINE

    ALLOCATE(gz(nrhomax,nchimax+1))
    ALLOCATE(gzl(nchimax+1),dgzl(nchimax+1),ugzl(4,nchimax+1))
    DO nchi=1,nchimax+1
       gz(1,nchi)=d1guv(inum)
    END DO
    DO nrho=2,nrhomax
       nl=nlnrho(nrho)
       nchimaxl=i0pdiv_number*2**(i1mlvl(nl)-1)
       IF(nchimaxl==nchimax) THEN
          DO nchi=1,nchimax
!             write(6,'(A,6I5)') &
!                  'nrho,nnnrho,nchi,i0vmax,i0,i=', &
!                  nrho,nnnrho(nrho),nchi,i0vmax, &
!                  i0vmax*(nnnrho(nrho)+nchi-2), &
!                  i0vmax*(nnnrho(nrho)+nchi-2)+inum
             gz(nrho,nchi)=d1guv(i0vmax*(nnnrho(nrho)+nchi-2)+inum)
          END DO
       ELSE
          DO nchi=1,nchimaxl
             gzl(nchi)=d1guv(i0vmax*(nnnrho(nrho)+nchi-2)+inum)
          END DO
          gzl(nchimaxl+1)=d1guv(i0vmax*(nnnrho(nrho)-1)+inum)
          CALL SPL1D(chinl(1:nchimaxl+1,nl),gzl,dgzl,ugzl,nchimaxl+1,4,ierr)
          DO nchi=1,nchimax
             CALL SPL1DF(chinchi(nchi),gz(nrho,nchi), &
                         chinl(1:nchimaxl+1,nl),ugzl,nchimax+1,ierr)
          END DO
       END IF
    END DO

    WRITE(LINE,'(A,I3,A)') '@diguv(',inum,')@'
    CALL PAGES
    CALL GRD1D(0,rhonrho,gz,nrhomax,nrhomax,nchimax, &
               TITLE=LINE,MODE_LS=0)
    CALL PAGEE

    RETURN
  END SUBROUTINE T2_GR

END Module T2GOUT
