! wfhout.f90

MODULE wfhout

  INTEGER,PARAMETER:: dp = KIND(1.0D0)
  INTEGER:: init=0
  INTEGER:: idata=0
  INTEGER:: nconf,n1max,n2max
  INTEGER:: nzm_min,nzm_max,nzm_step
  INTEGER:: iopt_aspect,id2
  REAL(dp):: xdmin,xdmax,ydmin,ydmax
  REAL(dp):: r_aspect
  REAL(dp),ALLOCATABLE:: xg(:),yg(:),fg(:,:)
  REAL(dp),ALLOCATABLE:: gpos(:,:)
  COMPLEX(dp),DIMENSION(:,:,:),ALLOCATABLE:: CAFA,CAFXA,CAFYA,CAFZA
  COMPLEX(dp),DIMENSION(:,:,:),ALLOCATABLE:: CARA,CARXA,CARYA,CARZA
  INTEGER,ALLOCATABLE:: nelmxy(:,:)
  
  PRIVATE
  PUBLIC wf_hout

CONTAINS

! New graphic module

  SUBROUTINE wf_hout(xd,yd,nnod,nsmax,nzmax)

    USE libgrf
    USE libchar
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nnod,nsmax,nzmax
    REAL(dp),INTENT(IN):: xd(nnod),yd(nnod)
    INTEGER,PARAMETER:: NCHMAX=80
    INTEGER,PARAMETER:: NWDM=12
    CHARACTER(LEN=NCHMAX):: KLINE,KWORD(NWDM),KWD,KWTEMP
    CHARACTER(LEN=1):: KID,KG1,KG2,KG3,KG4
    INTEGER:: NL,NWD,NCH,nfig,nwdmax,nconf_temp
    INTEGER:: nzm_temp,nfigxmax,nfigymax,nfigmax,nfigx,nfigy
    INTEGER:: ng_min,ng_max,ng_step,ng,nx,ny
    REAL(dp):: rvpage,fsizex,fsizey,chsize,chsize_real

! --- initializaion ---

    IF(init.EQ.0) THEN
       
! --- set boundary of x and y ---

       xdmin=MINVAL(xd)+1.D-6
       xdmax=MAXVAL(xd)-1.D-6
       ydmin=MINVAL(yd)+1.D-6
       ydmax=MAXVAL(yd)-1.D-6

! --- set aspect ratio ---

       iopt_aspect=0
       r_aspect=(ydmax-ydmin)/(xdmax-xdmin)

! --- set mode number range ---

       nzm_min=1
       IF(NZMAX.GT.1) THEN
          nzm_max=NZMAX
       ELSE
          nzm_max=1
       END IF
       nzm_step=1
       n1max=100
       n2max=100

! --- setup element table ---

       CALL setup_region_xy
       nconf=1
       
       init=1
    END IF
    
! --- setup variable array ---
    
    ALLOCATE(CAFA(4,n1max,n2max),CAFXA(4,n1max,n2max))
    ALLOCATE(CAFYA(4,n1max,n2max),CAFZA(4,n1max,n2max))
    ALLOCATE(CARA(4,n1max,n2max),CARXA(4,n1max,n2max))
    ALLOCATE(CARYA(4,n1max,n2max),CARZA(4,n1max,n2max))
    
! --- adjust for the change of NZMAX ---

    IF(NZMAX.EQ.1) THEN
       nzm_min=1
       nzm_max=1
       nzm_step=1
    END IF
    IF(nzm_max.GT.NZMAX) nzm_max=NZMAX
    IF(nzm_step.GT.nzm_max-nzm_min) nzm_step=MAX(nzm_max-nzm_min,1)
    
! --- input graphic commands ---

1   WRITE(6,'(A)') '# INPUT : A,X/Y/Z/F,F/R,R/I'
    WRITE(6,'(A)') '          R=region N=mode On=option X=exit'
    READ(5,'(A80)',ERR=1,END=9000) KLINE

! --- separate KLINE int KWORD ---

9   NL=0
    NWD=0
    NCH=0
10  IF(NL.GE.80) GOTO 20
    NL=NL+1
    KID=KLINE(NL:NL)
    CALL toupper(KID)
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

20  NWDMAX=NWD

! --- spectial commands ('X', 'R', 'N') at the top of line ---

    KWD=KWORD(1)
    KG1=KWD(1:1)
    KG2=KWD(2:2)

! --- KG1='X' for exit  ---

    IF(KG1.EQ.'X') GO TO 9000

! --- KG1='R' for setting plot area  ---

    SELECT CASE(KG1)
    CASE('R')
101    WRITE(6,'(A)') '# INPUT: nconf,xdmin,xdmax,ydmin,ydmax,n1max,n2max'
       WRITE(6,'(A)') '           nconf=1,4:XY, 2,5:XZ, 3,6:YZ, 0:end'
       WRITE(6,'(A,I3,1P4E12.4,2I5)') &
            '       ',nconf,xdmin,xdmax,ydmin,ydmax,n1max,n2max
       nconf_temp=nconf
       READ(5,*,ERR=101,END=1) nconf_temp,xdmin,xdmax,ydmin,ydmax,n1max,n2max
       IF(nconf_temp.EQ.0) GOTO 1 
       nconf=nconf_temp

       r_aspect=(ydmax-ydmin)/(xdmax-xdmin)
       IF(iopt_aspect.EQ.1) THEN
          IF(r_aspect.GT.1.D0) THEN
             r_aspect=3.D0/2.D0
          ELSE
             r_aspect=2.D0/3.D0
          END IF
       END IF
       GOTO 1

! --- KG1='N' for setup mode number  ---

    CASE('N')
102    WRITE(6,'(A,3I5)') '# INPUT: nzm_min,nzm_max,nzm_step:', &
                                    nzm_min,nzm_max,nzm_step
       nzm_temp=nzm_min
       READ(5,*,ERR=102,END=1) nzm_temp,nzm_max,nzm_step
       IF(nzm_temp.EQ.0) GOTO 1 
       nzm_min=nzm_temp
       GOTO 1

! --- KG1='O' for setup options  ---

    CASE('O')
       SELECT CASE(KG2)
       CASE('0')
          iopt_aspect=0
          WRITE(6,'(A,1PE12.4)') '== Real aspect ratio:', &
               (ydmax-ydmin)/(xdmax-xdmin)
       CASE('1')
          iopt_aspect=1
          WRITE(6,'(A)') '== Optimum aspect ratio =='
       END SELECT
       GO TO 1
    END SELECT

! --- setup nelmxy table for quick local data access ---

    SELECT CASE(nconf)
    CASE(1)
       CALL setup_region_xy
    CASE(2)
!       CALL setup_region_xz(xdmin,xdmax,n1max,ydmin,ydmax,n2max, &
!                            xg,yg,nelmxy)
    CASE(3)
!       CALL setup_region_yz(xdmin,xdmax,n1max,ydmin,ydmax,n2max, &
!                            xg,yg,nelmxy)
    END SELECT

! --- count number of figures ---

    nfig=0
    NWD=0
21  NWD=NWD+1
    IF(NWD.GT.NWDMAX) GOTO 30
    KWD=KWORD(NWD)
    KG1=KWD(1:1)
    SELECT CASE(KG1)
    CASE('A','E','B')
       nfig=nfig+(nzm_max-nzm_min+1)/nzm_step
    CASE('P')
       nfig=nfig+NSMAX
    CASE('M')
       nfig=nfig+1
    CASE('N','T')
       nfig=nfig+NSMAX
    END SELECT
    GO TO 21

! --- set figure size ---

30  IF(nfig.EQ.0) GOTO 1
    nfigmax=nfig
    rvpage=nfigmax/r_aspect
    nfigymax=INT(rvpage-0.001D0)+1
    nfigxmax=(nfigmax-1)/nfigymax+1
    write(6,'(A,3I5)') 'nfigmax,nfigxmax,nfigymax=',nfigmax,nfigxmax,nfigymax

! --- set figure position ---

    fsizex=25.6D0/nfigxmax
    fsizey=17.1D0/nfigymax
    chsize_real=0.2D0
    chsize=chsize_real*r_aspect
    IF(ALLOCATED(gpos)) DEALLOCATE(gpos)
    ALLOCATE(GPOS(1:4,nfigmax))
    DO nfig=1,nfigmax
       nfigy=(nfig-1)/nfigxmax+1
       nfigx=nfig-nfigxmax*(nfigy-1)
       gpos(1,nfig)=fsizex*(nfigx-1)          +6*chsize_real
       gpos(2,nfig)=fsizex* nfigx
       gpos(3,nfig)=fsizey*(nfigymax-nfigy)   +3*chsize_real
       gpos(4,nfig)=fsizey*(nfigymax-nfigy+1) -3*chsize_real
    END DO

! --- set first KWORD and extract KG1..4  ---

    CALL PAGES

    nwd=0
    nfig=0
40  NWD=NWD+1
    IF(NWD.GT.NWDMAX) GOTO 1000
    KWD=KWORD(NWD)
    KG1=KWD(1:1)
    KG2=KWD(2:2)
    KG3=KWD(3:3)
    KG4=KWD(4:4)

    SELECT CASE(KG1)
    CASE('A','E','B')
       ng_min=nzm_min
       ng_max=nzm_max
       ng_step=nzm_step
    CASE('P','N','T')
       SELECT CASE(KG2)
       CASE('1','2','3','4','5','6','7','8','9')
          READ(KG1,'(I1)') ng
          ng_min=ng
          ng_max=ng
          ng_step=1
       CASE DEFAULT
          ng_min=1
          ng_max=NSMAX
          ng_step=1
       END SELECT
    END SELECT
!    write(6,'(A,3I5)') 'ng_min,ng_max,ng_step=',ng_min,ng_max,ng_step
    
! --- setup field data CA and its derivatives for nx and ny' ---

    DO ng=ng_min,ng_max,ng_step

       nfig=nfig+1
       
!       write(6,'(A,2I5,A)') 'ng,nfig,KG1=',ng,nfig,KG1
       SELECT CASE(KG1)
       CASE('A','E','B')
          IF(IDATA.EQ.1.OR.IDATA.EQ.3) THEN
             DO ny=1,n2max
                DO nx=1,n1max
                   CALL FIELDA(xg(nx),yg(ny),nelmxy(nx,ny),ng,0, &
                               CAFA(1:4,nx,ny),CAFXA(1:4,nx,ny), &
                               CAFYA(1:4,nx,ny),CAFZA(1:4,nx,ny))
                   CALL FIELDA(xg(nx),yg(ny),nelmxy(nx,ny),ng,1, &
                               CARA(1:4,nx,ny),CARXA(1:4,nx,ny), &
                               CARYA(1:4,nx,ny),CARZA(1:4,nx,ny))
                END DO
             END DO
             write(6,'(1P4E12.4)') CAFA(1,2,2),CAFA(2,2,2)
             write(6,'(1P4E12.4)') CAFA(3,2,2),CAFA(4,2,2)
             IDATA=2
          END IF
       CASE('P','N','T')
          IF(IDATA.EQ.1.OR.IDATA.EQ.2) THEN
             IDATA=3
          END IF
       END SELECT

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
       END SELECT

       CALL GRD2D(0,xg,yg,fg,n1max,n1max,n2max, &
                  GPXMIN=gpos(1,nfig),GPXMAX=gpos(2,nfig), &
                  GPYMIN=gpos(3,nfig),GPYMAX=gpos(4,nfig), &
                  TITLE_SIZE=chsize, &
                  ASPECT=r_aspect)
    END DO
    GO TO 40
    
1000 CALL PAGEE
    GO TO 1

9000 CONTINUE
    DEALLOCATE(CAFA,CAFXA,CAFYA,CAFZA)
    DEALLOCATE(CARA,CARXA,CARYA,CARZA)
    DEALLOCATE(gpos,xg,yg,nelmxy)
    RETURN
  END SUBROUTINE wf_hout

! --- setup mesh data and element data of mesh points ---

    SUBROUTINE setup_region_xy
    IMPLICIT NONE
    REAL(dp):: delx,dely
    INTEGER:: nx,ny
    REAL(dp),SAVE:: xdmin_save=0.D0,xdmax_save=0.D0
    REAL(dp),SAVE:: ydmin_save=0.D0,ydmax_save=0.D0
    INTEGER,SAVE:: n1max_save=0,n2max_save=0

!    write(6,'(1P4E12.4)') xdmin,xdmin_save,xdmax,xdmax_save
!    write(6,'(1P4E12.4)') ydmin,ydmin_save,ydmax,ydmax_save
!    write(6,'(4I5)') n1max,n1max_save,n2max,n2max_save

    IF(xdmin.EQ.xdmin_save.AND. &
       xdmax.EQ.xdmax_save.AND. &
       ydmin.EQ.ydmin_save.AND. &
       ydmax.EQ.ydmax_save.AND. &
       n1max.EQ.n1max_save.AND. &
       n2max.EQ.n2max_save) RETURN
       
      IF(ALLOCATED(xg)) DEALLOCATE(xg)
      ALLOCATE(xg(n1max))
      delx=(xdmax-xdmin)/(n1max-1)
      DO nx=1,n1max
         xg(nx)=xdmin+delx*(nx-1)
      END DO

      IF(ALLOCATED(yg)) DEALLOCATE(yg)
      ALLOCATE(yg(n2max))
      dely=(ydmax-ydmin)/(n2max-1)
      DO ny=1,n2max
         yg(ny)=ydmin+delx*(ny-1)
      END DO

      IF(ALLOCATED(nelmxy)) DEALLOCATE(nelmxy)
      ALLOCATE(nelmxy(n1max,n2max))
      DO ny=1,n2max
         DO nx=1,n1max
            CALL WFFEP(xg(nx),yg(ny),nelmxy(nx,ny))
            WRITE(21,'(2I5,1P2E12.4,I5)') nx,ny,xg(nx),yg(ny),nelmxy(nx,ny)
         END DO
      END DO

      IF(ALLOCATED(fg)) DEALLOCATE(fg)
      ALLOCATE(fg(n1max,n2max))

      xdmin_save=xdmin
      xdmax_save=xdmax
      ydmin_save=ydmin
      ydmax_save=ydmax
      n1max_save=n1max
      n2max_save=n2max
      idata=1
    RETURN
  END SUBROUTINE setup_region_xy
END MODULE wfhout
