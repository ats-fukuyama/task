!  ***** TASK/XX GRAPHIC OUTPU *****

Module picgout
  PRIVATE
  PUBLIC pic_gout

CONTAINS

  SUBROUTINE pic_gout

    USE piccomm
    USE picparm,ONLY: pic_parm
    USE grd2dp_mod,ONLY: grd2dp
    USE libgrf
    IMPLICIT NONE
    INTEGER:: kid,mode,err,ich1,ich2,nx,ny,np,npo,nto
    CHARACTER(LEN=2):: kch
    CHARACTER(LEN=80):: line
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: work
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:,:):: workp
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: x,y,vtote,vtoti
    REAL(rkind):: WORK_RGB(3,1)
    INTEGER,PARAMETER:: nlmax=10
    INTEGER:: WORK_PAT(nlmax)
    REAL(rkind):: aspect,vtotemax,vtotimax
    aspect=DBLE(nymax)/DBLE(nxmax)

    ALLOCATE(x(nxmax1),y(nymax1))
    ALLOCATE(vtote(npmax),vtoti(npmax))
    DO nx=1,nxmax1
       x(nx)=DBLE(nx-1)
    END DO
    DO ny=1,nymax1
       y(ny)=DBLE(ny-1)
    END DO
    vtotemax=0.D0
    vtotimax=0.D0
    DO np=1,npmax
       vtote(np)=SQRT(vparae(np)**2+vperpe(np)**2)
       vtoti(np)=SQRT(vparai(np)**2+vperpi(np)**2)
       vtotemax=MAX(vtotemax,vtote(np))
       vtotimax=MAX(vtotimax,vtoti(np))
    END DO

1   CONTINUE
    err=0
    WRITE(6,'(A)') &
         '#### PIC GOUT: T1 E1-8 F1 O1-2 P1-4 X/exit'
    CALL TASK_KLIN(line,kid,mode,pic_parm)
    IF(mode == 2 .OR. mode == 3) GOTO 1

    ICH1=ICHAR(LINE(1:1))
    IF(ICH1.GE.97.AND.ICH1.LE.122) ICH1=ICH1-32
    ICH2=ICHAR(LINE(2:2))
    IF(ICH2.GE.97.AND.ICH2.LE.122) ICH2=ICH2-32
    KCH=CHAR(ICH1)//CHAR(ICH2)

    SELECT CASE(kch)
    CASE('T1')
       ALLOCATE(work(ntgmax,7))
       CALL PAGES
       !work(1:ntgmax,1)=aktott(1:ntgmax)
       !work(1:ntgmax,2)=akinit(1:ntgmax)
       !work(1:ntgmax,3)=akinet(1:ntgmax)
       !CALL GRD1D(1,timet,work,ntgmax,ntgmax,3,'@Ktot,Ki,Ke vs t@')
       !work(1:ntgmax,1)=aptott(1:ntgmax)
       !work(1:ntgmax,2)=apotet(1:ntgmax)
       !work(1:ntgmax,3)=apotmt(1:ntgmax)
       !CALL GRD1D(2,timet,work,ntgmax,ntgmax,3,'@Ptot,Pe,Pm vs t@')
       !work(1:ntgmax,1)=atott(1:ntgmax)
       !work(1:ntgmax,2)=aktott(1:ntgmax)
       !work(1:ntgmax,3)=aptott(1:ntgmax)
       !CALL GRD1D(3,timet,work,ntgmax,ntgmax,3,'@Wtot,Ktot,Ptot vs t@')
       work(1:ntgmax,1)=atott(1:ntgmax)
       work(1:ntgmax,2)=aktott(1:ntgmax)
       work(1:ntgmax,3)=akinet(1:ntgmax)
       work(1:ntgmax,4)=akinit(1:ntgmax)
       work(1:ntgmax,5)=aptott(1:ntgmax)
       work(1:ntgmax,6)=apotet(1:ntgmax)
       work(1:ntgmax,7)=apotmt(1:ntgmax)
       CALL GRD1D(0,timet,work,ntgmax,ntgmax,7, &
                  '@Wtot,Ktot,Ke,Ki,Ptot,Pe,Pm vs t@')
       CALL PAGEE
       DEALLOCATE(work)
    CASE('E1')
       CALL PAGES
       CALL GRD2D(1,x,y,ex,nxmax1,nxmax1,nymax1,'@Ex@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(2,x,y,ey,nxmax1,nxmax1,nymax1,'@Ey@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(3,x,y,ez,nxmax1,nxmax1,nymax1,'@Ez@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(4,x,y,rho,nxmax1,nxmax1,nymax1,'@rho@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
!       CALL GRD2D(4,x,y,phi,nxmax1,nxmax1,nymax1,'@phi@',ASPECT=aspect, &
!                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
!                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E2')
       CALL PAGES
       CALL GRD2D(1,x,y,bx,nxmax1,nxmax1,nymax1,'@Bx@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(2,x,y,by,nxmax1,nxmax1,nymax1,'@By@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(3,x,y,bz,nxmax1,nxmax1,nymax1,'@Bz@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E3')
       CALL PAGES
       CALL GRD2D(1,x,y,Ax,nxmax1,nxmax1,nymax1,'@Ax@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(2,x,y,Ay,nxmax1,nxmax1,nymax1,'@Ay@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(3,x,y,Az,nxmax1,nxmax1,nymax1,'@Az@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(4,x,y,AA,nxmax1,nxmax1,nymax1,'@AA@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E4')
       CALL PAGES
       CALL GRD2D(1,x,y,jx,nxmax1,nxmax1,nymax1,'@jx@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(2,x,y,jy,nxmax1,nxmax1,nymax1,'@jy@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(3,x,y,jz,nxmax1,nxmax1,nymax1,'@jz@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(4,x,y,rho,nxmax1,nxmax1,nymax1,'@rho@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E5')
       CALL PAGES
         CALL GRD2D( 1,x,y,ex,nxmax1,nxmax1,nymax1,'@Ex@',ASPECT=aspect, &
                      XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax),&
                      NLMAX=nlmax,LINE_PAT=WORK_PAT)
         CALL GRD2D( 2,x,y,ey,nxmax1,nxmax1,nymax1,'@Ey@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                       NLMAX=nlmax,LINE_PAT=WORK_PAT)
         CALL GRD2D( 3,x,y,ez,nxmax1,nxmax1,nymax1,'@Ez@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                      NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D( 8,x,y,esx,nxmax1,nxmax1,nymax1,'@ESx@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !               NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D( 9,x,y,esy,nxmax1,nxmax1,nymax1,'@ESy@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !               NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D(10,x,y,esz,nxmax1,nxmax1,nymax1,'@ESz@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !               NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D(11,x,y,emx,nxmax1,nxmax1,nymax1,'@EMx@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !               NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D(12,x,y,emy,nxmax1,nxmax1,nymax1,'@EMy@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !              NLMAX=nlmax,LINE_PAT=WORK_PAT)
       !  CALL GRD2D(13,x,y,emz,nxmax1,nxmax1,nymax1,'@EMz@',ASPECT=aspect, &
       !               XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
       !               NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E6')
       CALL PAGES
       CALL GRD1D( 1,x,ex,nxmax1,nxmax1,nymax1,'@Ex(x)@', &
                      XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 2,x,ey,nxmax1,nxmax1,nymax1,'@Ey(x)@', &
                      XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 3,x,ez,nxmax1,nxmax1,nymax1,'@Ez(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D( 8,x,esx,nxmax1,nxmax1,nymax1,'@ESx(x)@', &
       !               XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D( 9,x,esy,nxmax1,nxmax1,nymax1,'@ESy(x)@', &
       !              XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D(10,x,esz,nxmax1,nxmax1,nymax1,'@ESz(x)@', &
       !              XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D(11,x,emx,nxmax1,nxmax1,nymax1,'@EMx(x)@', &
       !              XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D(12,x,emy,nxmax1,nxmax1,nymax1,'@EMy(x)@', &
       !              XMIN=0.D0,XMAX=DBLE(nxmax))
       ! CALL GRD1D(13,x,emz,nxmax1,nxmax1,nymax1,'@EMz(x)@', &
       !              XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL PAGEE
    CASE('E7')
       CALL PAGES
       CALL GRD2D( 5,x,y,rho,nxmax1,nxmax1,nymax1,'@rho@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D( 6,x,y,phi,nxmax1,nxmax1,nymax1,'@phi@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D( 7,x,y,bb,nxmax1,nxmax1,nymax1,'@BB@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D( 8,x,y,jx,nxmax1,nxmax1,nymax1,'@jx@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D( 9,x,y,jy,nxmax1,nxmax1,nymax1,'@jy@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(10,x,y,jz,nxmax1,nxmax1,nymax1,'@jz@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(11,x,y,Bx,nxmax1,nxmax1,nymax1,'@Bx@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(12,x,y,By,nxmax1,nxmax1,nymax1,'@By@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL GRD2D(13,x,y,Bz,nxmax1,nxmax1,nymax1,'@Bz@',ASPECT=aspect, &
                    XMIN=0.D0,XMAX=DBLE(nxmax),YMIN=0.D0,YMAX=DBLE(nymax), &
                    NLMAX=nlmax,LINE_PAT=WORK_PAT)
       CALL PAGEE
    CASE('E8')
       CALL PAGES
       CALL GRD1D( 5,x,rho,nxmax1,nxmax1,nymax1,'@rho(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 6,x,phi,nxmax1,nxmax1,nymax1,'@phi(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 7,x,bb,nxmax1,nxmax1,nymax1,'@BB(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 8,x,jx,nxmax1,nxmax1,nymax1,'@jx(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 9,x,jy,nxmax1,nxmax1,nymax1,'@jy(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(10,x,jz,nxmax1,nxmax1,nymax1,'@jz(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(11,x,Bx,nxmax1,nymax1,nymax1,'@Bx(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(12,x,By,nxmax1,nxmax1,nymax1,'@By(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(13,x,Bz,nxmax1,nxmax1,nymax1,'@Bz(x)@', &
                    XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL PAGEE
    CASE('F1')
       CALL PAGES
       WORK_RGB(1,1)=0.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=1.D0
       CALL GRD2DP(1,vparae,vperpe,vtote,npmax,'@vpara-vperp:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,XMIN=-vtotemax,XMAX=vtotemax, &
                   YMIN=0.D0,YMAX=vtotemax,ASPECT=0.5D0)
       CALL GRD2DP(3,xe,ye,ze,npmax,'@x-y-z:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,ASPECT=aspect, &
                   XMIN=0.D0,XMAX=DBLE(nxmax), &
                   YMIN=0.D0,YMAX=DBLE(nymax))
       WORK_RGB(1,1)=1.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=0.D0
       CALL GRD2DP(2,vparai,vperpi,vtoti,npmax,'@vpara-vperp:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,XMIN=-vtotimax,XMAX=vtotimax, &
                   YMIN=0.D0,YMAX=vtotimax,ASPECT=0.5D0)
       CALL GRD2DP(4,xi,yi,zi,npmax,'@x-y-z:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,ASPECT=aspect, &
                   XMIN=0.D0,XMAX=DBLE(nxmax), &
                   YMIN=0.D0,YMAX=DBLE(nymax))
       CALL PAGEE
    CASE('O1')
       ALLOCATE(work(ntomax,3))
       DO npo=1,npomax
          work(1:ntomax,1)=xpo(npo,1:ntomax)
          work(1:ntomax,2)=ypo(npo,1:ntomax)
          work(1:ntomax,3)=zpo(npo,1:ntomax)
          CALL PAGES
          WORK_RGB(1,1)=1.D0
          WORK_RGB(2,1)=0.D0
          WORK_RGB(3,1)=0.D0
          CALL GRD2DP(1,work(1:ntomax,1),work(1:ntomax,2),work(1:ntomax,3), &
                      ntomax,'@xp-yp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=DBLE(nymax)/DBLE(nxmax), &
                      XMIN=0.D0,XMAX=DBLE(nxmax), &
                      YMIN=0.D0,YMAX=DBLE(nymax))
          CALL GRD2DP(2,work(1:ntomax,2),work(1:ntomax,3),work(1:ntomax,1), &
                      ntomax,'@yp-zp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=DBLE(nzmax)/DBLE(nymax), &
                      XMIN=0.D0,XMAX=DBLE(nymax), &
                      YMIN=0.D0,YMAX=DBLE(nzmax))
          CALL GRD2DP(3,work(1:ntomax,3),work(1:ntomax,1),work(1:ntomax,2), &
                      ntomax,'@zp-xp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=DBLE(nxmax)/DBLE(nzmax), &
                      XMIN=0.D0,XMAX=DBLE(nzmax), &
                      YMIN=0.D0,YMAX=DBLE(nxmax))
          CALL PAGEE
       END DO
       DEALLOCATE(work)
    CASE('O2')
       ALLOCATE(work(ntomax,3))
       DO npo=1,npomax
          work(1:ntomax,1)=xpo(npo,1:ntomax)
          work(1:ntomax,2)=ypo(npo,1:ntomax)
          work(1:ntomax,3)=zpo(npo,1:ntomax)
          CALL PAGES
          WORK_RGB(1,1)=1.D0
          WORK_RGB(2,1)=0.D0
          WORK_RGB(3,1)=0.D0
          CALL GRD2DP(1,work(1:ntomax,1),work(1:ntomax,2),work(1:ntomax,3), &
                      ntomax,'@xp-yp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=0.D0)
          CALL GRD2DP(2,work(1:ntomax,2),work(1:ntomax,3),work(1:ntomax,1), &
                      ntomax,'@yp-zp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=0.D0)
          CALL GRD2DP(3,work(1:ntomax,3),work(1:ntomax,1),work(1:ntomax,2), &
                      ntomax,'@zp-xp@',NMMAX=1, &
                      MARK_RGB=WORK_RGB,ASPECT=0.D0)
          CALL PAGEE
       END DO
       DEALLOCATE(work)
    CASE('P1')
       ALLOCATE(workp(0:nxmax,ntpmax,9))
       CALL PAGES
       CALL sum_over_y(nxmax,nymax,ntpmax,profilee,workp)
        CALL GRD1D( 5,x,workp(0:nxmax,1:ntpmax,1),nxmax1,nxmax1,ntpmax, &
                        '@ne(x)@',FMIN=0.D0, &
                        XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 8,x,workp(0:nxmax,1:ntpmax,2),nxmax1,nxmax1,ntpmax, &
                        '@vxe(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 9,x,workp(0:nxmax,1:ntpmax,3),nxmax1,nxmax1,ntpmax, &
                        '@vye(x)@', &
                        XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D(10,x,workp(0:nxmax,1:ntpmax,4),nxmax1,nxmax1,ntpmax, &
                        '@vze(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 6,x,workp(0:nxmax,1:ntpmax,5),nxmax1,nxmax1,ntpmax, &
                        '@vparae(x)@', &
                        XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D( 7,x,workp(0:nxmax,1:ntpmax,6),nxmax1,nxmax1,ntpmax, &
                        '@vperpe(x)@', &
                        XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D(11,x,workp(0:nxmax,1:ntpmax,7),nxmax1,nxmax1,ntpmax, &
                       '@Tparae(x)@',FMIN=0.D0, &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D(12,x,workp(0:nxmax,1:ntpmax,8),nxmax1,nxmax1,ntpmax, &
                         '@Tperpe(x)@',FMIN=0.D0, &
                         XMIN=0.D0,XMAX=DBLE(nxmax))
        CALL GRD1D(13,x,workp(0:nxmax,1:ntpmax,9),nxmax1,nxmax1,ntpmax, &
                          '@Te(x)@',FMIN=0.5D0, FMAX=2.5d0,&
                          XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL PAGEE
       DEALLOCATE(workp)
    CASE('P2')
       CALL PAGES
       CALL GRD2D( 5,x,y,profilee(0:nxmax,0:nymax,1,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@ne(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 8,x,y,profilee(0:nxmax,0:nymax,2,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vxe(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 9,x,y,profilee(0:nxmax,0:nymax,3,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vye(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(10,x,y,profilee(0:nxmax,0:nymax,4,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vze(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 6,x,y,profilee(0:nxmax,0:nymax,5,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vparae(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 7,x,y,profilee(0:nxmax,0:nymax,6,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vperpe(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(11,x,y,profilee(0:nxmax,0:nymax,7,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Tparae(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(12,x,y,profilee(0:nxmax,0:nymax,8,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Tperpe(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(13,x,y,profilee(0:nxmax,0:nymax,9,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Te(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL PAGEE
    CASE('P3')
       ALLOCATE(workp(0:nxmax,ntpmax,9))
       CALL PAGES
       CALL sum_over_y(nxmax,nymax,ntpmax,profilei,workp)
       CALL GRD1D( 5,x,workp(0:nxmax,1:ntpmax,1),nxmax1,nxmax1,ntpmax, &
                       '@ni(x)@',FMIN=0.D0, &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 8,x,workp(0:nxmax,1:ntpmax,2),nxmax1,nxmax1,ntpmax, &
                       '@vxi(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 9,x,workp(0:nxmax,1:ntpmax,3),nxmax1,nxmax1,ntpmax, &
                       '@vyi(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(10,x,workp(0:nxmax,1:ntpmax,4),nxmax1,nxmax1,ntpmax, &
                       '@vzi(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 6,x,workp(0:nxmax,1:ntpmax,5),nxmax1,nxmax1,ntpmax, &
                       '@vparai(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D( 7,x,workp(0:nxmax,1:ntpmax,6),nxmax1,nxmax1,ntpmax, &
                       '@vperpi(x)@', &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(11,x,workp(0:nxmax,1:ntpmax,7),nxmax1,nxmax1,ntpmax, &
                       '@Tparai(x)@',FMIN=0.D0, &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(12,x,workp(0:nxmax,1:ntpmax,8),nxmax1,nxmax1,ntpmax, &
                       '@Tperpi(x)@',FMIN=0.D0, &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL GRD1D(13,x,workp(0:nxmax,1:ntpmax,9),nxmax1,nxmax1,ntpmax, &
                       '@Ti(x)@',FMIN=0.D0, &
                       XMIN=0.D0,XMAX=DBLE(nxmax))
       CALL PAGEE
       DEALLOCATE(workp)
    CASE('P4')
       CALL PAGES
       CALL GRD2D( 5,x,y,profilei(0:nxmax,0:nymax,1,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@ni(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 8,x,y,profilei(0:nxmax,0:nymax,2,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vxi(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 9,x,y,profilei(0:nxmax,0:nymax,3,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vyi(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(10,x,y,profilei(0:nxmax,0:nymax,4,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vzi(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 6,x,y,profilei(0:nxmax,0:nymax,5,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vparai(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D( 7,x,y,profilei(0:nxmax,0:nymax,6,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@vperpi(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(11,x,y,profilei(0:nxmax,0:nymax,7,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Tparai(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(12,x,y,profilei(0:nxmax,0:nymax,8,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Tperpi(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL GRD2D(13,x,y,profilei(0:nxmax,0:nymax,9,ntpmax), &
                       nxmax1,nxmax1,nymax1,'@Ti(x,y)@',ASPECT=aspect, &
                       XMIN=0.D0,XMAX=DBLE(nxmax), &
                       YMIN=0.D0,YMAX=DBLE(nymax))
       CALL PAGEE
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    DEALLOCATE(x,y)
    RETURN
  END SUBROUTINE pic_gout

!-----

  SUBROUTINE sum_over_y(nxmax,nymax,ntpmax,fxy,fx)

    IMPLICIT NONE
    INTEGER:: nxmax,nymax,ntpmax
    REAL(8):: fxy(0:nxmax,0:nymax,9,ntpmax),fx(0:nxmax,ntpmax,9)
    INTEGER:: nx,ny,ntp,i

    DO i=1,9
       DO ntp=1,ntpmax
          DO nx=0,nxmax
             fx(nx,ntp,i)=0.5D0*(fxy(nx,0,i,ntp) &
                                +fxy(nx,nymax,i,ntp))
             DO ny=1,nymax-1
                fx(nx,ntp,i)=fx(nx,ntp,i)+fxy(nx,ny,i,ntp)
             END DO
             fx(nx,ntp,i)=fx(nx,ntp,i)/DBLE(nymax)
          END DO
       END DO
    END DO
  END SUBROUTINE sum_over_y

END Module picgout
