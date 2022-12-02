! wrgout.f90

MODULE WRGOUT

  USE wrcomm,ONLY: rkind

  PRIVATE
  PUBLIC wr_gout

  INTEGER,PARAMETER:: nlmax_p=10
!  REAL(rkind),PARAMETER:: line_rgb_nlmax(3,nlmax_p) &
!       =RESHAPE([1.0D0,0.0D0,0.0D0, &
!                 1.0D0,0.5D0,0.0D0, &
!                 1.0D0,0.8D0,0.0D0, &
!                 0.5D0,0.8D0,0.0D0, &
!                 0.2D0,0.8D0,0.0D0, &
!                 0.0D0,0.8D0,0.2D0, &
!                 0.0D0,0.8D0,0.5D0, &
!                 0.0D0,0.8D0,1.0D0, &
!                 0.0D0,0.5D0,1.0D0, &
!                 0.0D0,0.0D0,1.0D0],[3,10])
!  INTEGER,PARAMETER:: line_pat_nlmax(nlmax_p)=[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    
  
CONTAINS

!*********************** GRAPHIC DATA OUTPUT ************************

  SUBROUTINE WR_GOUT(NSTAT)

    USE wrcomm
    USE plgout
    USE libchar
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSTAT
    CHARACTER(LEN=1)::  KID
    INTEGER:: ns,nmax,nmax_save,n
    REAL(rkind),ALLOCATABLE:: rhona(:)

    ns=1
    nmax=0
    nmax_save=-1
    
    IF(MODELG.EQ.11) THEN
       CALL WRGRF11A
       CALL WRGRF11B
       RETURN
    END IF

1   WRITE(6,*) &
         '## INPUT GRAPH TYPE : ray:1,2,3,4,5,6,7 beam:A,B,C,D  P:prof  end:X'
    READ(5,'(A1)',ERR=1,END=9000) KID
    CALL toupper(KID)

    IF(KID.EQ.'1'.AND.NSTAT.GE.1) CALL WRGRF1
    IF(KID.EQ.'2'.AND.NSTAT.GE.1) CALL WRGRF2
    IF(KID.EQ.'3'.AND.NSTAT.GE.1) CALL WRGRF3
    IF(KID.EQ.'4'.AND.NSTAT.GE.1) CALL WRGRF4
    IF(KID.EQ.'5'.AND.NSTAT.GE.1) CALL WRGRF5
    IF(KID.EQ.'6'.AND.NSTAT.GE.1) CALL WRGRF6
    IF(KID.EQ.'7'.AND.NSTAT.GE.1) THEN
2      CONTINUE
       WRITE(6,'(A,2I6)') '## ns,nmax=',ns,nmax
       READ(5,*,ERR=2,END=1) ns,nmax
       IF(nmax.NE.nmax_save) THEN
          IF(ALLOCATED(rhona)) DEALLOCATE(rhona)
          IF(nmax.GE.1) ALLOCATE(rhona(nmax))
          nmax_save=nmax
          rhona(1:nmax)=0.D0
       END IF
3      CONTINUE
       WRITE(6,'(A)') '## rhona='
       WRITE(6,'(10F8.4)') (rhona(n),n=1,nmax)
       READ(5,*,ERR=3,END=2) (rhona(n),n=1,nmax)
       CALL WRGRF7(ns,nmax,rhona)
    END IF
    IF(KID.EQ.'A'.AND.NSTAT.GE.2) CALL WRGRFB1
    IF(KID.EQ.'B'.AND.NSTAT.GE.2) CALL WRGRFB2
    IF(KID.EQ.'C'.AND.NSTAT.GE.2) CALL WRGRFB3
    IF(KID.EQ.'D'.AND.NSTAT.GE.2) CALL WRGRFB4
    IF(KID.EQ.'P') CALL pl_gout
    IF(KID.EQ.'X') GOTO 9000

    GOTO 1

9000 RETURN
  END SUBROUTINE WR_GOUT


!     ***** POLOIDAL TRAJECTORY AND POWER *****

  SUBROUTINE WRGRF1

    USE wrcomm
    USE plprof,ONLY:PL_MAG_OLD,PL_MAG_TYPE,PL_MAG
    IMPLICIT NONE
    INTEGER,PARAMETER:: NSUM=501
    INTEGER,PARAMETER:: NGXL=101
    INTEGER,PARAMETER:: NGYL=101

    REAL,ALLOCATABLE:: GRS(:),GZS(:)
    REAL:: GX(NSTPMAX+1),GY(NSTPMAX+1)
    REAL:: GUX(NSTPMAX+1),GUY(NSTPMAX+1)
    REAL:: GPX(nrsmax),GPY(nrsmax,NRAYMAX)
    REAL:: GCF(NGXL,NGYL),GCX(NGXL),GCY(NGYL)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    TYPE(pl_mag_type):: MAG
    INTEGER:: NSU,NGX,NGY,NRAY,NSTP
    INTEGER:: ID,NRS
    REAL:: GRSMIN,GRSMAX,GRMIN,GRMAX,GRSTEP,GRORG,GRLEN
    REAL:: GZSMIN,GZSMAX,GZMIN,GZMAX,GZSTEP,GZLEN
    REAL:: GRMID,GZMID
    REAL:: GYSMIN,GYSMAX,GYSCAL
    REAL:: GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GXORG,GZSCAL
    REAL(rkind):: DGX,DGY,ZL,RL,RHON,XL,YL
    REAL(rkind):: RFL,wwl,wce
    EXTERNAL GMNMX1,GMNMX2,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT
    INTEGER:: NGULEN
    REAL:: GUCLIP

    ALLOCATE(KA(4,NGXL,NGYL))

    !  ----- PLASMA BOUNDARY -----

    ALLOCATE(GRS(nsumax+1),GZS(nsumax+1))
    DO NSU=1,NSUMAX+1
       GRS(NSU)=GUCLIP(RSU(NSU))
       GZS(NSU)=GUCLIP(ZSU(NSU))
    ENDDO

    CALL GMNMX1(GRS,1,NSUMAX,1,GRSMIN,GRSMAX)
    CALL GMNMX1(GZS,1,NSUMAX,1,GZSMIN,GZSMAX)
    CALL GQSCAL(GRSMIN,GRSMAX,GRMIN,GRMAX,GRSTEP)
    CALL GQSCAL(GZSMIN,GZSMAX,GZMIN,GZMAX,GZSTEP)

    GRORG=(INT(GRMIN/(2*GRSTEP))+1)*2*GRSTEP

    GRLEN=GRMAX-GRMIN
    GZLEN=GZMAX-GZMIN
    IF(1.5*GRLEN.GT.GZLEN) THEN
       GZMID=0.5*(GZMIN+GZMAX)
       GZMIN=GZMID-0.75*GRLEN
       GZMAX=GZMID+0.75*GRLEN
    ELSE
       GRMID=0.5*(GRMIN+GRMAX)
       GRMIN=GRMID-0.333*GZLEN
       GRMAX=GRMID+0.333*GZLEN
    ENDIF

    CALL PAGES
    CALL SETCHS(0.25,0.)
    CALL SETFNT(32)
    CALL SETLIN(0,2,7)

    CALL GDEFIN(1.5,11.5,1.0,16.0,GRMIN,GRMAX,GZMIN,GZMAX)
    CALL GFRAME
    CALL GSCALE(GRORG,  GRSTEP,0.0,  GZSTEP,0.1,9)
    CALL GVALUE(GRORG,2*GRSTEP,0.0,0.0,NGULEN(2*GRSTEP))
    CALL GVALUE(0.0,0.0,0.0,2*GZSTEP,NGULEN(2*GZSTEP))
    CALL SETLIN(0,2,4)
    CALL GPLOTP(GRS,GZS,1,NSUMAX+1,1,0,0,0)

!     ----- magnetic surface -----

    IF(MODELG.EQ.3.OR.MODELG.EQ.5.OR.MODELG.EQ.8) THEN
       DGX=DBLE(GRMAX-GRMIN)/(NGXL-1)
       DO NGX=1,NGXL
          GCX(NGX)=GRMIN+(NGX-1)*GUCLIP(DGX)
       ENDDO
       DGY=DBLE(GZMAX-GZMIN)/(NGYL-1)
       DO NGY=1,NGYL
          GCY(NGY)=GZMIN+(NGY-1)*GUCLIP(DGY)
       ENDDO
       DO NGY=1,NGYL
          ZL=DBLE(GCY(NGY))
          DO NGX=1,NGXL
             RL=DBLE(GCX(NGX))
             CALL PL_MAG_OLD(RL,0.D0,ZL,RHON)
             GCF(NGX,NGY)=GUCLIP(RHON)
          ENDDO
       ENDDO
       CALL CONTP2(GCF,GCX,GCY,NGXL,NGXL,NGYL,0.1,0.1,9,0,1,KA)
    ENDIF

!     ----- cyclotron resonnce -----

    RFL=RFIN(1)
    ID=0
    DO NRAY=2,NRAYMAX
       IF(RFIN(NRAY).NE.RFL) ID=1
    END DO
    IF(ID.EQ.0) THEN
       DGX=DBLE(GRMAX-GRMIN)/(NGXL-1)
       DO NGX=1,NGXL
          GCX(NGX)=GRMIN+(NGX-1)*GUCLIP(DGX)
       ENDDO
       DGY=DBLE(GZMAX-GZMIN)/(NGYL-1)
       DO NGY=1,NGYL
          GCY(NGY)=GZMIN+(NGY-1)*GUCLIP(DGY)
       ENDDO
       DO NGY=1,NGYL
          ZL=DBLE(GCY(NGY))
          DO NGX=1,NGXL
             RL=DBLE(GCX(NGX))
             CALL PL_MAG(RL,0.D0,ZL,mag)
             wce=AEE*mag%BABS/AME
             wwl=2.D0*PI*RFL*1.D6
             GCF(NGX,NGY)=GUCLIP(wwl/wce)
             IF(NGX.EQ.(NGXL+1)/2.AND.NGY.EQ.NGYL) THEN
                WRITE(6,'(A,1P5E12.4)') 'www/wce:',RL,ZL,wwl,wce,wwl/wce
             END IF
          ENDDO
       ENDDO
       CALL SETRGB(0.7,0.7,0.0)
       CALL CONTP2(GCF,GCX,GCY,NGXL,NGXL,NGYL,1.0,1.0,5,0,3,KA)
       CALL SETRGB(0.0,1.0,0.0)
    ENDIF

!  ----- RAY TRAJECTORY -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RL=SQRT(RAYS(1,NSTP,NRAY)**2+RAYS(2,NSTP,NRAY)**2)
          GX(NSTP+1)=GUCLIP(RL)
          GY(NSTP+1)=GUCLIP(RAYS(3,NSTP,NRAY))
       ENDDO
       CALL SETLIN(0,2,7-MOD(NRAY-1,5))
       CALL GPLOTP(GX,GY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO


!   ----- draw deposition profile vs minor radius -----

    DO nray=1,nraymax
       DO nrs=1,nrsmax
          GPY(nrs,nray)=GUCLIP(pwr_nrs_nray(nrs,nray)*pos_nrs(nrs))
       ENDDO
    ENDDO
    DO nrs=1,nrsmax
       GPX(nrs)=GUCLIP(pos_nrs(nrs))
    END DO

    CALL GQSCAL(GUCLIP(RHOGMN),GUCLIP(RHOGMX),GXMIN,GXMAX,GXSTEP)

    CALL GMNMX2(GPY,nrsmax,1,nrsmax,1,1,NRAYMAX,1,GYMIN,GYMAX)
    CALL GQSCAL(GYMIN,GYMAX,GYSMIN,GYSMAX,GYSCAL)
!    GYSMIN=0.0
!    GYSMAX2=0.1 
!    GYSCAL2=2.E-2
!    GYSMAX=MAX(GYSMAX1,GYSMAX2)
!    GYSCAL=MAX(GYSCAL1,GYSCAL2)

    CALL SETLIN(0,2,7)
    CALL MOVE(13.5,16.1)
    CALL TEXT('ABS POWER*RS',12)
    IF(MOD(MDLWRG/2,2).EQ.0) THEN
       CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,0.0,GYSMAX)
       CALL SETLIN(0,2,7)
       CALL GFRAME
       CALL GSCALE(0.0,  GXSTEP,0.0,  GYSCAL,0.1,9)
       CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
       CALL GVALUE(0.0,0.0,0.0,GYSCAL,NGULEN(GYSCAL))
    ELSE
       CALL GDEFIN(13.5,23.5,9.0,16.0,GXMIN,GXMAX,0.0,GYSMAX)
       CALL SETLIN(0,2,7)
       CALL GFRAME
       GXORG=(INT(GXMIN/(2*GXSTEP))+1)*2*GXSTEP
       CALL GSCALE(GXORG,  GXSTEP,0.0,  GYSCAL,0.1,9)
       CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
       CALL GVALUE(0.0,0.0,0.0,GYSCAL,NGULEN(GYSCAL))
    ENDIF

    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,2,7-MOD(NRAY-1,5))
       CALL GPLOTP(GPX,GPY(1:nrsmax,nray),1,nrsmax,1,0,0,0)
    ENDDO

!  ----- draw power flux vs minor radius -----
      
    GZSMIN=0.0
    GZSMAX=1.1
    GZSCAL=0.2

    CALL SETLIN(0,2,7)
    CALL MOVE(13.5,8.1)
    CALL TEXT('POWER FLUX',10)
    IF(MOD(MDLWRG/2,2).EQ.0) THEN
       CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,0.0,GZSMAX)
       CALL SETLIN(0,2,7)
       CALL GFRAME
       CALL GSCALE(0.0,  GXSTEP,0.0,  GZSCAL,0.1,9)
       CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
       CALL GVALUE(0.0,0.0,0.0,GZSCAL,1)
    ELSE
       CALL GDEFIN(13.5,23.5,1.0,8.0,GXMIN,GXMAX,0.0,GZSMAX)
       CALL SETLIN(0,2,7)
       CALL GFRAME
       GXORG=(INT(GXMIN/(2*GXSTEP))+1)*2*GXSTEP
       CALL GSCALE(GXORG,  GXSTEP,0.0,  GZSCAL,0.1,9)
       CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
       CALL GVALUE(0.0,0.0,0.0,GZSCAL,1)
    ENDIF
    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          XL=RAYS(1,NSTP,NRAY)
          YL=RAYS(2,NSTP,NRAY)
          ZL=RAYS(3,NSTP,NRAY)
          CALL PL_MAG_OLD(XL,YL,ZL,RHON)
          GUX(NSTP+1)=GUCLIP(RHON)
          GUY(NSTP+1)=GUCLIP(RAYS(7,NSTP,NRAY))
       ENDDO
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GUX,GUY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)      
    ENDDO

    CALL WRGPRM(1)
    CALL PAGEE
    RETURN
  END SUBROUTINE WRGRF1

!  ***** TOROIDAL TRAJECTORY AND POWER *****

  SUBROUTINE WRGRF2

    USE wrcomm
    USE plprof,ONLY: pl_mag,pl_mag_old
    IMPLICIT NONE
    INTEGER,PARAMETER:: NSUM=201
    INTEGER,PARAMETER:: NGXL=101
    INTEGER,PARAMETER:: NGYL=101
    REAL,ALLOCATABLE:: GLCX(:),GLCY(:),GSCX(:),GSCY(:)
    REAL:: GX(NSTPMAX+1),GY(NSTPMAX+1)
    INTEGER:: NSR,NSRMAX,NRL,NRAY,NSTP
    REAL:: GRMIN,GRMAX,GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GYSTEP,GXORG
    REAL(rkind):: RMAX,RMIN,DTH,TH,XL,YL,ZL,RHON
    REAL:: GUX(NSTPMAX+1),GUY(NSTPMAX+1)
    REAL:: GPX(nrlmax),GPY(nrlmax,nraymax)
    REAL:: GYSMIN,GYSMAX,GYSCAL,GRSMIN,GRSMAX,GRSCAL
    REAL:: GZSMIN,GZSMAX,GZSCAL
    EXTERNAL GMNMX1,GMNMX2,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT
    INTEGER:: NGULEN
    REAL:: GUCLIP

!   --- ray trajectory on toroidal cross section -----

    NSRMAX=101
    ALLOCATE(GLCX(NSRMAX),GLCY(NSRMAX),GSCX(NSRMAX),GSCY(NSRMAX))
    DTH=2*PI/(NSRMAX-1)
    DO NSR=1,NSRMAX
       TH=DTH*(NSR-1)
       GLCX(NSR)=GUCLIP(RMAX*COS(TH))
       GLCY(NSR)=GUCLIP(RMAX*SIN(TH))
       GSCX(NSR)=GUCLIP(RMIN*COS(TH))
       GSCY(NSR)=GUCLIP(RMIN*SIN(TH))
    ENDDO
    GRMIN=GUCLIP(RMIN)
    GRMAX=GUCLIP(RMAX)
    CALL GQSCAL(-GRMAX,GRMAX,GXMIN,GXMAX,GXSTEP)
    CALL GQSCAL(-GRMAX,GRMAX,GYMIN,GYMAX,GYSTEP)

    CALL PAGES
    CALL SETCHS(0.25,0.)
    CALL SETFNT(32)

    CALL SETLIN(0,2,7)
    CALL GDEFIN(1.7,11.7,1.0,11.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME

!    GXORG=(INT(GXMIN/(2*GXSTEP))+1)*2*GXSTEP
    GXORG=0.0
    IF(idebug_wr(92).NE.0) THEN
       WRITE(6,'(A,1P3E12.4)') 'GX:',GXORG,GXSTEP,GYSTEP
    END IF
    CALL GSCALE(GXORG,GXSTEP,0.0,GYSTEP,0.1,9)
    CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,0.0,2*GYSTEP,NGULEN(2*GYSTEP))

    CALL SETLIN(0,2,4)
    CALL GPLOTP(GLCX,GLCY,1,101,1,0,0,0)
    CALL GPLOTP(GSCX,GSCY,1,101,1,0,0,0)

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          GX(NSTP+1)=GUCLIP(RAYS(1,NSTP,NRAY))
          GY(NSTP+1)=GUCLIP(RAYS(2,NSTP,NRAY))
       ENDDO
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GX,GY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

!   ----- draw deposition profile vs minor radius -----

    DO nray=1,nraymax
       DO nrl=1,nrlmax
          GPY(nrl,nray)=GUCLIP(pwr_nrl_nray(nrl,nray))
       ENDDO
    ENDDO
    DO nrl=1,nrlmax
       GPX(nrl)=GUCLIP(pos_nrl(nrl))
    END DO

    CALL GMNMX2(GPY,NRLMAX,1,NRLMAX,1,1,NRAYMAX,1,GYMIN,GYMAX)
    CALL GQSCAL(GYMIN,GYMAX,GYSMIN,GYSMAX,GYSCAL)
    CALL GQSCAL(GRMIN,GRMAX,GRSMIN,GRSMAX,GRSCAL)
    
    CALL SETLIN(0,2,7)
    CALL MOVE(13.5,16.1)
    CALL TEXT('ABS POWER vs R',10)
    CALL GDEFIN(13.5,23.5,9.0,16.0,GRSMIN,GRSMAX,0.0,GYSMAX)
    CALL GFRAME
    CALL GSCALE(GRSMIN,  GRSCAL,0.0,GYSCAL,0.1,9)
    CALL GVALUE(GRSMIN,2*GRSCAL,0.0,   0.0,NGULEN(2*GRSCAL))
    CALL GVALUE(GRSMIN,     0.0,0.0,GYSCAL,NGULEN(GYSCAL))
    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,2,7-MOD(NRAY-1,5))
       CALL GPLOTP(GPX,GPY(1:nrlmax,nray),1,nrlmax,1,0,0,0)
    ENDDO

!     ----- draw power flux vs major radius ---

    GZSMIN=0.0
    GZSMAX=1.1
    GZSCAL=0.2

    CALL SETLIN(0,2,7)
    CALL MOVE(13.5,8.1)
    CALL TEXT('POWER FLUX',10)
    CALL GDEFIN(13.5,23.5,1.0,8.0,GRSMIN,GRSMAX,0.0,GZSMAX)
    CALL GFRAME
    CALL GSCALE(GRSMIN,  GRSCAL,0.0,GZSCAL,0.1,9)
    CALL GVALUE(GRSMIN,2*GRSCAL,0.0,   0.0,NGULEN(2*GRSCAL))
    CALL GVALUE(0.0,0.0,0.0,GZSCAL,1)
    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          XL=RAYS(1,NSTP,NRAY)
          YL=RAYS(2,NSTP,NRAY)
          ZL=RAYS(3,NSTP,NRAY)
          CALL PL_MAG_OLD(XL,YL,ZL,RHON)
          GUX(NSTP+1)=GUCLIP(SQRT(XL**2+YL**2))
          GUY(NSTP+1)=GUCLIP(RAYS(7,NSTP,NRAY))
!          IF(MOD(NSTP,100).EQ.0) &
!               WRITE(6,'(A,I5,1P5E12.4)') &
!               'PF:',NSTP,XL,YL,ZL,RHON,RAYS(7,NSTP,NRAY)
       ENDDO
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GUX,GUY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)      
    ENDDO

    CALL WRGPRM(1)
    CALL PAGEE

    RETURN
  END SUBROUTINE WRGRF2

!     ***** Parallel and perpendicular wave number *****

  SUBROUTINE WRGRF3

    USE wrcomm
    USE plprof,ONLY: pl_mag_old,pl_rzsu
    USE wrsub,ONLY: wrcalk
    IMPLICIT NONE
    INTEGER,PARAMETER:: NSUM=201
    REAL:: GKX(NSTPMAX+1,NRAYMAX)
    REAL:: GKY1(NSTPMAX+1,NRAYMAX),GKY2(NSTPMAX+1,NRAYMAX)
    INTEGER:: NRAY,NSTP,NSU
    REAL(rkind):: XL,YL,ZL,RHON,RKPARA,RKPERP,RMIN,RMAX
    REAL:: GYMIN1,GYMAX1,GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP,GYORG
    REAL:: GXMIn,GXMAX,GXSTEP
    EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT
    INTEGER:: NGULEN
    REAL:: GUCLIP

    CALL PAGES
    CALL SETFNT(32)
    CALL SETCHS(0.25,0.)
    CALL SETLIN(0,2,7)

!     ----- minor radius dependence of wave number -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          XL=RAYS(1,NSTP,NRAY)
          YL=RAYS(2,NSTP,NRAY)
          ZL=RAYS(3,NSTP,NRAY)
          CALL PL_MAG_OLD(XL,YL,ZL,RHON)
          CALL WRCALK(NSTP,NRAY,RKPARA,RKPERP)
          GKX(NSTP+1,NRAY)=GUCLIP(RHON)
          GKY1(NSTP+1,NRAY)=GUCLIP(RKPARA*VC &
               /(2*PI*RAYIN(1,NRAY)*1.D6))
          GKY2(NSTP+1,NRAY)=GUCLIP(RKPERP)
       ENDDO
    ENDDO

    CALL GQSCAL(0.0,1.0,GXMIN,GXMAX,GXSTEP)
    NRAY=1
       CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN1,GYMAX1)
       CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN2,GYMAX2)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
       CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
       GYMIN2=MIN(GYMIN2,GYMIN)
       GYMAX2=MAX(GYMAX2,GYMAX)
    ENDDO
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(1.7,16.1)
    CALL TEXT('n-para vs rhon',6)
    CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

    CALL GQSCAL(GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(1.7,8.1)
    CALL TEXT('k-perp vs rhon',6)
    CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

!     ----- major radius dependence of wave number -----

    RMIN=RSU(1)
    RMAX=RSU(1)
    DO NSU=2,NSUMAX
       RMIN=MIN(RMIN,RSU(NSU))
       RMAX=MAX(RMAX,RSU(NSU))
    END DO

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          XL=RAYS(1,NSTP,NRAY)
          YL=RAYS(2,NSTP,NRAY)
          ZL=RAYS(3,NSTP,NRAY)
          CALL WRCALK(NSTP,NRAY,RKPARA,RKPERP)
          GKX( NSTP+1,NRAY)=GUCLIP(SQRT(XL**2+YL**2))
          GKY1(NSTP+1,NRAY)=GUCLIP(RKPARA*VC/(2*PI*RAYIN(1,NRAY)*1.D6))
          GKY2(NSTP+1,NRAY)=GUCLIP(RKPERP)
       ENDDO
    ENDDO

    CALL GQSCAL(GUCLIP(RMIN),GUCLIP(RMAX),GXMIN,GXMAX,GXSTEP)

    NRAY=1
       CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN1,GYMAX1)
       CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN2,GYMAX2)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
       CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
       GYMIN2=MIN(GYMIN2,GYMIN)
       GYMAX2=MAX(GYMAX2,GYMAX)
    ENDDO

    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(13.7,16.1)
    CALL TEXT('n-para vs R',11)
    CALL GDEFIN(13.7,23.7,9.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(GXMIN,GXSTEP,GYMIN,GYSTEP,0.1,9)
    CALL GVALUE(GXMIN,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

    CALL GQSCAL(GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(13.7,8.1)
    CALL TEXT('k-perp vs R',11)
    CALL GDEFIN(13.7,23.7,1.0,8.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(GXMIN,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(GXMIN,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

    CALL WRGPRM(1)
    CALL PAGEE

    RETURN
  END SUBROUTINE WRGRF3

!     ***** RADIAL DEPENDENCE 2 *****

  SUBROUTINE WRGRF4

    USE wrcomm
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE

    REAL:: GKX(NSTPMAX+1,NRAYMAX)
    REAL:: GKY(NSTPMAX+1,NRAYMAX)
    INTEGER:: NRAY,NSTP
    REAL:: GXMIN,GXMAX,GXSTEP
    REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP,GYORG
    REAL(rkind):: XL,YL,ZL,RHON,RL,RKR,RKZ,RKPH,RNPHI
    EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT
    INTEGER:: NGULEN
    REAL:: GUCLIP

    CALL PAGES
    CALL SETFNT(32)
    CALL SETCHS(0.25,0.)

!     ----- X AXIS -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          XL=RAYS(1,NSTP,NRAY)
          YL=RAYS(2,NSTP,NRAY)
          ZL=RAYS(3,NSTP,NRAY)
          CALL PL_MAG_OLD(XL,YL,ZL,RHON)
          GKX( NSTP+1,NRAY)=GUCLIP(RHON)
       ENDDO
    ENDDO
    CALL GQSCAL(0.0,1.0,GXMIN,GXMAX,GXSTEP)

!     ----- Fig.1 -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RL=SQRT(RAYS(1,NSTP,NRAY)**2+RAYS(2,NSTP,NRAY)**2)
          RKR  =( RAYS(4,NSTP,NRAY)*RAYS(1,NSTP,NRAY) &
                 +RAYS(5,NSTP,NRAY)*RAYS(2,NSTP,NRAY))/RL
          GKY( NSTP+1,NRAY)=GUCLIP(RKR)
       ENDDO
    ENDDO
    NRAY=1
    CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN1,GYMAX1)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
    ENDDO
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(1.7,16.1)
    CALL TEXT('k-R',3)
    CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

!     ----- Fig.2 -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RKZ  =  RAYS(6,NSTP,NRAY)
          GKY( NSTP+1,NRAY)=GUCLIP(RKZ)
       ENDDO
    ENDDO
    NRAY=1
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN1,GYMAX1)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
    ENDDO
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(13.7,16.1)
    CALL TEXT('k-Z',3)
    CALL GDEFIN(13.7,23.7,9.0,16.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

!     ----- Fig.3 -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RL=SQRT(RAYS(1,NSTP,NRAY)**2+RAYS(2,NSTP,NRAY)**2)
          RKPH=(-RAYS(4,NSTP,NRAY)*RAYS(2,NSTP,NRAY) &
               +RAYS(5,NSTP,NRAY)*RAYS(1,NSTP,NRAY))/RL
          GKY( NSTP+1,NRAY)=GUCLIP(RKPH)
       ENDDO
    ENDDO
    NRAY=1
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN1,GYMAX1)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
    ENDDO
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(1.7,8.1)
    CALL TEXT('k-phi',5)
    CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

!     ----- Fig.4 -----

    DO NRAY=1,NRAYMAX
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RNPHI=(-RAYS(4,NSTP,NRAY)*RAYS(2,NSTP,NRAY) &
                 +RAYS(5,NSTP,NRAY)*RAYS(1,NSTP,NRAY))
          GKY( NSTP+1,NRAY)=GUCLIP(RNPHI)
       ENDDO
    ENDDO
    NRAY=1
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN1,GYMAX1)
    DO NRAY=2,NRAYMAX
       CALL GMNMX1(GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
    ENDDO
    CALL GQSCAL(GYMIN1-1.0,GYMAX1+1.0,GYMIN,GYMAX,GYSTEP)

    CALL SETLIN(0,2,7)
    CALL MOVE(13.7,8.1)
    CALL TEXT('R*k-phi',7)
    CALL GDEFIN(13.7,23.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
    CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
    CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
    DO NRAY=1,NRAYMAX
       CALL SETLIN(0,0,7-MOD(NRAY-1,5))
       CALL GPLOTP(GKX(1,NRAY),GKY(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
    ENDDO

    CALL WRGPRM(1)
    CALL PAGEE
    RETURN
  END SUBROUTINE WRGRF4

!     ***** wave field polarization *****

  SUBROUTINE WRGRF5

    USE wrcomm
    USE wrsub,ONLY: wrcalep
    USE libgrf
    IMPLICIT NONE
    INTEGER:: NRAY,NSTP,NSTP_PLMAX,NSTP1,NSTP2,NSTP3,NSTP4,i
    REAL(rkind):: PLMAX,PL,err,errmax
    COMPLEX(rkind):: cepola(3),cenorm(3)
    INTEGER:: nxmax,nymax,ntmax(3*nraymax)
    REAL(rkind):: gx(2),gy(2),gepola(2,nstpmax+1,3*nraymax)
    REAL(rkind):: LINE_RGB(3,3)
    DATA LINE_RGB/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/
    EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT

    DO NRAY=1,NRAYMAX
       NSTP_PLMAX=0
       PLMAX=RAYS(8,0,NRAY)
       DO NSTP=1,NSTPMAX_NRAY(NRAY)
          PL=RAYS(8,NSTP,NRAY)
          IF(PL.GT.PLMAX) THEN
             PLMAX=PL
             NSTP_PLMAX=NSTP
          END IF
       END DO
       NSTP1=0
       NSTP2=NSTP_PLMAX/2
       NSTP3=NSTP_PLMAX
       NSTP4=NSTPMAX_NRAY(NRAY)
       
       CALL wrcalep(nstp1,nray,cepola,cenorm,err)
       IF(idebug_wr(90).NE.0) THEN
          WRITE(6,'(A,I4,I8,1PE12.4)') 'nray,nstp,R=',nray,nstp1, &
               SQRT(rays(1,nstp1,nray)**2+rays(2,nstp1,nray)**2)
          WRITE(6,'(A,1P6E12.4)')    'cexyz=',(cepola(i),i=1,3)
          WRITE(6,'(A,1P6E12.4)')    'ceoxp=',(cenorm(i),i=1,3)
       END IF

       CALL wrcalep(nstp2,nray,cepola,cenorm,err)
       IF(idebug_wr(90).NE.0) THEN
          WRITE(6,'(A,I4,I8,1PE12.4)') 'nray,nstp,R=',nray,nstp2, &
               SQRT(rays(1,nstp2,nray)**2+rays(2,nstp2,nray)**2)
          WRITE(6,'(A,1P6E12.4)')    'cexyz=',(cepola(i),i=1,3)
          WRITE(6,'(A,1P6E12.4)')    'ceoxp=',(cenorm(i),i=1,3)
       END IF
       
       CALL wrcalep(nstp3,nray,cepola,cenorm,err)
       IF(idebug_wr(90).NE.0) THEN
          WRITE(6,'(A,I4,I8,1PE12.4)') 'nray,nstp,R=',nray,nstp3, &
               SQRT(rays(1,nstp3,nray)**2+rays(2,nstp3,nray)**2)
          WRITE(6,'(A,1P6E12.4)')    'cexyz=',(cepola(i),i=1,3)
          WRITE(6,'(A,1P6E12.4)')    'ceoxp=',(cenorm(i),i=1,3)
       END IF
       
       CALL wrcalep(nstp4,nray,cepola,cenorm,err)
       IF(idebug_wr(90).NE.0) THEN
          WRITE(6,'(A,I4,I8,1PE12.4)') 'nray,nstp,R=',nray,nstp4, &
               SQRT(rays(1,nstp4,nray)**2+rays(2,nstp4,nray)**2)
          WRITE(6,'(A,1P6E12.4)')    'cexyz=',(cepola(i),i=1,3)
          WRITE(6,'(A,1P6E12.4)')    'ceoxp=',(cenorm(i),i=1,3)
       END IF
    END DO

    CALL PAGES
    CALL SETFNT(32)
    CALL SETCHS(0.25,0.)
    CALL SETLIN(0,2,7)

    NXMAX=2
    GX(1)=-1.D0
    GX(2)= 1.D0
    NYMAX=2
    GY(1)=-1.D0
    GY(2)= 1.D0

    DO NRAY=1,NRAYMAX
       errmax=0.D0
       NTMAX(3*(nray-1)+1)=NSTPMAX_NRAY(NRAY)+1
       NTMAX(3*(nray-1)+2)=NSTPMAX_NRAY(NRAY)+1
       NTMAX(3*(nray-1)+3)=NSTPMAX_NRAY(NRAY)+1
       IF(idebug_wr(90).NE.0) THEN
          WRITE(6,'(A,3I8)') &
               '## nray,ntmax:',nray,nstpmax+1,NSTPMAX_NRAY(NRAY)+1
       END IF
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          CALL wrcalep(nstp,nray,cepola,cenorm,err)
          gepola(1,nstp+1,3*(nray-1)+1)= REAL(cenorm(1))
          gepola(2,nstp+1,3*(nray-1)+1)=AIMAG(cenorm(1))
          gepola(1,nstp+1,3*(nray-1)+2)= REAL(cenorm(2))
          gepola(2,nstp+1,3*(nray-1)+2)=AIMAG(cenorm(2))
          gepola(1,nstp+1,3*(nray-1)+3)= REAL(cenorm(3))
          gepola(2,nstp+1,3*(nray-1)+3)=AIMAG(cenorm(3))
          errmax=MAX(err,errmax)
       END DO
    END DO
    IF(idebug_wr(90).NE.0) THEN
       WRITE(6,'(A,I4,1PE12.4)') 'nray,errmax=',nray,errmax
    END IF

    CALL grdxy(0,gepola,2,nstpmax+1,ntmax,3*nraymax, &
         '@Polarization@',NLMAX=3,XSCALE_ZERO=0,YSCALE_ZERO=0, &
         LINE_RGB=LINE_RGB)
    CALL PAGEE
  END SUBROUTINE WRGRF5
    
!     ***** relativistic cyclotron resonance condition *****

  SUBROUTINE WRGRF6

    USE wrcomm
    USE wrsub,ONLY: wrcalk
    USE plprof
    USE libgrf
    IMPLICIT NONE
    INTEGER,PARAMETER:: nang_max=100
    TYPE(pl_mag_type):: mag
    TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
    INTEGER:: nray,nres,nstp,nc,nang,line_pat,ngid,ierr
    INTEGER:: nstp_nres(nres_max)
    REAL:: rgb_nres(3,nres_max)
    REAL(rkind):: pvtmax,dang,omega,xl,yl,zl,rhon,babs,omegace
    REAL(rkind):: ptpr_e,ptpp_e,rkpara,rkperp,rnpara,pc_org,temp
    REAL(rkind):: pc_para,pc_perp,pte,vte,pvt_org,pvt_para,pvt_perp,ang
    CHARACTER(LEN=46):: title
    EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT,MOVE2D,DRAW2D
    REAL:: GUCLIP

    pvtmax=6.D0
    dang=PI/nang_max
    
    CALL pages
    DO nray=1,MIN(nraymax,25)
       SELECT CASE(nraymax)
       CASE(1)
          ngid=0
       CASE(2:4)
          ngid=nray    ! ngid=1:4
       CASE(5:9)
          ngid=nray+4  ! ngid=9:13  (5:13)
       CASE(10:16)
          ngid=nray+13 ! ngid=23:29 (14:29)
       CASE(17:25)
          ngid=nray+29 ! ngid=46:54 (30:54)
       END SELECT

       CALL setup_nres(nray,nstp_nres,rgb_nres,ierr)
       WRITE(title,'(A,ES12.4,A,F8.4,A,F8.4,A)') &
            '@RF:',RAYIN(1,nray), &
            ' angPH:',RAYIN(7,nray),' angZ:',RAYIN(6,nray),'@'
       CALL grd2d_frame_start(ngid,-pvtmax,pvtmax,0.D0,pvtmax, &
            title,ASPECT=0.5D0,NOINFO=1)

       omega=2*PI*RAYIN(1,NRAY)*1.D6

       DO nres=1,nres_max
          CALL SETRGB(rgb_nres(1,nres),rgb_nres(2,nres),rgb_nres(3,nres))
          nstp=nstp_nres(nres)
          IF(nstp.LE.0.OR.nstp.GT.nstpmax_nray(nray)) CYCLE
          IF(idebug_wr(91).NE.0) THEN
             WRITE(6,'(A,3I8)') 'nres,nstp,nres_max=',nres,nstp,nres_max
          END IF
          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          zl=rays(3,nstp,nray)
          CALL pl_mag(xl,yl,zl,mag)
          rhon=mag%rhon
          babs=mag%babs
          omegace=AEE*babs/AME
          CALL pl_prof(rhon,plf)
          ptpr_e=plf(1)%rtpr
          ptpp_e=plf(1)%rtpp
          pte=(ptpr_e+2.D0*ptpp_e)/3.D0
          vte=SQRT(pte*AEE*1.D3/(AME*VC**2))
          CALL wrcalk(nstp,nray,rkpara,rkperp)
          rnpara=rkpara*VC/omega
          IF(ABS(rnpara).LT.1.D0) THEN
             DO nc=0,ncmax(1)  ! for electron only
                SELECT CASE(ABS(nc))
                CASE(1)
                   line_pat=0
                CASE(2)
                   line_pat=2
                CASE(3)
                   line_pat=4
                CASE DEFAULT
                   line_pat=6
                END SELECT
                pc_org=rnpara/(1.D0-rnpara**2)*nc*omegace/omega
                temp=((nc*omegace/omega)**2-(1.D0-rnpara**2))/(1-rnpara**2)**2
                IF(temp.GT.0.D0) THEN
                   pc_para= SQRT(temp)
                   pc_perp=pc_para*SQRT(1.D0-rnpara**2)
                   pvt_org=pc_org/vte
                   pvt_para=pc_para/vte
                   pvt_perp=pc_perp/vte
                   ang=0.D0
                   CALL MOVEPT2D(GUCLIP(pvt_org+pvt_para*cos(ang)), &
                                 GUCLIP(pvt_perp*sin(ang)),line_pat)
                   DO nang=2,nang_max
                      ang=(nang-1)*dang
                      CALL DRAWPT2D(GUCLIP(pvt_org+pvt_para*cos(ang)), &
                                    GUCLIP(pvt_perp*sin(ang)))
                   END DO
                END IF
             END DO
          END IF
       END DO
       CALL grd2d_frame_end
    END DO

    CALL pagee
  END SUBROUTINE WRGRF6

  SUBROUTINE setup_nres(nray,nstp_nres,rgb_nres,ierr)
    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nray
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nstp_nres(nres_max)
    REAL:: rgb_nres(3,nres_max)
    REAL(rkind):: level_nres(nres_max),level_nstp(nstpmax)
    REAL(rkind):: pwmax,pw,rlen,del_nres,pf_init,pf_last,pf_abs
    INTEGER:: nstp_pwmax,nstp,nres
    REAL(rkind):: level,factor
    REAL:: GUCLIP

    ierr=0

    ! --- range of power flux ---

    pf_init=rays(7,0,nray)
    pf_last=rays(7,nstpmax_nray(nray),nray)
    pf_abs=pf_init-pf_last
    IF(pf_abs.LE.0.D0) THEN
       ierr=1
       RETURN
    END IF

    ! --- power abs density ---

    SELECT CASE(nres_type)
    CASE(0)

       ! --- find pabs max ---
       
       pwmax=rays(8,0,nray)
       nstp_pwmax=0
       DO nstp=0,nstpmax_nray(nray)
          pw=rays(8,nstp,nray)
          IF(pw.GT.pwmax) THEN
             pwmax=pw
             nstp_pwmax=nstp
          END IF
       END DO

       ! --- set level for nstp

       DO nstp=1,nstpmax_nray(nray)
          IF(nstp.LT.nstp_pwmax) THEN
             level_nstp(nstp)=0.5D0*rays(8,nstp,nray)/pwmax
          ELSE IF(nstp.GT.nstp_pwmax) THEN
             level_nstp(nstp)=1.D0-0.5D0*rays(8,nstp,nray)/pwmax
          ELSE
             level_nstp(nstp)=0.5D0
          END IF
       END DO
       
       ! --- set levels for nres points ---
       
       del_nres=1.D0/(nres_max+1)
       DO nres=1,nres_max
          DO nstp=1,nstpmax_nray(nray)
             IF(level_nstp(nstp).GT.del_nres*nres) THEN
                nstp_nres(nres)=nstp
                level_nres(nres)=level_nstp(nstp)
                EXIT
             END IF
          END DO
       END DO

    ! --- power flux level ---

    CASE(1)

       del_nres=rays(7,0,nray)/DBLE(nres_max+1)
       DO nres=1,nres_max
          DO nstp=1,nstpmax_nray(nray)
             level=rays(7,0,nray)-rays(7,nstp,nray)
             IF(level.GT.del_nres*nres) THEN
                nstp_nres(nres)=nstp
                level_nres(nres)=level
                EXIT
             END IF
          END DO
       END DO

    ! --- trajectory length ---

    CASE(2)   

       rlen=0.D0
       DO nstp=1,nstpmax_nray(nray)
          rlen=rlen+SQRT((rays(1,nstp,nray)-rays(1,nstp-1,nray))**2 &
                        +(rays(2,nstp,nray)-rays(2,nstp-1,nray))**2 &
                        +(rays(3,nstp,nray)-rays(3,nstp-1,nray))**2)
          level_nstp(nstp)=rlen
       END DO
       DO nstp=1,nstpmax_nray(nray)
          level_nstp(nstp)=level_nstp(nstp)/rlen
       END DO
       
       del_nres=1.D0/nres_max
       DO nres=1,nres_max
          DO nstp=1,nstpmax_nray(nray)
             IF(level_nstp(nstp).GE.del_nres*nres) THEN
                nstp_nres(nres)=nstp
                level_nres(nres)=level_nstp(nstp)
                EXIT
             END IF
          END DO
       END DO
    END SELECT

    IF(idebug_wr(91).NE.0) THEN
       WRITE(6,'(A,I8)') 'nres_max=',nres_max
    END IF
    
    DO nres=1,nres_max
       level=GUCLIP(level_nres(nres))
       
       IF(idebug_wr(91).NE.0) THEN
          WRITE(6,'(A,2I4,I8,3ES12.4)') &
               'nray,nres,nstp,level,pf,pw=',nray,nres,nstp,level, &
               rays(7,nstp,nray),rays(8,nstp,nray)
       end IF
       
       IF(level.LT.0.1) THEN
          factor=level/0.1
          rgb_nres(1,nres)=0.0
          rgb_nres(2,nres)=0.0
          rgb_nres(3,nres)=GUCLIP(factor)
       ELSE IF(level.LT.0.5) THEN
          factor=(level-0.1)/0.4
          rgb_nres(1,nres)=GUCLIP(factor)
          rgb_nres(2,nres)=0.0
          rgb_nres(3,nres)=1.0-GUCLIP(factor)
       ELSE IF(level.LE.0.9D0) THEN
          factor=(level-0.5)/0.4
          rgb_nres(1,nres)=1.0-GUCLIP(factor)
          rgb_nres(2,nres)=GUCLIP(factor)
          rgb_nres(3,nres)=0.0
       ELSE
          factor=(GUCLIP(level)-0.9)/0.1
          rgb_nres(1,nres)=0.0
          rgb_nres(2,nres)=1.0-GUCLIP(factor)
          rgb_nres(3,nres)=0.0
       END IF
    END DO
    IF(idebug_wr(91).NE.0) THEN
       WRITE(6,'(A,2I8)') 'nres:',nray,nstpmax_nray(nray)
       DO nres=1,nres_max
          WRITE(6,'(A,2I8,3ES12.4)') 'nstp:',nstp,nstp_nres(nres), &
               rgb_nres(1,nres),rgb_nres(2,nres),rgb_nres(3,nres)
       END DO
    END IF
    RETURN
  END SUBROUTINE setup_nres
       
!     ***** relativistic cyclotron resonance condition *****

  SUBROUTINE WRGRF7(ns,nmax,rhona)

    USE wrcomm
    USE wrsub,ONLY: wrcalk
    USE plprof
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ns,nmax
    REAL(rkind),INTENT(IN):: rhona(nmax)
    TYPE(pl_mag_type):: mag
    CHARACTER(LEN=20):: title
    INTEGER,parameter:: nang_max=50
    INTEGER:: n,ngid,nray,i,nstp,nc,line_pat,nang,ierr
    REAL(rkind):: dang,rhon,pt0,vt0c,omega,xl,yl,zl,omegac
    REAL(rkind):: rkpara,rkperp,rnpara,pc_org,pc_rad2,pc_rad,pv_org,pv_rad
    REAL(rkind):: ang,x,y
    INTEGER:: nstp_nray(2,nraymax)
    EXTERNAL PAGES,PAGEE,MOVE2D,DRAW2D

    IF(nmax.LE.4) THEN
       WRITE(6,'(A,2I6,4ES12.4)') 'ns,nmax,rhona=',ns,nmax,(rhona(n),n=1,nmax)
    ELSE
       WRITE(6,'(A,2I6,4ES12.4)') 'ns,nmax,rhona=',ns,nmax,(rhona(n),n=1,4)
       WRITE(6,'(14X   5ES12.4)')                          (rhona(n),n=5,nmax)
    END IF
       
    dang=PI/nang_max
    CALL pages

    DO n=1,nmax
       SELECT CASE(nmax)
       CASE(1)
          ngid=0
       CASE(2:4)
          ngid=n    ! ngid=1:4
       CASE(5:9)
          ngid=n+4  ! ngid=9:13  (5:13)
       CASE(10:16)
          ngid=n+13 ! ngid=23:29 (14:29)
       CASE(17:25)
          ngid=n+29 ! ngid=46:54 (30:54)
       END SELECT

       rhon=rhona(n)
       CALL setup_nray(rhon,nstp_nray,ierr)
       WRITE(title,'(A,F8.4,A)') '@Res: rhon=',rhon,'@'
       
       CALL grd2d_frame_start(ngid,-pmax_dp(ns),pmax_dp(ns),0.D0,pmax_dp(ns), &
                              title,ASPECT=0.5D0,NOINFO=1)

       pt0=(ptpr(ns)+2.D0*ptpp(ns))/3.D0          ! axis temperature [J]
       vt0c=SQRT(pt0*AEE*1.D3/(PA(ns)*AMP*VC**2))  ! vt0/c
       DO nray=1,nraymax
!          CALL SETRGBDA(line_rgb_nlmax(1:3,MOD(2*(nray-1)+(i-1),nlmax_p)+1))
          CALL SETLIN(0,2,7-MOD(nray-1,5))
          omega=2*PI*RAYIN(1,nray)*1.D6
          DO i=1,2
             nstp=nstp_nray(i,nray)
             IF(nstp.NE.0) THEN
                xl=rays(1,nstp,nray)
                yl=rays(2,nstp,nray)
                zl=rays(3,nstp,nray)
                CALL pl_mag(xl,yl,zl,mag)
                omegac=PZ(ns)*AEE*mag%babs/(PA(ns)*AMP)
                CALL wrcalk(nstp,nray,rkpara,rkperp)
                rnpara=rkpara*VC/omega
                IF(ABS(rnpara).LT.1.D0) THEN
                   DO nc=ncmin(ns),ncmax(ns)
                      SELECT CASE(ABS(nc))
                      CASE(1)
                         line_pat=0
                      CASE(2)
                         line_pat=2
                      CASE(3)
                         line_pat=4
                      CASE DEFAULT
                         line_pat=6
                      END SELECT
                      pc_org=rnpara/(1.D0-rnpara**2)*nc*omegac/omega
                      IF(nc*omegac/omega+rnpara*pc_org.GT.1.D0) THEN
                         pc_rad2=((nc*omegac/omega)**2-(1.D0-rnpara**2)) &
                                 /(1-rnpara**2)**2
                         IF(pc_rad2.GT.0.D0) THEN
                            pc_rad= SQRT(pc_rad2)
                            pv_org=pc_org/vt0c
                            pv_rad=pc_rad/vt0c
                            ang=0.D0
                            x=pv_org+                     pv_rad*cos(ang)
                            y=       SQRT(1.D0-rnpara**2)*pv_rad*sin(ang)
                            CALL MOVEPT2D(gdclip(x),gdclip(y),line_pat)
                            DO nang=2,nang_max+1
                               ang=(nang-1)*dang
                               x=pv_org+                     pv_rad*cos(ang)
                               y=       SQRT(1.D0-rnpara**2)*pv_rad*sin(ang)
                               CALL DRAWPT2D(gdclip(x),gdclip(y))
                            END DO
                         END IF ! pc_rad2>0
                      END IF ! cyclotron resonance condition
                   END DO ! nc
                END IF ! anpara
             END IF ! nstp
          END DO ! i
       END DO ! nray
       CALL grd2d_frame_end
    END DO ! n

    CALL pagee
  END SUBROUTINE WRGRF7

  SUBROUTINE setup_nray(rhon,nstp_nray,ierr)
    USE wrcomm
    USE plprof
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    INTEGER,INTENT(OUT):: ierr
    INTEGER,INTENT(OUT):: nstp_nray(2,nraymax)
    REAL(rkind):: pwmax,pw,x,y,z,rhon1,rhon2
    INTEGER:: nstp_pwmax,nstp,nray,mode
    TYPE(pl_mag_type):: mag

    ierr=0

    DO nray=1,nraymax

       ! --- pabs max for rhon=0.D0 ---

       IF(rhon.LE.0.D0) THEN
          
          pwmax=rays(8,0,nray)
          nstp_pwmax=0
          DO nstp=1,nstpmax_nray(nray)
             pw=rays(8,nstp,nray)
             IF(pw.GT.pwmax) THEN
                pwmax=pw
                nstp_pwmax=nstp
             END IF
          END DO
          nstp_nray(1,nray)=nstp_pwmax
          nstp_nray(2,nray)=nstp_pwmax
          mode=1

       ! --- rhon surface (passin max twice) ---

       else
          mode=0
          nstp=0
          nstp_nray(1,nray)=0
          nstp_nray(2,nray)=0
          x=rays(1,nstp,nray)
          y=rays(2,nstp,nray)
          z=rays(3,nstp,nray)
          CALL pl_mag(x,y,z,mag)
          rhon1=mag%rhon
          DO nstp=1,nstpmax_nray(nray)
             x=rays(1,nstp,nray)
             y=rays(2,nstp,nray)
             z=rays(3,nstp,nray)
             CALL pl_mag(x,y,z,mag)
             rhon2=mag%rhon
             SELECT CASE(mode)
             CASE(0)
                IF(rhon1.LT.rhon) THEN
                   mode=2
                ELSE
                   IF(rhon2.LT.rhon) THEN
                      IF(rhon-rhon2.GT.rhon1-rhon) THEN
                         nstp_nray(1,nray)=nstp-1
                      ELSE
                         nstp_nray(1,nray)=nstp
                      END IF
                      mode=1
                   END IF
                END IF
             CASE(1,2)
                IF(rhon2.GT.rhon) THEN
                   IF(rhon2-rhon.GT.rhon-rhon1) THEN
                      nstp_nray(2,nray)=nstp-1
                   ELSE
                      nstp_nray(2,nray)=nstp
                   END IF
                   IF(mode.EQ.1) THEN
                      mode=3
                   END IF
                   EXIT
                END IF
             END SELECT
             rhon1=rhon2
          END DO
       END IF
    END DO
  END SUBROUTINE setup_nray
       
!     ***** BEAM AND POWER *****

  SUBROUTINE WRGRFB1

    USE wrcomm
    USE libspf,ONLY: erf0
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL:: GPX(NSTPMAX+1),GPY(NSTPMAX+1,NRAYMAX)
    REAL:: GKX(NSTPMAX+1,NRAYMAX)
    REAL:: GKY(NSTPMAX+1,NRAYMAX)
    REAL(rkind):: FRLRO1(0:nrsmax+1),FZLRO1(0:nrsmax+1)
    REAL(rkind):: FRLRO2(0:nrsmax+1),FZLRO2(0:nrsmax+1)
    REAL:: GPAY(NSTPMAX+1,NRAYMAX)
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
         RLMA1,ZLMA1,FASSX1,FASSZ1,DELP1, &
         RLMA2,ZLMA2,FASSX2,FASSZ2,DELP2
    INTEGER:: NRAY,NSTP,NRDIV,NRS1,NDR,NR,NRS2
    REAL:: GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP
    REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP
    REAL:: GYMINA,GYMAXA,GYMINB,GYMAXB,GYSMIN,GYSMAX1,GYSCAL1
    REAL:: GYSMAX2,GYSMAX,GYSCAL,GYSCAL2
    REAL(rkind):: XL,YL,ZL,RHON,DRHO,RHON1,RHON2,SDR,DELPWR
    REAL(rkind):: DRAD,RLA,RLAD,ZLA,ZLAD,DRL,DZL,DRLS,DZLS,DABSS
    REAL(rkind):: RLDMAX1,ZLDMAX1,PHIL,BR,BZ,BPHI
    REAL(rkind):: RHOMIN1,RHOMAX1,RL1,ZL1,RLD1,ZLD1,RHOD1,DRLSN
    REAL(rkind):: SIGMA,DRLT,DZLT,DABST,DRLTN,DZLTN,RLDMAX2,ZLDMAX2
    REAL(rkind):: RHOMAX2,RHOMIN2,RL2,ZL2,RLD2,ZLD2,RHOD2,DZLSN
    INTEGER:: NRSMIN1,NRSMAX1,NDBRD1,NABSM1,NRSA1,NCS,NRSMAX2,NRSMIN2
    INTEGER:: NDBRD2,NABSM2,NRSA2,NRSDA1,NRSDA2,NCS2
    EXTERNAL GMNMX1,GMNMX2,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
    EXTERNAL MOVE,TEXT,SETLNW,NUMBD
    EXTERNAL GETRZ
    REAL:: GUCLIP

      ALLOCATE(RLMA1(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(RLMA2(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(ZLMA1(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(ZLMA2(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(FASSX1(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(FASSX2(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(FASSZ1(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(FASSZ2(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(DELP1(0:NSTPMAX+1,0:nrsmax))
      ALLOCATE(DELP2(0:NSTPMAX+1,0:nrsmax))

      CALL PAGES
      CALL SETFNT(32)
      CALL SETCHS(0.25,0.)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!     ----- X AXIS -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
            GKX( NSTP+1,NRAY)=GUCLIP(RHON)
         ENDDO
      ENDDO

      CALL GMNMX1(GKX,1,NSTPMAX_NRAY(1),1,GXMIN1,GXMAX1)
      CALL GQSCAL(0.0,GXMAX1,GXMIN,GXMAX,GXSTEP)

!     ----- Fig.1  R-CURV-----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY( NSTP+1,1)=GUCLIP(RAYB(21,NSTP))
            GKY( NSTP+1,2)=GUCLIP(RAYB(22,NSTP))
         ENDDO
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
!     CALL SETLNW(0.035)
      CALL SETLIN(0,2,7)
      CALL MOVE(1.7,16.1)
      CALL TEXT('R-curv',6)
      CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,-6.0,6.0)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,2.0,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,2.0,1)
      CALL SETRGB(0.0,0.0,1.0)
      CALL GPLOTP(GKX,GKY(1,1),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(1.0,0.0,0.0)
      CALL GPLOTP(GKX,GKY(1,2),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(0.0,0.0,0.0)

!     ----- Fig.2   R-BEAM-----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY( NSTP+1,1)=GUCLIP(RAYB(23,NSTP))
            GKY( NSTP+1,2)=GUCLIP(RAYB(24,NSTP))
         ENDDO
      ENDDO

      CALL MOVE(1.7,8.1)
      CALL SETRGB(0.0,0.0,0.0)
      CALL TEXT('R-beam',6)
      CALL GMNMX1(GKY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
      CALL GMNMX1(GKY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)
      GYMIN1=MIN(GYMINA,GYMINB)
      GYMAX1=MAX(GYMAXA,GYMAXB)     
      CALL GQSCAL (GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
!      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,0.0,0.10)
      CALL SETLIN(0,2,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.01,0.1,9)
!      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
!      CALL GVALUE(0.0,2*GXSTEP,GYMIN,GYSTEP,2)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)
      CALL GPLOTP(GKX,GKY(1,1),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(1.0,0.0,0.0)
      CALL GPLOTP(GKX,GKY(1,2),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(0.0,0.0,0.0)

!      CALL MOVE(9.0,16.4)
!      CALL TEXT ('X1= ',3)
!      CALL NUMBR(GKX(K,1),'(F5.3)',5)
!      CALL MOVE(5.0,16.4)
!      CALL TEXT ('dMIN1= ',6)
!      CALL NUMBR(GKY(K,1),'(F9.7)',9)

!      CALL MOVE(9.0,16.0)
!      CALL TEXT ('X2= ',3)
!      CALL NUMBR(GKX(K,1),'(F5.3)',5)
!      CALL MOVE(5.0,16.0)
!      CALL TEXT ('dMIN2= ',6)
!      CALL NUMBR(GKY(K,2),'(F9.7)',9)

      CALL MOVE(1.0,16.4)
      CALL TEXT ('SMAX= ',6)
      CALL NUMBD(smax,'(F5.2)',5)
      CALL MOVE(11.0,16.0)
      CALL TEXT ('RCURV= ',6)
      CALL NUMBD(wr_nray_status%RCURVA,'(F5.2)',5)
      CALL MOVE(11.0,16.4)
      CALL TEXT ('RBRAD= ',6)
      CALL NUMBD(wr_nray_status%RBRADA,'(F5.2)',5)

!     ----- Fig.3  Beam-Rot -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY(NSTP+1,1)=GUCLIP(RAYB(25,NSTP))
            GKY(NSTP+1,2)=GUCLIP(RAYB(26,NSTP))
         ENDDO
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLIN(0,2,7)
      CALL MOVE(13.5,16.1)
      CALL TEXT('Beam-rot',8) 
      CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,0.0,180.0)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.0,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,1)
      CALL GSCALE(0.0,0.0,0.0,30.0,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,30.0,0)
      CALL SETRGB(0.0,0.0,1.0)
      CALL GPLOTP(GKX,GKY(1,1),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(1.0,0.0,0.0)
      CALL GPLOTP(GKX,GKY(1,2),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(0.0,0.0,0.0)

!     -------------------------------------------------------------------
!     ----- CALCULATE RADIAL DEPOSITION PROFILE (without beam radial)----
!     -------------------------------------------------------------------
      DRHO=1.D0/nrsmax
      DO NRDIV=1,nrsmax
         GPX(NRDIV)=GUCLIP(NRDIV*DRHO)
      ENDDO
      DO NRAY=1,NRAYMAX
         DO NRDIV=1,nrsmax
            GPY(NRDIV,NRAY)=0.0
         ENDDO
      ENDDO

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)-1
            CALL pl_mag_old(RAYS(1,NSTP,NRAY),RAYS(2,NSTP,NRAY)**2, &
                            RAYS(3,NSTP,NRAY),RHON1)
            NRS1=INT(RHON1/DRHO)+1
            CALL pl_mag_old(RAYS(1,NSTP+1,NRAY),RAYS(2,NSTP+1,NRAY)**2, &
                            RAYS(3,NSTP+1,NRAY),RHON2)
            NRS2=INT(RHON2/DRHO)+1
            NDR=ABS(NRS2-NRS1)

            IF(NDR.EQ.0) THEN
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY)+GUCLIP(RAYS(8,NSTP+1,NRAY))
            ELSE IF(NRS1.LT.NRS2) THEN
               SDR=(RHON2-RHON1)/DRHO
               DELPWR=RAYS(8,NSTP+1,NRAY)/SDR
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY) &
                             +GUCLIP((REAL(NRS1)-RHON1/DRHO)*DELPWR)
               DO NR=NRS1+1,NRS2-1
                  GPY(NR,NRAY)=GPY(NR,NRAY)+GUCLIP(DELPWR)
               ENDDO
               GPY(NRS2,NRAY)=GPY(NRS2,NRAY) &
                             +GUCLIP((RHON2/DRHO-DBLE(NRS2-1))*DELPWR)
            ELSE
               SDR=(RHON1-RHON2)/DRHO
               DELPWR=RAYS(8,NSTP+1,NRAY)/SDR
               GPY(NRS2,NRAY)=GPY(NRS2,NRAY) &
                             +GUCLIP((REAL(NRS2)-RHON2/DRHO)*DELPWR)
               DO NR=NRS2+1,NRS1-1
                  GPY(NR,NRAY)=GPY(NR,NRAY)+GUCLIP(DELPWR)
               ENDDO
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY) &
                            +GUCLIP((RHON1/DRHO-DBLE(NRS1-1))*DELPWR)
            ENDIF
         ENDDO
      ENDDO

      DO NRAY=1,NRAYMAX
      DO NRDIV=1,nrsmax
         GPY(NRDIV,NRAY)=GPY(NRDIV,NRAY) &
                        /GUCLIP(2*PI*(DBLE(NRDIV)-0.5D0)*DRHO*DRHO)
      ENDDO
!         WRITE(6,'(5(I3,1PE12.4))') (NRDIV,GPY(NRDIV,NRAY),NRDIV=1,nrsmax)
      ENDDO

!     -------------------------------------------------------------------
!     ----- CALCULATE RADIAL DEPOSITION PROFILE (with beam radial)-------
!     -------------------------------------------------------------------

      IF (MODELG.EQ.3.OR.MODELG.EQ.5.OR.MODELG.EQ.8) THEN
         DRHO=1.D0/nrsmax  
         DO NRAY=1,NRAYMAX
            DO NRDIV=1,nrsmax
               GPX (NRDIV)=GUCLIP(NRDIV*DRHO)
               GPAY(NRDIV,NRAY)=0.0
            ENDDO
         ENDDO

         DO NRAY=1,NRAYMAX
            DO NSTP=0,NSTPMAX_NRAY(NRAY)-1
               DRAD=1.D0/DBLE(nrsmax)

!     ----- CALC SLOPE OF RAY TRAJECTORY-----------------------

               RLA =SQRT(RAYS(1,NSTP,NRAY)**2+RAYS(2,NSTP,NRAY)**2)
               RLAD=SQRT(RAYS(1,NSTP+1,NRAY)**2+RAYS(2,NSTP+1,NRAY)**2)
               ZLA =RAYS(3,NSTP ,NRAY)
               ZLAD=RAYS(3,NSTP+1,NRAY)
               DRL=RLAD-RLA
               DZL=ZLAD-ZLA
!        ---------------------------------------------------------

!        -----(i)ROTATE DRL and DZL (90 Degree)-------------------
               DRLS=-DZL
               DZLS= DRL
               DABSS=SQRT(DRLS**2+DZLS**2)
               DRLSN=DRLS/DABSS
               DZLSN=DZLS/DABSS
!        ---------------------------------------------------------

!        ----- JUDGE MAGNETIC SURFACE-----------------------------
               DO NRDIV=0,nrsmax
                  FRLRO1(NRDIV)=RLA+NRDIV*DRAD*DRLSN*RAYB(23,NSTP)
                  FZLRO1(NRDIV)=ZLA+NRDIV*DRAD*DZLSN*RAYB(23,NSTP)
               ENDDO

!            ------ BEAM RADIAL DIRECTION ----------------

               RLDMAX1=FRLRO1(nrsmax)
               ZLDMAX1=FZLRO1(nrsmax)
               CALL GETRZ(RLDMAX1,ZLDMAX1,PHIL,BR,BZ,BPHI,RHOMAX1)
               NRSMAX1=INT(RHOMAX1/DRHO)+1
               CALL GETRZ(RLA,ZLA,PHIL,BR,BZ,BPHI,RHOMIN1)
               NRSMIN1=INT(RHOMIN1/DRHO)+1
               NDBRD1=(NRSMAX1-NRSMIN1)
               NABSM1=ABS(NDBRD1)

               
           DO NRDIV=0,nrsmax-1
               RL1=FRLRO1(NRDIV)
               ZL1=FZLRO1(NRDIV)
               CALL GETRZ(RL1,ZL1,PHIL,BR,BZ,BPHI,RHON1)
!               WRITE(6,'(A,1P3E12.4)') 'RL1: ',RL1,ZL1,RHON1
               NRSA1=INT(RHON1/DRHO)+1

               RLD1=FRLRO1(NRDIV+1)
               ZLD1=FZLRO1(NRDIV+1)
               CALL GETRZ(RLD1,ZLD1,PHIL,BR,BZ,BPHI,RHOD1)
!               WRITE(6,'(A,1P3E12.4)') 'RLD1: ',RLD1,ZLD1,RHOD1
               NRSDA1=INT(RHOD1/DRHO)+1

                IF (NRSDA1.NE.NRSA1) THEN
                   RLMA1(NSTP,ABS(NRSDA1-NRSMIN1))=RL1
                   ZLMA1(NSTP,ABS(NRSDA1-NRSMIN1))=ZL1
                ENDIF
            ENDDO
!           
            IF (NDBRD1.EQ.0) THEN
               GPAY(NRSMIN1,NRAY)=GPAY(NRSMIN1,NRAY) &
                                 +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY))
            ELSE IF (NDBRD1.LT.0) THEN 
               RLMA1(NSTP,0)=RLA
               ZLMA1(NSTP,0)=ZLA
               DO NCS=1,NABSM1
               SIGMA=RAYB(23,NSTP)/3.D0
!               SIGMA=5.D0*RAYB(23,NSTP)
               FASSX1(NSTP,NCS)=SQRT((RLMA1(NSTP,NCS)-RLA)**2 &
                               +(ZLMA1(NSTP,NCS)-ZLA)**2)
!               FASSZ1(NSTP,NCS)=FASSX1(NSTP,NCS)/RAYB(23,NSTP)
               FASSZ1(NSTP,NCS)=FASSX1(NSTP,NCS)/SIGMA
               DELP1(NSTP,NCS)=ERF0(FASSZ1(NSTP,NCS)) &
                              -ERF0(FASSZ1(NSTP,NCS-1))
!               WRITE(6,*) DELP1(NSTP,NCS)
               GPAY(NRSMIN1-NCS+1,NRAY)=GPAY(NRSMIN1-NCS+1,NRAY) &
                    +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY)*DELP1(NSTP,NCS))
               ENDDO
               GPAY(NRSMAX1,NRAY)=GPAY(NRSMAX1,NRAY) &
                    +GUCLIP(0.5D0-0.5D0*ERF0(FASSZ1(NSTP,NABSM1))) &
                    *GUCLIP(RAYS(8,NSTP+1,NRAY))
!             
            ELSE 
               RLMA1(NSTP,0)=RLA
               ZLMA1(NSTP,0)=ZLA
               DO NCS=1,NABSM1
               SIGMA=RAYB(23,NSTP)/3.D0
               FASSX1(NSTP,NCS)=SQRT((RLMA1(NSTP,NCS)-RLA)**2 &
                               +(ZLMA1(NSTP,NCS)-ZLA)**2)
!               FASSZ1(NSTP,NCS)=FASSX1(NSTP,NCS)/RAYB(23,NSTP)
               FASSZ1(NSTP,NCS)=FASSX1(NSTP,NCS)/SIGMA
               DELP1(NSTP,NCS)=ERF0(FASSZ1(NSTP,NCS)) &
                              -ERF0(FASSZ1(NSTP,NCS-1))
!               WRITE(6,*) DELP1(NSTP,NCS)
               GPAY(NRSMIN1+NCS-1,NRAY)=GPAY(NRSMIN1+NCS-1,NRAY) &
                    +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY)*DELP1(NSTP,NCS))
               ENDDO
               GPAY(NRSMAX1,NRAY)=GPAY(NRSMAX1,NRAY) &
                    +0.5*GUCLIP(1.0D0-ERF0(FASSZ1(NSTP,NABSM1))) &
                    *GUCLIP(RAYS(8,NSTP+1,NRAY))
            ENDIF
            
!      ------------------------------------------------------

!      -----(ii)ROTATE DRLS and DZLS (-90 Degree)
            DRLT= DZL
            DZLT=-DRL
            DABST=SQRT(DRLT**2+DZLT**2)
            DRLTN=DRLT/DABST
            DZLTN=DZLT/DABST
!      ------------------------------------------------------

!      ------ JUDGE MAGNETIC SURFACE ------------------------
            DO NRDIV=0,nrsmax
               FRLRO2(NRDIV)=RLA+NRDIV*DRAD*DRLTN*RAYB(23,NSTP)
               FZLRO2(NRDIV)=ZLA+NRDIV*DRAD*DZLTN*RAYB(23,NSTP)
            ENDDO

               RLDMAX2=FRLRO2(nrsmax)
               ZLDMAX2=FZLRO2(nrsmax)
               CALL GETRZ(RLDMAX2,ZLDMAX2,PHIL,BR,BZ,BPHI,RHOMAX2)
!               WRITE(6,'(A,1P3E12.4)') 'RLD: ',RLDMAX2,ZLDMAX2,RHOMAX2
               NRSMAX2=INT(RHOMAX2/DRHO)+1
               CALL GETRZ(RLA,ZLA,PHIL,BR,BZ,BPHI,RHOMIN2)
!               WRITE(6,'(A,1P3E12.4)') 'RLA: ',RLA,ZLA,RHOMIN2
               NRSMIN2=INT(RHOMIN2/DRHO)+1
               NDBRD2=(NRSMAX2-NRSMIN2)
               NABSM2=ABS(NDBRD2)

            DO NRDIV=0,nrsmax-1
               RL2=FRLRO2(NRDIV)
               ZL2=FZLRO2(NRDIV)
               CALL GETRZ(RL2,ZL2,PHIL,BR,BZ,BPHI,RHON2)
!              WRITE(6,'(A,1P3E12.4)') 'RL2: ',RL2,ZL2,RHON2
               NRSA2=INT(RHON2/DRHO)+1

               RLD2=FRLRO2(NRDIV+1)
               ZLD2=FZLRO2(NRDIV+1)
               CALL GETRZ(RLD2,ZLD2,PHIL,BR,BZ,BPHI,RHOD2)
!               WRITE(6,'(A,1P3E12.4)') 'RL2D:',RLD2,ZLD2,RHOD2
               NRSDA2=INT(RHOD2/DRHO)+1
!            
               IF (NRSDA2.NE.NRSA2) THEN
                     RLMA2(NSTP,ABS(NRSDA2-NRSMIN2))=RL2
                     ZLMA2(NSTP,ABS(NRSDA2-NRSMIN2))=ZL2
               ENDIF
            ENDDO
! 
            IF (NDBRD2.EQ.0) THEN
               GPAY(NRSMIN2,NRAY)=GPAY(NRSMIN2,NRAY) &
                                 +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY))
            ELSE IF (NDBRD2.LT.0) THEN
               RLMA2(NSTP,0)=RLA
               ZLMA2(NSTP,0)=ZLA
               DO NCS2=1,NABSM2
               SIGMA=RAYB(23,NSTP)/3.D0
!               SIGMA=5.D0*RAYB(23,NSTP)
               FASSX2(NSTP,NCS2)=SQRT((RLMA2(NSTP,NCS2)-RLA)**2 &
                                +(ZLMA2(NSTP,NCS2)-ZLA)**2)
!               FASSZ2(NSTP,NCS2)=FASSX2(NSTP,NCS2)/RAYB(23,NSTP)
               FASSZ2(NSTP,NCS2)=FASSX2(NSTP,NCS2)/SIGMA
               DELP2(NSTP,NCS2)=ERF0(FASSZ2(NSTP,NCS2)) &
                               -ERF0(FASSZ2(NSTP,NCS2-1))
               GPAY(NRSMIN2-NCS2+1,NRAY)=GPAY(NRSMIN2-NCS2+1,NRAY) &
                    +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY)*DELP2(NSTP,NCS2))
               ENDDO
               GPAY(NRSMAX2,NRAY)=GPAY(NRSMAX2,NRAY) &
                    +(0.5-0.5*GUCLIP(ERF0(FASSZ2(NSTP,NABSM2)))) &
                    *GUCLIP(RAYS(8,NSTP+1,NRAY))
            ELSE
               RLMA2(NSTP,0)=RLA
               ZLMA2(NSTP,0)=ZLA
               DO NCS2=1,NABSM2
               SIGMA=RAYB(23,NSTP)/3.D0
!                SIGMA=5.D0*RAYB(23,NSTP)
               FASSX2(NSTP,NCS2)=SQRT((RLMA2(NSTP,NCS2)-RLA)**2 &
                                +(ZLMA2(NSTP,NCS2)-ZLA)**2)
!               FASSZ2(NSTP,NCS2)=FASSX2(NSTP,NCS2)/RAYB(23,NSTP)
                FASSZ2(NSTP,NCS2)=FASSX2(NSTP,NCS2)/SIGMA
               DELP2(NSTP,NCS2)=ERF0(FASSZ2(NSTP,NCS2)) &
                               -ERF0(FASSZ2(NSTP,NCS2-1))
               GPAY(NRSMIN2+NCS2-1,NRAY)=GPAY(NRSMIN2+NCS2-1,NRAY) &
                    +0.5*GUCLIP(RAYS(8,NSTP+1,NRAY)*DELP2(NSTP,ncs2))
               ENDDO
               GPAY(NRSMAX2,NRAY)=GPAY(NRSMAX2,NRAY) &
                    +(0.5-0.5*GUCLIP(ERF0(FASSZ2(NSTP,NABSM2)))) &
                    *GUCLIP(RAYS(8,NSTP+1,NRAY))
            ENDIF
         ENDDO
      ENDDO


      DO NRAY=1,NRAYMAX
         DO NRDIV=1,nrsmax
            GPAY(NRDIV,NRAY)=GPAY(NRDIV,NRAY) &
                            /GUCLIP(2*PI*(DBLE(NRDIV)-0.5D0)*DRHO*DRHO)
         ENDDO
!         WRITE(6,'(5(I3,1PE12.4))') (NRDIV,GPAY(NRDIV,NRAY),NRDIV=1,nrsmax)
      ENDDO
      ENDIF
!     -----Fig.4 draw deposition profile -----

      IF(MOD(MDLWRG/2,2).EQ.0) THEN
         CALL GQSCAL(0.0,1.0,GXMIN,GXMAX,GXSTEP)
      ELSE
         CALL GQSCAL(GUCLIP(RHOGMN),GUCLIP(RHOGMX),GXMIN,GXMAX,GXSTEP)
      ENDIF
!
      CALL GMNMX2(GPY,NSTPMAX+1,1,nrsmax,1,1,NRAYMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMIN,GYSMAX1,GYSCAL1)
      GYSMIN=0.0
      GYSMAX2=0.1 
      GYSCAL2=2.E-2
      GYSMAX=MAX(GYSMAX1,GYSMAX2)
      GYSCAL=MAX(GYSCAL1,GYSCAL2)

      CALL MOVE(13.5,8.1)
      CALL TEXT('ABS POWER',10)
!      CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,GYSMIN,GYSMAX)
      CALL GDEFIN(13.5,23.5,1.0,8.0,GXMIN,GXMAX,GYSMIN,GYSMAX)
      CALL SETRGB(0.0,0.0,0.0)
      CALL GFRAME
      CALL GSCALE(GXMIN,GXSTEP,GYSMIN,GYSCAL,0.1,9)      
      CALL GVALUE(GXMIN,2*GXSTEP,GYSMIN,GYSCAL,2)
!      CALL SETLIN(0,0,4)
      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,2,7-MOD(NRAY-1,5))
         CALL GPLOTP(GPX,GPY(1,NRAY),1,nrsmax,1,0,0,0)
         CALL SETRGB(0.0,0.0,0.0)
      ENDDO

      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,2,7-MOD(NRAY,5))
         CALL GPLOTP(GPX,GPAY(1,NRAY),1,nrsmax,1,0,0,0)
         CALL SETRGB(1.0,0.0,0.0)
      ENDDO     
      DEALLOCATE(RLMA1)
      DEALLOCATE(RLMA2)
      DEALLOCATE(ZLMA1)
      DEALLOCATE(ZLMA2)
      DEALLOCATE(FASSX1)
      DEALLOCATE(FASSX2)
      DEALLOCATE(FASSZ1)
      DEALLOCATE(FASSZ2)
      DEALLOCATE(DELP1)
      DEALLOCATE(DELP2)

      CALL WRGPRM(1)
      CALL PAGEE
      RETURN
    END SUBROUTINE WRGRFB1

    !     ***** Rays ***********
    
      SUBROUTINE WRGRFB2

      USE wrcomm
      USE plprof,ONLY:PL_MAG_OLD
      IMPLICIT NONE
      REAL:: GSX (NSTPMAX+1,NRAYM),GRX (NSTPMAX+1,NRAYM)
      REAL:: GKY1(NSTPMAX+1,NRAYM),GKY2(NSTPMAX+1,NRAYM)
      INTEGER:: NRAY,NSTP
      REAL:: GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP
      REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP,GYMIN2,GYMAX2
      REAL(rkind):: XL,YL,ZL,RHON
      EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
      EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
      EXTERNAL MOVE,TEXT,SETLNW
      REAL:: GUCLIP

      CALL PAGES
      CALL SETFNT(32)
      CALL SETCHS(0.25,0.)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!    ------- X Axis (direction along ray reference)-------

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GSX(NSTP+1,NRAY)=GUCLIP(RAYS(0,NSTP,NRAY))
         ENDDO
         CALL GMNMX1(GSX(1,1),1,NSTPMAX_NRAY(1),1,GXMIN1,GXMAX1)
         CALL GMNMX1(GSX(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GXMIN,GXMAX)
         GXMAX1=MAX(GXMAX1,GXMAX)
      ENDDO

      CALL GQSCAL(0.0,GXMAX1,GXMIN,GXMAX,GXSTEP)

!      ------ NRAY BEAM WIDTH 1 -------
      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY1(NSTP+1,NRAY)=GUCLIP(RAYRB1(NSTP,NRAY))
         ENDDO
         CALL GMNMX1(GKY1(1,1),1,NSTPMAX_NRAY(1),1,GYMIN1,GYMAX1)
         CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMAX1=MAX(GYMAX1,GYMAX)
         GYMIN1=MIN(GYMIN1,GYMIN)
      ENDDO

      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(1.7,16.1)
      CALL TEXT('S-Beam1',7)
      CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,0.0,GYMAX)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,2*GYSTEP,2)


      DO NRAY=1,NRAYM
      CALL SETLIN(0,2,7-MOD(NRAY-1,5))
      CALL GPLOTP(GSX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO

!            ------ NRAY BEAM WIDTH 2 -------
      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY2(NSTP+1,NRAY)=GUCLIP(RAYRB2(NSTP,NRAY))
         ENDDO
         CALL GMNMX1(GKY2(1,1),1,NSTPMAX_NRAY(1),1,GYMIN2,GYMAX2)
         CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMAX2=MAX(GYMAX2,GYMAX)
         GYMIN2=MIN(GYMIN2,GYMIN)
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
      CALL GQSCAL(GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.5,16.1)
      CALL TEXT('S-Beam2',7)
      CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,0.0,GYMAX)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,2*GYSTEP,2)

      DO NRAY=1,NRAYM
      CALL SETLIN(0,2,7-MOD(NRAY-1,5))
      CALL GPLOTP(GSX(1,NRAY),GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
!   
!    ------- X Axis (radial direction)-------

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
            GRX(NSTP+1,NRAY)=GUCLIP(RHON)
         ENDDO
         CALL GMNMX1(GRX(1,1),1,NSTPMAX_NRAY(1),1,GXMIN1,GXMAX1)
         CALL GMNMX1(GRX(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GXMIN,GXMAX)
         GXMAX1=MAX(GXMAX1,GXMAX)
      ENDDO

      CALL GQSCAL(0.0,GXMAX1,GXMIN,GXMAX,GXSTEP)
      

!     ------ NRAY BEAM WIDTH1 ----------------

      CALL SETRGB(0.0,0.0,0.0)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(1.7,8.1)
      CALL TEXT('R-Beam1',7)
      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,0.0,GYMAX)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,2*GYSTEP,2)

      DO NRAY=1,NRAYM
      CALL SETLIN(0,2,7-MOD(NRAY-1,5))
      CALL GPLOTP(GRX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO

!     ------ NRAY BEAM WIDTH2 ----------------

      CALL SETRGB(0.0,0.0,0.0)
      CALL GQSCAL(GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.5,8.1)
      CALL TEXT('R-Beam2',7)
      CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,0.0,GYMAX)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,2*GYSTEP,2)

      DO NRAY=1,NRAYM
      CALL SETLIN(0,2,7-MOD(NRAY-1,5))
      CALL GPLOTP(GRX(1,NRAY),GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO

      CALL PAGEE
      RETURN
    END SUBROUTINE WRGRFB2

!     ***** ANGLE AND WS(I,J) *****

      SUBROUTINE WRGRFB3

      USE wrcomm
      USE plprof,ONLY:PL_MAG_OLD
      IMPLICIT NONE
      REAL:: GKX(NSTPMAX+1,1),GKY(NSTPMAX+1,2)
      REAL:: GKSY(NSTPMAX+1,6),GKPY(NSTPMAX+1,6)
      REAL:: GTHY(NSTPMAX+1,1)
      INTEGER:: NRAY,NSTP,N
      REAL(rkind):: Xl,YL,ZL,RHON
      REAL:: GXMIN,GXMAX,GXSTEP
      REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP
      EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
      EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
      EXTERNAL MOVE,TEXT,SETLNW
      REAL:: GUCLIP

      CALL PAGES
      CALL SETFNT(32)
      CALL SETCHS(0.25,0.)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!     ----- X AXIS -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
            GKX( NSTP+1,NRAY)=GUCLIP(RHON)
         ENDDO
      ENDDO

!      CALL GMNMX1(GKX,1,NSTPMAX_NRAY(1),1,GXMIN1,GXMAX1)
!      CALL GQSCAL(GKMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP) 
!      CALL GQSCAL(3.0,4.0,GXMIN,GXMAX,GXSTEP)
       CALL GQSCAL(0.0,4.0,GXMIN,GXMAX,GXSTEP)

!     ----- Fig.1  WS-----

      DO NRAY=1,NRAYMAX         
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            DO N=1,6
            GKSY(NSTP+1,N)=GUCLIP(RAYB(25+N,NSTP))
            ENDDO
         ENDDO
      ENDDO

      CALL GMNMX1(GKSY(1,1),1,NSTPMAX_NRAY(1),1,GYMIN1,GYMAX1)

      DO N=2,6
       CALL GMNMX1(GKSY(1,N),1,NSTPMAX_NRAY(1),1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)
      CALL GQSCAL (GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(1.7,16.1)
      CALL TEXT('X-WS',6)
!      CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,-40.0,40.0)
      CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,GYMIN,GYMAX)
      CALL GFRAME
!      CALL GSCALE(0.0,GXSTEP,0.0,5.0,0.1,9)
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,GYMIN,GYSTEP,1)
      CALL SETRGB(0.0,0.0,1.0)

      DO N=1,6
      CALL GPLOTP(GKX,GKSY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETLIN(0,0,N)
      ENDDO
      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!    
!     
!     
!      


!     ------Fig2  WP-------
!   
!    
      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            DO N=1,6
            GKPY( NSTP+1,N)=GUCLIP(RAYB(31+N,NSTP))
            ENDDO
         ENDDO
      ENDDO

      CALL GMNMX1(GKPY(1,1),1,NSTPMAX_NRAY(1),1,GYMIN1,GYMAX1)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

      DO N=2,6
       CALL GMNMX1(GKPY(1,N),1,NSTPMAX_NRAY(1),1,GYMIN,GYMAX)
       GYMIN1=MIN(GYMIN1,GYMIN)
       GYMAX1=MAX(GYMAX1,GYMAX)
      ENDDO

      CALL MOVE(1.7,8.1)
      CALL TEXT('X-WP',6)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
!       CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,GYMIN,GYSTEP,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)
!      
      DO N=1,6
      CALL GPLOTP(GKX,GKPY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETLIN(0,0,N)
      ENDDO

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!     ----- Fig.3   Angle-----

      DO NRAY=1,NRAYMAX         
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GTHY( NSTP+1,1)=GUCLIP(RAYB(25,NSTP))
         ENDDO
      ENDDO

!      CALL GMNMX1(GTHY(1,1),1,NSTPMAX_NRAY(1),1,GYMIN1,GYMAX1)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)
!      CALL GQSCAL (GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.5,16.1)
      CALL TEXT('Angle',6)
      CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,-40.0,40.0)
!      CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,GYMIN1,GYMAX1)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,5.0,0.1,9)
!      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
!      CALL GVALUE(0.0,2*GXSTEP,GYMIN,GYSTEP,1)
       CALL GVALUE(0.0,2*GXSTEP,0.0,10.0,1)
      CALL SETRGB(0.0,0.0,1.0)

      CALL GPLOTP(GKX,GTHY(1,1),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!    
!     
!     
!      



!     -------fig4  LAMDA------------------------

            DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
!            GKY( NSTP+1,1)=GUCLIP(RAYB(47,NSTP))
!            GKY( NSTP+1,2)=GUCLIP(RAYB(48,NSTP))
            GKY( NSTP+1,1)=GUCLIP(RAYB(20,NSTP))
         ENDDO
      ENDDO

      CALL GMNMX1(GKY(1,1),1,NSTPMAX_NRAY(1),1,GYMIN1,GYMAX1)
!      CALL GMNMX1(GKY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
!      CALL GMNMX1(GKY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)

!      GYMIN1=MIN(GYMINA,GYMINB)
!      GYMAX1=MAX(GYMAXA,GYMAXB)
      CALL MOVE(13.5,8.1)
!      CALL TEXT('LAMDA',6)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
!      CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
!      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
      CALL GDEFIN(13.5,23.5,1.0,8.0,GXMIN,GXMAX,GYMIN,GYMAX)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(GXMIN,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(GXMIN,2*GXSTEP,GYMIN,GYSTEP,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)
!      
      CALL GPLOTP(GKX,GKY(1,1),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL GPLOTP(GKX,GKY(1,2),1,NSTPMAX_NRAY(1)+1,1,0,0,0)
      CALL SETLIN(0,0,7)

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)
      CALL WRGPRM(1)
      CALL PAGEE
      RETURN
    END SUBROUTINE WRGRFB3

!     ***** K,(K*B),K*(K*B),VG *****

      SUBROUTINE WRGRFB4

      USE wrcomm
      IMPLICIT NONE
      REAL:: GKX(NSTPMAX+1,1)
      REAL:: GKY(NSTPMAX+1,3),GKBY(NSTPMAX+1,3),GKKBY(NSTPMAX+1,3)
      REAL:: GVGY(NSTPMAX+1,3)
      INTEGER:: NRAY,NSTP,N
      REAL:: GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP
!      REAL:: GYMNINA,GYMAXA,GYMINB,GYMAXB,GYMINC,GYMAXC
!      REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP
      EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
      EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
      EXTERNAL MOVE,TEXT,SETLNW
      REAL:: GUCLIP

      CALL PAGES
      CALL SETFNT(32)
      CALL SETCHS(0.25,0.)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!     ----- X AXIS -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKX( NSTP+1,NRAY)=GUCLIP(RAYB(0,NSTP))
         ENDDO
      ENDDO

      CALL GMNMX1(GKX,1,NSTPMAX_NRAY(1),1,GXMIN1,GXMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP) 

!     ------fig 1  K----------

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKY( NSTP+1,1)=GUCLIP(RAYB(38,NSTP))
            GKY( NSTP+1,2)=GUCLIP(RAYB(39,NSTP))
            GKY( NSTP+1,3)=GUCLIP(RAYB(40,NSTP))
         ENDDO
      ENDDO
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!      CALL GMNMX1(GKY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
!      CALL GMNMX1(GKY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)
!      CALL GMNMX1(GKY(1,3),1,NSTPMAX_NRAY(1),1,GYMINC,GYMAXC)

!      GYMIN1=MIN(GYMINA,GYMINB,GYMINC)
!      GYMAX1=MAX(GYMAXA,GYMAXB,GYMAXC)
      CALL MOVE(1.7,16.1)
      CALL TEXT('K',6)
!      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL GDEFIN(1.7,11.7,9.0,16.0,0.0,GXMAX,-1.2,1.2)
!       CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.1,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.2,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)

      DO N=1,3
        CALL GPLOTP(GKX,GKY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)      
        CALL SETLIN(0,0,N+2)
      ENDDO  

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!     ------fig 2 K*B---------

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKBY( NSTP+1,1)=GUCLIP(RAYB(41,NSTP))
            GKBY( NSTP+1,2)=GUCLIP(RAYB(42,NSTP))
            GKBY( NSTP+1,3)=GUCLIP(RAYB(43,NSTP))
         ENDDO
      ENDDO
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!      CALL GMNMX1(GKBY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
!      CALL GMNMX1(GKBY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)
!      CALL GMNMX1(GKBY(1,3),1,NSTPMAX_NRAY(1),1,GYMINC,GYMAXC)

!      GYMIN1=MIN(GYMINA,GYMINB,GYMINC)
!      GYMAX1=MAX(GYMAXA,GYMAXB,GYMAXC)
      CALL MOVE(1.7,8.1)
      CALL TEXT('K*B',6)
!      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-1.2,1.2)
!       CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.1,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.2,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)

      DO N=1,3
        CALL GPLOTP(GKX,GKBY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)      
        CALL SETLIN(0,0,N+2)
      ENDDO  

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!     ------fig 3 K*(K*B)-----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GKKBY( NSTP+1,1)=GUCLIP(RAYB(44,NSTP))
            GKKBY( NSTP+1,2)=GUCLIP(RAYB(45,NSTP))
            GKKBY( NSTP+1,3)=GUCLIP(RAYB(46,NSTP))
         ENDDO
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)

!      CALL GMNMX1(GKKBY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
!      CALL GMNMX1(GKKBY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)
!      CALL GMNMX1(GKKBY(1,3),1,NSTPMAX_NRAY(1),1,GYMINC,GYMAXC)

!      GYMIN1=MIN(GYMINA,GYMINB,GYMINC)
!      GYMAX1=MAX(GYMAXA,GYMAXB,GYMAXC)
      CALL MOVE(13.5,16.1)
      CALL TEXT('K*(K*B)',6)
!      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL GDEFIN(13.5,23.5,9.0,16.0,0.0,GXMAX,-1.2,1.2)
!       CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.1,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.2,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)

      DO N=1,3
        CALL GPLOTP(GKX,GKKBY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)      
        CALL SETLIN(0,0,N+2)
      ENDDO  

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)

!     ------fig 4 VG ---------

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GVGY( NSTP+1,1)=GUCLIP(RAYB(49,NSTP))
            GVGY( NSTP+1,2)=GUCLIP(RAYB(50,NSTP))
            GVGY( NSTP+1,3)=GUCLIP(RAYB(51,NSTP))
         ENDDO
      ENDDO

      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLNW(0.035)
!      CALL GMNMX1(GVGY(1,1),1,NSTPMAX_NRAY(1),1,GYMINA,GYMAXA)
!      CALL GMNMX1(GVGY(1,2),1,NSTPMAX_NRAY(1),1,GYMINB,GYMAXB)
!      CALL GMNMX1(GVGY(1,3),1,NSTPMAX_NRAY(1),1,GYMINC,GYMAXC)

!      GYMIN1=MIN(GYMINA,GYMINB,GYMINC)
!      GYMAX1=MAX(GYMAXA,GYMAXB,GYMAXC)
      CALL MOVE(13.5,8.1)
      CALL TEXT('VG',6)
!      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,-1.2,1.2)
!       CALL GDEFIN(1.7,11.7,1.0,8.0,0.0,GXMAX,-50.0,50.0)
!      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(0.0,GXSTEP,0.0,0.1,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.2,1)
!      CALL GVALUE(0.0,0.0,0.0,0.01,2)
      CALL SETRGB(0.0,0.0,1.0)

      DO N=1,3
        CALL GPLOTP(GKX,GVGY(1,N),1,NSTPMAX_NRAY(1)+1,1,0,0,0)      
        CALL SETLIN(0,0,N+2)
      ENDDO  

      CALL SETRGB(1.0,0.0,0.0)
      CALL SETRGB(0.0,0.0,0.0)
      CALL WRGPRM(1)
      CALL PAGEE 
      RETURN
    END SUBROUTINE WRGRFB4

!     ***** DRAW PARAMETERS *****

    SUBROUTINE WRGPRM(nray)

      USE wrcomm
      IMPLICIT NONE
      INTEGER nray
      EXTERNAL SETLIN

!     ----- DRAW PARAMETERS -----

      CALL SETLIN(0,0,7)
      CALL WRPRMT(1.0,18.0, 'BB   =',6, BB  ,'(F8.3)',8)
      CALL WRPRMT(1.0,17.6, 'RR   =',6, RR  ,'(F8.3)',8)
      CALL WRPRMT(1.0,17.2, 'RA   =',6, RA  ,'(F8.3)',8)
      CALL WRPRMT(1.0,16.8, 'RKAP =',6, RKAP,'(F8.3)',8)
      CALL WRPRMT(5.0,18.0, 'Q0   =',6, MIN(Q0,1000.D0),'(F8.3)',8)
      CALL WRPRMT(5.0,17.6, 'QA   =',6, MIN(QA,1000.D0),'(F8.3)',8)
      CALL WRPRMT(5.0,17.2, 'RF   =',6, RFIN(nray),'(1PE10.3)',10)
      CALL WRPRMT(5.0,16.8, 'RNPH =',6, RNPHIN(nray),'(F8.4)',8)
      CALL WRPRMT(9.0,18.0, 'PA(1)  =',8,PA(1),'(F8.4)',8)
      CALL WRPRMT(9.0,17.6, 'PZ(1)  =',8,PZ(1),'(F8.4)',8)
      CALL WRPRMT(9.0,17.2, 'PN(1)  =',8,PN(1),'(F8.4)',8)
      CALL WRPRMT(9.0,16.8, 'PNS(1) =',8,PNS(1),'(F8.4)',8)
      CALL WRPRMT(13.0,18.0,'PTPR(1)=',8,PTPR(1),'(F8.4)',8)
      CALL WRPRMT(13.0,17.6,'PTPP(1)=',8,PTPP(1),'(F8.4)',8)
      CALL WRPRMT(13.0,17.2,'PTS(1) =',8,PTS(1),'(F8.4)',8)
      CALL WRPRMT(13.0,16.8,'MODELP1=',8,DBLE(MODELP(1)),'(F4.0)',4)
      CALL WRPRMT(17.0,18.0,'PA(2)  =',8,PA(2),'(F8.4)',8)
      CALL WRPRMT(17.0,17.6,'PZ(2)  =',8,PZ(2),'(F8.4)',8)
      CALL WRPRMT(17.0,17.2,'PN(2)  =',8,PN(2),'(F8.4)',8)
      CALL WRPRMT(17.0,16.8,'PNS(2)= ',8,PNS(2),'(F8.4)',8)
      CALL WRPRMT(21.0,18.0,'PTPR(2)=',8,PTPR(2),'(F8.4)',8)
      CALL WRPRMT(21.0,17.6,'PTPP(2)=',8,PTPP(2),'(F8.4)',8)
      CALL WRPRMT(21.0,17.2,'PTS(2) =',8,PTS(2),'(F8.4)',8)
      CALL WRPRMT(21.0,16.8,'MODELP2=',8,DBLE(MODELP(2)),'(F4.0)',4)
      RETURN
    END SUBROUTINE WRGPRM

!***********************************************************************

    SUBROUTINE WRPRMT(X,Y,KTEXT,ITEXT,D,KFORM,IFORM)

      USE bpsd_kinds,ONLY: rkind
      REAL,INTENT(IN):: X,Y
      INTEGER,INTENT(IN):: ITEXT,IFORM
      CHARACTER(LEN=*),INTENT(IN):: KTEXT,KFORM
      REAL(rkind),INTENT(IN):: D
      EXTERNAL MOVE,TEXT,NUMBD

      CALL MOVE(X,Y)
      CALL TEXT(KTEXT,ITEXT)
      CALL NUMBD(D,KFORM,IFORM)
      RETURN
    END SUBROUTINE WRPRMT


!     ***** POLOIDAL TRAJECTORY AND POWER *****

    SUBROUTINE WRGRF11A

      USE wrcomm
      USE plprof,ONLY: PL_RZSU,PL_MAG_OLD
      USE wrsub,ONLY: wrcalk
      IMPLICIT NONE
      INTEGER,PARAMETER:: NSUM=501
      INTEGER,PARAMETER:: NGXL=101
      INTEGER,PARAMETER:: NGYL=101
      REAL:: GX(NSTPMAX+1),GY(NSTPMAX+1)
      REAL:: GPX(NSTPMAX+1),GPY(NSTPMAX+1,NRAYM)
      REAL:: GPRY(NSTPMAX+1,NRAYM)
      REAL:: GKX(NSTPMAX+1,NRAYM)
      REAL:: GKY1(NSTPMAX+1,NRAYM),GKY2(NSTPMAX+1,NRAYM)
      REAL:: GRSMIN,GRSMAX,GZSMIN,GZSMAX
      REAL:: GRMIN,GRMAX,GRSTEP,GZMIN,GZMAX,GZSTEP
      REAL:: GRORG,GZORG,GYMAXT,GYORG
      REAL:: GXMIN,GXMAX,GXSTEP
      REAL:: GYMIN1,GYMAX1,GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP
      INTEGER:: NRAY,NSTP,NRDIV,NRS1,NRS2,NDR,NR
      REAL(rkind):: DRHO,RHO0,XL1,XL2,SDR,DELPWR,XL,YL,ZL,RHON,RKPARA,RKPERP
      EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
      EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
      EXTERNAL MOVE,TEXT
      REAL:: GUCLIP
      INTEGER:: NGULEN

      GRSMIN=GUCLIP(r_corner(1))
      GRSMAX=GUCLIP(r_corner(2))
      GZSMIN=GUCLIP(z_corner(1))
      GZSMAX=GUCLIP(z_corner(3))
      CALL GQSCAL(GRSMIN,GRSMAX,GRMIN,GRMAX,GRSTEP)
      CALL GQSCAL(GZSMIN,GZSMAX,GZMIN,GZMAX,GZSTEP)
      IF(GRMIN*GRMAX.LT.0.D0) THEN
         GRORG=0.0
      ELSE
         GRORG=(INT(GRMIN/(2*GRSTEP))+1)*2*GRSTEP
      END IF
      IF(GZMIN*GZMAX.LT.0.D0) THEN
         GZORG=0.0
      ELSE
         GZORG=(INT(GZMIN/(2*GZSTEP))+1)*2*GZSTEP
      END IF


      CALL PAGES
      CALL SETCHS(0.25,0.)
      CALL SETFNT(32)
      CALL SETLIN(0,2,7)
      CALL GDEFIN(1.5,11.5,1.0,16.0,GRMIN,GRMAX,GZMIN,GZMAX)
      CALL GFRAME
      CALL GSCALE(GRORG,  GRSTEP,0.0,  GZSTEP,0.1,9)
      CALL GVALUE(GRORG,2*GRSTEP,0.0,0.0,NGULEN(2*GRSTEP))
      CALL GVALUE(0.0,0.0,0.0,2*GZSTEP,NGULEN(2*GZSTEP))

!     ----- RAY TRAJECTORY -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GX(NSTP+1)=GUCLIP(RAYS(1,NSTP,NRAY))
            GY(NSTP+1)=GUCLIP(RAYS(3,NSTP,NRAY))
!            WRNSTPE(6,'(A,I5,1P2E12.4)') 
!     &           'NNSTP,GX,GY=',NSTP+1,GX(NSTP+1),GY(NSTP+1)
         ENDDO
         CALL SETLIN(0,2,7-MOD(NRAY-1,5))
         CALL GPLOTP(GX,GY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

!     ----- CALCULATE RADIAL DEPOSITION PROFILE -----

      DRHO=(r_corner(2)-r_corner(1))/nrsmax
      RHO0=r_corner(1)
      DO NRDIV=1,nrsmax
         GPX(NRDIV)=GUCLIP(RHO0+(NRDIV-0.5D0)*DRHO)
      ENDDO
      DO NRAY=1,NRAYMAX
         DO NRDIV=1,nrsmax
             GPY(NRDIV,NRAY)=0.0
             GPRY(NRDIV,NRAY)=0.0
         ENDDO
      ENDDO

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)-1
            XL1=RAYS(1,NSTP,NRAY)
            NRS1=INT((XL1-RHO0)/DRHO)+1
            XL2=RAYS(1,NSTP+1,NRAY)
            NRS2=INT((XL2-RHO0)/DRHO)+1
            NDR=ABS(NRS2-NRS1)
            IF(MIN(NRS1,NRS2).GE.1.AND.MAX(NRS1,NRS2).LE.nrsmax) THEN
            IF(NDR.EQ.0) THEN
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY)+GUCLIP(RAYS(8,NSTP+1,NRAY))
            ELSE IF(NRS1.LT.NRS2) THEN
               SDR=(XL2-XL1)/DRHO
               DELPWR=RAYS(8,NSTP+1,NRAY)/SDR
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY) &
                    +GUCLIP((DBLE(NRS1)-XL1/DRHO)*DELPWR)
               DO NR=NRS1+1,NRS2-1
                  GPY(NR,NRAY)=GPY(NR,NRAY)+GUCLIP(DELPWR)
               ENDDO
               GPY(NRS2,NRAY)=GPY(NRS2,NRAY) &
                    +GUCLIP((XL2/DRHO-DBLE(NRS2-1))*DELPWR)
            ELSE
               SDR=(XL1-XL2)/DRHO
               DELPWR=RAYS(8,NSTP+1,NRAY)/SDR
               GPY(NRS2,NRAY)=GPY(NRS2,NRAY) &
                    +GUCLIP((DBLE(NRS2)-XL2/DRHO)*DELPWR)
               DO NR=NRS2+1,NRS1-1
                  GPY(NR,NRAY)=GPY(NR,NRAY)+GUCLIP(DELPWR)
               ENDDO
               GPY(NRS1,NRAY)=GPY(NRS1,NRAY) &
                    +GUCLIP((XL1/DRHO-DBLE(NRS1-1))*DELPWR)
            ENDIF
            ENDIF
!            WRITE(6,'(3I5,1P5E12.4))') &
!                NSTP,NRS1,NRS2,RHON1,RHON2,SDR,DELPWR,RAYS(8,NSTP+1,NRAY)
         ENDDO
!         WRITE(6,'(5(I3,1PE12.4))') (NRDIV,GPY(NRDIV,NRAY),NRDIV=1,nrsmax)
      ENDDO

!     ----- draw deposition profile -----

      CALL GQSCAL(GRSMIN,GRSMAX,GXMIN,GXMAX,GXSTEP)

      GYMAXT=0.0
      DO NRAY=1,NRAYMAX
         CALL GMNMX1(GPY(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN1,GYMAX1)
         GYMAXT=MAX(GYMAXT,GYMAX1)
      END DO
      CALL GQSCAL(0.0,GYMAXT,GYMIN,GYMAX,GYSTEP)
      GYMIN=0.0

      CALL MOVE(13.5,8.1)
      CALL TEXT('ABS POWER',10)
      IF(MOD(MDLWRG/2,2).EQ.0) THEN
         CALL GDEFIN(13.5,23.5,1.0,8.0,GXMIN,GXMAX,GYMIN,GYMAX)
         CALL SETLIN(0,2,7)
         CALL GFRAME
         CALL GSCALE(0.0,  GXSTEP,0.0,  GYSTEP,0.1,9)
         CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GYSTEP))
         CALL GVALUE(0.0,0.0,0.0,GYSTEP,NGULEN(GYSTEP))
      ENDIF
!      CALL SETLIN(0,0,4)
      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,2,7-MOD(NRAY-1,5))
         CALL GPLOTP(GPX,GPY(1,NRAY),1,nrsmax,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

!     ----- DRAW POWER FLUX -----
!      
!      GZSMIN=0.0
!      GZSMAX=1.1
!      GZSCAL=0.2

!      CALL MOVE(13.5,8.1)
!      CALL TEXT('POWER FLUX',10)
!      IF(MOD(MDLWRG/2,2).EQ.0) THEN
!         CALL GDEFIN(13.5,23.5,1.0,8.0,0.0,GXMAX,0.0,GZSMAX)
!         CALL SETLIN(0,2,7)
!         CALL GFRAME
!         CALL GSCALE(0.0,  GXSTEP,0.0,  GZSCAL,0.1,9)
!         CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
!         CALL GVALUE(0.0,0.0,0.0,GZSCAL,1)
!      ELSE
!         CALL GDEFIN(13.5,23.5,1.0,8.0,GXMIN,GXMAX,0.0,GXSMAX)
!         CALL SETLIN(0,2,7)
!         CALL GFRAME
!         GXORG=(INT(GXMIN/(2*GXSTEP))+1)*2*GXSTEP
!         CALL GSCALE(GXORG,  GXSTEP,0.0,  GZSCAL,0.1,9)
!         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
!         CALL GVALUE(0.0,0.0,0.0,GZSCAL,1)
!      ENDIF
!      DO NRAY=1,NRAYMAX
!         DO NSTP=0,NSTPMAX_NRAY(NRAY)
!            XL=RAYS(1,NSTP,NRAY)
!            YL=RAYS(2,NSTP,NRAY)
!            ZL=RAYS(3,NSTP,NRAY)
!            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
!	    GUX(NSTP+1)=GUCLIP(RHON)
!            GUY(NSTP+1)=GUCLIP(RAYS(7,NSTP,NRAY))
!         ENDDO
!         CALL SETLIN(0,0,7-MOD(NRAY-1,5))
!         CALL GPLOTP(GUX,GUY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)      
!      ENDDO


!     ----- WAVE NUMBER -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
            CALL WRCALK(NSTP,NRAY,RKPARA,RKPERP)
            GKX( NSTP+1,NRAY)=GUCLIP(XL)
            GKY1(NSTP+1,NRAY)=GUCLIP(RKPARA)
            GKY2(NSTP+1,NRAY)=GUCLIP(RKPERP)
         ENDDO
      ENDDO

      NRAY=1
         CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN1,GYMAX1)
         CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN2,GYMAX2)
      DO NRAY=2,NRAYMAX
         CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMIN1=MIN(GYMIN1,GYMIN)
         GYMAX1=MAX(GYMAX1,GYMAX)
         CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMIN2=MIN(GYMIN2,GYMIN)
         GYMAX2=MAX(GYMAX2,GYMAX)
      ENDDO

      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.7,16.1)
      CALL TEXT('k-para',6)
      CALL GDEFIN(13.7,23.7,9.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,2,7)
      CALL GFRAME
      GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,0,7-MOD(NRAY-1,5))
         CALL GPLOTP(GKX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

      CALL WRGPRM(1)
      CALL PAGEE
      RETURN
    END SUBROUTINE WRGRF11A


!     ***** RADIAL DEPENDENCE 1 *****

    SUBROUTINE WRGRF11B

      USE wrcomm
      USE plprof,ONLY: PL_MAG_OLD
      USE wrsub,ONLY: wrcalk
      IMPLICIT NONE
      INTEGER,PARAMETER:: NSRM=101
      REAL:: GX(NSTPMAX+1),GY(NSTPMAX+1)
      REAL:: GKX(NSTPMAX+1,NRAYM)
      REAL:: GKY1(NSTPMAX+1,NRAYM),GKY2(NSTPMAX+1,NRAYM)
      REAL:: GLCX(NSRM),GLCY(NSRM),GSCX(NSRM),GSCY(NSRM)
      INTEGER:: NSRMAX,NSU,NRAY,NSTP
      REAL(rkind):: RMAXL,RMINL,DTH,TH,XL,YL,ZL,RHON,RKPARA,RKPERP
      REAL:: GRRMAX,GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GYSTEP,GXORG
      REAL:: GYMIN1,GYMIN2,GYMAX1,GYMAX2,GYORG
      EXTERNAL GMNMX1,GQSCAL,PAGES,SETCHS,SETFNT,SETLIN,SETRGB,PAGEE
      EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP2
      EXTERNAL MOVE,TEXT
      REAL:: GUCLIP
      INTEGER:: NGULEN

      CALL PAGES
      CALL SETFNT(32)
      CALL SETCHS(0.25,0.)
      CALL SETLIN(0,2,7)

!     ----- TOROIDAL CROSS SECTION -----

      NSRMAX=101
      RMAXL=RR+RA
      RMINL=RR-RA
      DTH=2*PI/(NSRMAX-1)
      DO NSU=1,NSRMAX
         TH=DTH*(NSU-1)
         GLCX(NSU)=GUCLIP(RMAXL*COS(TH))
         GLCY(NSU)=GUCLIP(RMAXL*SIN(TH))
         GSCX(NSU)=GUCLIP(RMINL*COS(TH))
         GSCY(NSU)=GUCLIP(RMINL*SIN(TH))
      ENDDO
      IF(MOD(MDLWRG,2).EQ.0) THEN
         GRRMAX=GUCLIP(RMAXL)
         CALL GQSCAL(-GRRMAX,GRRMAX,GXMIN,GXMAX,GXSTEP)
         CALL GQSCAL(-GRRMAX,GRRMAX,GYMIN,GYMAX,GYSTEP)
      ELSE
         CALL GQSCAL(GUCLIP(RR-RA),GUCLIP(RR+RA),GXMIN,GXMAX,GXSTEP)
         CALL GQSCAL(GUCLIP(  -RA),GUCLIP(   RA),GYMIN,GYMAX,GYSTEP)
         GXORG=(INT(GXMIN/(2*GXSTEP))+1)*2*GXSTEP
      ENDIF

      CALL GDEFIN(1.7,11.7,1.0,11.0,GXMIN,GXMAX,GYMIN,GYMAX)

      CALL SETLIN(0,2,7)
      CALL GFRAME
      IF(MOD(MDLWRG,2).EQ.0) THEN
         CALL GSCALE(0.0,GXSTEP,0.0,GYSTEP,0.1,9)
         CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
         CALL GVALUE(0.0,0.0,0.0,2*GYSTEP,NGULEN(2*GYSTEP))
      ELSE
         CALL GSCALE(GXORG,GXSTEP,0.0,GYSTEP,0.1,9)
         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
         CALL GVALUE(0.0,0.0,0.0,2*GYSTEP,NGULEN(2*GYSTEP))
      ENDIF

      CALL SETLIN(0,2,4)
      CALL GPLOTP(GLCX,GLCY,1,NSRMAX,1,0,0,0)
      CALL GPLOTP(GSCX,GSCY,1,NSRMAX,1,0,0,0)

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            GX(NSTP+1)=GUCLIP(RAYS(1,NSTP,NRAY))
            GY(NSTP+1)=GUCLIP(RAYS(2,NSTP,NRAY))
         ENDDO
         CALL SETLIN(0,0,7-MOD(NRAY-1,5))
         CALL GPLOTP(GX,GY,1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

!     ----- WAVE NUMBER -----

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON)
            CALL WRCALK(NSTP,NRAY,RKPARA,RKPERP)
            GKX( NSTP+1,NRAY)=GUCLIP(RHON)
            GKY1(NSTP+1,NRAY)=GUCLIP(RKPARA*VC &
                                  /(2*PI*RAYIN(1,NRAY)*1.D6))
            GKY2(NSTP+1,NRAY)=GUCLIP(RKPERP)
         ENDDO
      ENDDO

      CALL GQSCAL(0.0,1.0,GXMIN,GXMAX,GXSTEP)
      NRAY=1
         CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN1,GYMAX1)
         CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN2,GYMAX2)
      DO NRAY=2,NRAYMAX
         CALL GMNMX1(GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMIN1=MIN(GYMIN1,GYMIN)
         GYMAX1=MAX(GYMAX1,GYMAX)
         CALL GMNMX1(GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY),1,GYMIN,GYMAX)
         GYMIN2=MIN(GYMIN2,GYMIN)
         GYMAX2=MAX(GYMAX2,GYMAX)
      ENDDO

      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.7,16.1)
      CALL TEXT('n-para',6)
      CALL GDEFIN(13.7,23.7,9.0,16.0,0.0,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,2,7)
      CALL GFRAME
      GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
      CALL GSCALE(0.0,GXSTEP,GYMIN,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,0,7-MOD(NRAY-1,5))
         CALL GPLOTP(GKX(1,NRAY),GKY1(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

      CALL GQSCAL(GYMIN2,GYMAX2,GYMIN,GYMAX,GYSTEP)
      CALL MOVE(13.7,8.1)
      CALL TEXT('k-perp',6)
      CALL GDEFIN(13.7,23.7,1.0,8.0,0.0,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,2,7)
      CALL GFRAME
      GYORG=(INT(GYMIN/(2*GYSTEP))+1)*2*GYSTEP
      CALL GSCALE(0.0,GXSTEP,GYORG,GYSTEP,0.1,9)
      CALL GVALUE(0.0,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))

      DO NRAY=1,NRAYMAX
         CALL SETLIN(0,0,7-MOD(NRAY-1,5))
         CALL GPLOTP(GKX(1,NRAY),GKY2(1,NRAY),1,NSTPMAX_NRAY(NRAY)+1,1,0,0,0)
      ENDDO
      CALL SETLIN(0,2,7)

      CALL WRGPRM(1)
      CALL PAGEE

      RETURN
    END SUBROUTINE WRGRF11B
END MODULE WRGOUT
