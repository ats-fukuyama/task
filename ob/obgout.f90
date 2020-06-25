! obgout.f90

MODULE obgout

  PRIVATE
  PUBLIC ob_gout

CONTAINS

!*********************** GRAPHIC DATA OUTPUT ************************

  SUBROUTINE ob_gout

    USE obcomm
    USE plgout
    IMPLICIT NONE
    CHARACTER(LEN=1)::  kid

1   WRITE(6,*) &
         '## INPUT GRAPH TYPE : 1   E:equilibrium P:profile X:end'
    READ(5,'(A1)',ERR=1,END=9000) kid
    CALL GUCPTL(kid)

    IF(KID.EQ.'1') CALL ob_grf1
    IF(KID.EQ.'E') CALL ob_gsube
    IF(KID.EQ.'P') CALL pl_gout
    IF(KID.EQ.'X') GOTO 9000

    GOTO 1

9000 RETURN
  END SUBROUTINE ob_gout

!     ***** POLOIDAL TRAJECTORY AND POWER *****

  SUBROUTINE ob_grf1

    USE obcomm
    USE plprof,ONLY: pl_rzsu,pl_mag_old,pl_mag_type,pl_mag
    USE libgrf
    IMPLICIT NONE
    INTEGER,PARAMETER:: nsu_m=501
    INTEGER,PARAMETER:: nx_max=101
    INTEGER,PARAMETER:: ny_max=101

    REAL(rkind),ALLOCATABLE:: rsu_temp(:),zsu_temp(:)
    REAL(rkind),ALLOCATABLE:: rsu(:),zsu(:)
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:),f(:,:)
    REAL(rkind):: dx,dy,x,y,xmin,xmax,ymin,ymax,dlx,dly,factor
    REAL(rkind),PARAMETER:: aspect_ob=0.D0
    REAL:: line_width
    REAL,ALLOCATABLE:: reda(:),greena(:),bluea(:)
    INTEGER:: nsu,nsu_max,nstp,nobt,nx,ny

!  ----- PLASMA BOUNDARY -----

    dlx=0.02D0*(xmax-xmin)
    dly=0.02D0*(ymax-ymin)

    ALLOCATE(rsu_temp(nsu_m),zsu_temp(nsu_m))
    CALL pl_rzsu(rsu_temp,zsu_temp,nsu_m,nsu_max)
    ALLOCATE(rsu(nsu_max+1),zsu(nsu_max+1))
    DO nsu=1,nsu_max
       rsu(nsu)=rsu_temp(nsu)
       zsu(nsu)=zsu_temp(nsu)
    END DO
    rsu(nsu_max+1)=rsu(1)
    zsu(nsu_max+1)=zsu(1)
    DEALLOCATE(rsu_temp,zsu_temp)

    xmin=rsu(1)
    xmax=rsu(1)
    ymin=zsu(1)
    ymax=zsu(1)
    DO nsu=2,nsu_max
       xmin=MIN(xmin,rsu(nsu))
       xmax=MAX(xmax,rsu(nsu))
       ymin=MIN(ymin,zsu(nsu))
       ymax=MAX(ymax,zsu(nsu))
    END DO
       
    dlx=0.02D0*(xmax-xmin)
    dly=0.02D0*(ymax-ymin)

    CALL pages
    CALL INQLNW(line_width)
    CALL SETLNW(0.0)

    CALL GRD2D_FRAME_START(0,xmin-dlx,xmax+dlx,ymin-dly,ymax+dly, &
         '@orbit in (r,z)@', ASPECT=aspect_ob,NOINFO=1, & !
         GPXMIN=2.D0,XSCALE_ZERO=0,YSCALE_ZERO=0)

!     ----- magnetic surface -----

    ALLOCATE(xg(1:nx_max),yg(1:ny_max),f(1:nx_max,1:ny_max))
    dx=(xmax-xmin)/nx_max
    dy=(ymax-ymin)/ny_max
    DO nx=1,nx_max
       xg(nx)=xmin+(nx-1)*dx
    ENDDO
    DO ny=1,ny_max
       yg(ny)=ymin+(ny-1)*dx
    ENDDO
    DO ny=1,ny_max
       y=yg(ny)
       DO nx=1,nx_max
          x=xg(nx)
          CALL PL_MAG_OLD(x,0.D0,y,f(nx,ny))
       ENDDO
    ENDDO
    CALL GRD2D(0,xg,yg,f,nx_max,nx_max,ny_max, &
         XMIN=xmin-dlx,XMAX=xmax+dlx,YMIN=ymin-dly,YMAX=ymax+dly, &
         NOTITLE=1,NOFRAME=1,NOINFO=1)

    !  ----- PLASMA SURFACE -----

    CALL SETRGB(0.0,0.0,1.0)
    CALL MOVE2D(rsu(1),zsu(1))
    DO nsu=1,nsu_max
       CALL DRAW2D(rsu(nsu+1),zsu(nsu+1))
    ENDDO
    ALLOCATE(rsu(nsu_max+1),zsu(nsu_max+1))

    !  ----- OBT TRAJECTORY -----

    ALLOCATE(reda(nobt_max),greena(nobt_max),bluea(nobt_max))
    IF(nobt_max.EQ.1) THEN
       reda(1)=1.0; greena(1)=0.0; bluea(1)=0.0
    ELSE
       DO nobt=1,nobt_max/2
          factor=REAL(nobt-1)/REAL(nobt_max/2)
          reda(1)=1.0-factor; greena(1)=factor; bluea(1)=0.0
       END DO
       DO nobt=nobt_max/2+1,nobt_max
          factor=REAL(nobt-nobt_max/2-1)/REAL(nobt_max/2)
          reda(1)=0.0; greena(1)=1.0-factor; bluea(1)=factor
       END DO
    END IF
    DO nobt=1,nobt_max
       CALL SETRGB(reda(nobt),greena(nobt),bluea(nobt))
       CALL MOVE2D(rr_ob(0,nobt),zz_ob(0,nobt))
       DO nstp=1,nstp_max_nobt(nobt)
          CALL DRAW2D(rr_ob(nstp,nobt),zz_ob(nstp,nobt))
       ENDDO
    ENDDO
    DEALLOCATE(reda,greena,bluea)

    CALL GRD2D_FRAME_END

    CALL SETLNW(line_width)
    CALL SETRGB(0.0,0.0,0.0)

    CALL OBGPRM
    CALL PAGEE
    RETURN
  END SUBROUTINE ob_grf1

!     ***** DRAW PARAMETERS *****

    SUBROUTINE OBGPRM

      USE obcomm
      IMPLICIT NONE

!     ----- DRAW PARAMETERS -----

      CALL SETLIN(0,0,7)
      CALL OBPRMT(1.0,18.0, 'BB   =',6, BB  ,'(F8.3)',8)
      CALL OBPRMT(1.0,17.6, 'RR   =',6, RR  ,'(F8.3)',8)
      CALL OBPRMT(1.0,17.2, 'RA   =',6, RA  ,'(F8.3)',8)
      CALL OBPRMT(1.0,16.8, 'RKAP =',6, RKAP,'(F8.3)',8)
      CALL OBPRMT(5.0,18.0, 'Q0   =',6, Q0  ,'(F8.3)',8)
      CALL OBPRMT(5.0,17.6, 'QA   =',6, QA  ,'(F8.3)',8)
      CALL OBPRMT(9.0,18.0, 'PA(1)  =',8,PA(1),'(F8.4)',8)
      CALL OBPRMT(9.0,17.6, 'PZ(1)  =',8,PZ(1),'(F8.4)',8)
      CALL OBPRMT(9.0,17.2, 'PN(1)  =',8,PN(1),'(F8.4)',8)
      CALL OBPRMT(9.0,16.8, 'PNS(1) =',8,PNS(1),'(F8.4)',8)
      CALL OBPRMT(13.0,18.0,'PTPR(1)=',8,PTPR(1),'(F8.4)',8)
      CALL OBPRMT(13.0,17.6,'PTPP(1)=',8,PTPP(1),'(F8.4)',8)
      CALL OBPRMT(13.0,17.2,'PTS(1) =',8,PTS(1),'(F8.4)',8)
      CALL OBPRMT(17.0,18.0,'PA(2)  =',8,PA(2),'(F8.4)',8)
      CALL OBPRMT(17.0,17.6,'PZ(2)  =',8,PZ(2),'(F8.4)',8)
      CALL OBPRMT(17.0,17.2,'PN(2)  =',8,PN(2),'(F8.4)',8)
      CALL OBPRMT(17.0,16.8,'PNS(2)= ',8,PNS(2),'(F8.4)',8)
      CALL OBPRMT(21.0,18.0,'PTPR(2)=',8,PTPR(2),'(F8.4)',8)
      CALL OBPRMT(21.0,17.6,'PTPP(2)=',8,PTPP(2),'(F8.4)',8)
      CALL OBPRMT(21.0,17.2,'PTS(2) =',8,PTS(2),'(F8.4)',8)
      CALL PAGEE
      RETURN
    END SUBROUTINE OBGPRM

!***********************************************************************

    SUBROUTINE OBPRMT(X,Y,KTEXT,ITEXT,D,KFORM,IFORM)

      REAL(4),INTENT(IN):: X,Y
      INTEGER,INTENT(IN):: ITEXT,IFORM
      CHARACTER(LEN=*),INTENT(IN):: KTEXT,KFORM
      REAL(8),INTENT(IN):: D

      CALL MOVE(X,Y)
      CALL TEXT(KTEXT,ITEXT)
      CALL NUMBD(D,KFORM,IFORM)
      RETURN
    END SUBROUTINE OBPRMT

! *** equilibrimu quantities plot ***
    
    SUBROUTINE ob_gsube
      USE obcomm
      USE obprep
      USE libgrf
      IMPLICIT NONE
      INTEGER,PARAMETER:: npsip_max=50
      INTEGER,PARAMETER:: nchi_max=65
      REAL(rkind),ALLOCATABLE:: psip_l(:),chi_l(:),gx(:)
      REAL(rkind),ALLOCATABLE:: f1(:,:),f2a(:,:),f2b(:,:),f2c(:,:)
      INTEGER:: nchi,npsip,ierr
      REAL(rkind):: dchi,dpsip,psip0

      ALLOCATE(psip_l(npsip_max),chi_l(nchi_max),gx(npsip_max))
      ALLOCATE(f1(npsip_max,2),f2a(nchi_max,npsip_max))
      ALLOCATE(f2b(nchi_max,npsip_max),f2c(nchi_max,npsip_max))

      psip0=0.D0
      dpsip=(psipa-psip0)/npsip_max
      dchi=2.D0*Pi/(nchi_max-1)
      DO npsip=1,npsip_max
         psip_l(npsip)=(psip0+dpsip*npsip)
         gx(npsip)=psip_l(npsip)/psipa
      END DO
      DO nchi=1,nchi_max
         chi_l(nchi)=dchi*(nchi-1)
      END DO

      CALL pages

      DO npsip=1,npsip_max
         DO nchi=1,nchi_max
            CALL cal_rr_pos(chi_l(nchi),psip_l(npsip),f2a(nchi,npsip),ierr)
            CALL cal_zz_pos(chi_l(nchi),psip_l(npsip),f2b(nchi,npsip),ierr)
         END DO
      END DO
      
      CALL grd2d(1,chi_l,psip_l,f2a,nchi_max,nchi_max,npsip_max, &
                 '@rr_pos(chi,psip)@',0,XMAX=chi_l(nchi_max))

      CALL grd2d(2,chi_l,psip_l,f2b,nchi_max,nchi_max,npsip_max, &
                 '@zz_pos(chi,psip)@',0,XMAX=chi_l(nchi_max))

      DO npsip=1,npsip_max
         DO nchi=1,nchi_max
            CALL cal_bdb_pos(chi_l(nchi),psip_l(npsip),f2a(nchi,npsip), &
                 f2b(nchi,npsip),f2c(nchi,npsip),ierr)
         END DO
      END DO
      
      CALL grd2d(3,chi_l,psip_l,f2a,nchi_max,nchi_max,npsip_max, &
                 '@bb_pos(chi,psip)@',0,XMAX=chi_l(nchi_max))

      CALL pagee

      CALL pages

      DO npsip=1,npsip_max
         CALL cal_qps_pos(psip_l(npsip),f1(npsip,1),f1(npsip,2),ierr)
         f1(npsip,2)=f1(npsip,2)*psip_l(npsip)/f1(npsip,1)
      END DO
      CALL grd1d(1,gx,f1,npsip_max,npsip_max,2,'@qps(psipn)@',0)
      
      DO npsip=1,npsip_max
         CALL cal_rbps_pos(psip_l(npsip),f1(npsip,1),f1(npsip,2),ierr)
         f1(npsip,2)=f1(npsip,2)*psipa
      END DO
      CALL grd1d(2,gx,f1,npsip_max,npsip_max,2,'@rbps(psipn)@',0)
      
      DO npsip=1,npsip_max
         CALL cal_ritps_pos(psip_l(npsip),f1(npsip,1),f1(npsip,2),ierr)
         f1(npsip,2)=f1(npsip,2)*psipa
      END DO
      CALL grd1d(3,gx,f1,npsip_max,npsip_max,2,'@ritps(psipn)@',0)
      
      DO npsip=1,npsip_max
         CALL cal_psit_pos(psip_l(npsip),f1(npsip,1),f1(npsip,2),ierr)
         f1(npsip,1)=f1(npsip,1)/psita
      END DO
      CALL grd1d(4,gx,f1,npsip_max,npsip_max,2,'@psit(psipn)@',0)
      
      CALL pagee

      ! page 3
      
      CALL pages
      
      CALL grd2d(1,chi_l,psip_l,f2a,nchi_max,nchi_max,npsip_max, &
                 '@bb_pos(chi,psip)@',0,XMAX=chi_l(nchi_max))
      CALL grd2d(3,chi_l,psip_l,f2b,nchi_max,nchi_max,npsip_max, &
                 '@db_dchi(chi,psip)@',0,XMAX=chi_l(nchi_max))
      CALL grd2d(4,chi_l,psip_l,f2c,nchi_max,nchi_max,npsip_max, &
                 '@db_dpsip(chi,psip)@',0,XMAX=chi_l(nchi_max))
      CALL pagee
      
    END SUBROUTINE ob_gsube
END MODULE OBGOUT
