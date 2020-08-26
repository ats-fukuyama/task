! obgout.f90

MODULE obgout

  PRIVATE
  PUBLIC ob_gout

  INTERFACE
     FUNCTION GUCLIP(X)
       REAL(8):: X
       REAL(4):: GUCLIP
     END FUNCTION GUCLIP
     FUNCTION NGULEN(Y)
       REAL(4):: Y
       INTEGER:: NGULEN
     END FUNCTION NGULEN
  END INTERFACE

CONTAINS

!*********************** GRAPHIC DATA OUTPUT ************************

  SUBROUTINE ob_gout

    USE obcomm
    USE plgout
    IMPLICIT NONE
    CHARACTER(LEN=1)::  kid

1   WRITE(6,*) &
         '## INPUT GRAPH TYPE : 1,2   L,E:equilibrium P:profile X:end'
    READ(5,'(A1)',ERR=1,END=9000) kid
    CALL GUCPTL(kid)

    IF(KID.EQ.'1') CALL ob_grf1
    IF(KID.EQ.'2') CALL ob_grf2
    IF(KID.EQ.'L') CALL ob_gsube
    IF(KID.EQ.'E') CALL eqgsdd
    IF(KID.EQ.'P') CALL pl_gout
    IF(KID.EQ.'X') GOTO 9000

    GOTO 1

9000 RETURN
  END SUBROUTINE ob_gout

!     ***** POLOIDAL TRAJECTORY AND POWER *****

  SUBROUTINE ob_grf1

    USE obcomm
    USE obprep,ONLY: psipa
    USE plprof,ONLY: pl_rzsu,pl_mag_old,pl_mag_type,pl_mag
    USE libgrf
    IMPLICIT NONE
    INTEGER,PARAMETER:: nsu_m=501
    INTEGER,PARAMETER:: nx_max=101
    INTEGER,PARAMETER:: ny_max=101
    REAL(rkind),PARAMETER:: aspect_ob=0.D0
    REAL(rkind):: rgba(3,1),ltha(1)
    REAL(rkind),ALLOCATABLE:: rsu_temp(:),zsu_temp(:)
    REAL(rkind),ALLOCATABLE:: rsu(:),zsu(:)
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:),f(:,:)
    REAL(rkind):: dx,dy,x,y,xmin,xmax,ymin,ymax,dlx,dly,factor
    REAL:: line_width
    REAL,ALLOCATABLE:: reda(:),greena(:),bluea(:)
    INTEGER:: nsu,nsu_max,nstp,nobt,nx,ny

    rgba(1,1)=0.0D0
    rgba(2,1)=0.5D0
    rgba(3,1)=0.0D0
    ltha(1)=0.0D0

    !  ----- PLASMA BOUNDARY -----

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
         XSCALE_ZERO=0,YSCALE_ZERO=0)

!     ----- magnetic surface -----

    ALLOCATE(xg(1:nx_max),yg(1:ny_max),f(1:nx_max,1:ny_max))
    dx=(xmax-xmin)/(nx_max-1)
    dy=(ymax-ymin)/(ny_max-1)
    DO nx=1,nx_max
       xg(nx)=xmin+(nx-1)*dx
    ENDDO
    DO ny=1,ny_max
       yg(ny)=ymin+(ny-1)*dy
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
         NOTITLE=1,NOFRAME=1,NOINFO=1,ASPECT=aspect_ob, &
         LINE_RGB=rgba,LINE_THICKNESS=ltha)

    !  ----- PLASMA SURFACE -----

    CALL SETRGB(0.0,0.0,1.0)
    CALL MOVE2D(REAL(rsu(1)),REAL(zsu(1)))
    DO nsu=2,nsu_max
       CALL DRAW2D(REAL(rsu(nsu)),REAL(zsu(nsu)))
    ENDDO
    CALL DRAW2D(REAL(rsu(1)),REAL(zsu(1)))

    !  ----- OBT TRAJECTORY -----

    ALLOCATE(reda(nobt_max),greena(nobt_max),bluea(nobt_max))
    IF(nobt_max.EQ.1) THEN
       reda(1)=1.0; greena(1)=0.0; bluea(1)=0.0
    ELSE
       DO nobt=1,nobt_max/2
          factor=REAL(nobt-1)/REAL(nobt_max/2)
          reda(nobt)=1.0-factor; greena(nobt)=factor; bluea(nobt)=0.0
       END DO
       DO nobt=nobt_max/2+1,nobt_max
          factor=REAL(nobt-nobt_max/2-1)/REAL(nobt_max/2)
          reda(nobt)=0.0; greena(nobt)=1.0-factor; bluea(nobt)=factor
       END DO
    END IF
    CALL setmks(3,0.2)
    DO nobt=1,nobt_max
!       WRITE(6,'(A,I5,1P3E12.4)') &
!            'rgb:',nobt,reda(nobt),greena(nobt),bluea(nobt)
       CALL SETRGB(reda(nobt),greena(nobt),bluea(nobt))
       CALL MARK2D(REAL(rr_ob(0,nobt)),REAL(zz_ob(0,nobt)))
       DO nstp=1,nstp_max_nobt(nobt)
          CALL DRAW2D(REAL(rr_ob(nstp,nobt)),REAL(zz_ob(nstp,nobt)))
       ENDDO
    ENDDO
    DEALLOCATE(reda,greena,bluea)

    CALL GRD2D_FRAME_END

    CALL SETLNW(line_width)
    CALL SETRGB(0.0,0.0,0.0)

!    CALL OBGPRM
    CALL PAGEE
    RETURN
  END SUBROUTINE ob_grf1

!     ***** 1D plots *****

  SUBROUTINE ob_grf2

    USE obcomm
    USE libgrf
    IMPLICIT NONE
    REAL(rkind):: fx(nstp_max+1),fy(nstp_max+1,1)
    INTEGER:: nstp,nobt,nxmax
    REAL(rkind):: line_rgb(3,1),line_mark_size(1)
    INTEGER:: line_mark(1),line_mark_step(1)
    REAL(rkind):: xmin,xmax,fmin,fmax,x,f

    line_rgb(1,1)=1.D0
    line_rgb(2,1)=0.D0
    line_rgb(3,1)=0.D0
    line_mark_size(1)=0.3
    line_mark(1)=-3

    xmin=time_ob(0,1)
    xmax=time_ob(nstp_max_nobt(1),1)
    DO nobt=2,nobt_max
       xmin=MIN(xmin,time_ob(0,nobt))
       xmax=MAX(xmax,time_ob(nstp_max_nobt(nobt),nobt))
    END DO

    CALL PAGES
    fmin=zetab_ob(0,1)
    fmax=zetab_ob(0,1)
    DO nobt=1,nobt_max
       DO nstp=1,nstp_max_nobt(nobt)
          fmin=MIN(fmin,zetab_ob(nstp,nobt))
          fmax=MAX(fmax,zetab_ob(nstp,nobt))
       END DO
    END DO
    fmin=fmin*180.D0/Pi
    fmax=fmax*180.D0/Pi
    WRITE(6,'(A,4ES12.4)') 'xf:',xmin,xmax,fmin,fmax
    CALL GRD1D_FRAME_START(1,xmin,xmax,fmin,fmax, &
                          '@zetab(s) [deg]@')
    DO nobt=1,nobt_max
       CALL SETMKS(3,0.2)
       CALL SETRGB(1.0-REAL(nobt)/REAL(nobt_max),0.0,REAL(nobt)/REAL(nobt_max))
       x=time_ob(0,nobt)
       f=zetab_ob(0,nobt)*180.D0/Pi
       CALL MARK2D(GUCLIP(x),GUCLIP(f))
       DO nstp=1,nstp_max_nobt(nobt)
          x=time_ob(nstp,nobt)
          f=zetab_ob(nstp,nobt)*180.D0/Pi
          CALL DRAW2D(GUCLIP(x),GUCLIP(f))
       END DO
    END DO
    CALL GRD1D_FRAME_END

    DO nobt=1,nobt_max
       CALL PAGES
       nxmax=nstp_max_nobt(nobt)+1
       line_mark_step(1)=nxmax ! mark only at first point
       fx(1:nxmax)=time_ob(0:nxmax-1,nobt)

       fy(1:nxmax,1)=zetab_ob(0:nxmax-1,nobt)*180.D0/Pi
       CALL GRD1D(1,fx,fy,nstp_max,nxmax,1,'@zetab(s) [deg]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=thetab_ob(0:nxmax-1,nobt)*180.D0/Pi
       CALL GRD1D(2,fx,fy,nstp_max,nxmax,1,'@thetab(s) [deg]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=psip_ob(0:nxmax-1,nobt)
       CALL GRD1D(3,fx,fy,nstp_max,nxmax,1,'@psip(s) [Bm^2]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=rhopara_ob(0:nxmax-1,nobt)
       CALL GRD1D(4,fx,fy,nstp_max,nxmax,1,'@rhopara(s) [m]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       CALL PAGES

       fy(1:nxmax,1)=pzeta_ob(0:nxmax-1,nobt)/AEE
       CALL GRD1D(1,fx,fy,nstp_max,nxmax,1,'@pzeta(s) [eV]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=ptheta_ob(0:nxmax-1,nobt)
       CALL GRD1D(2,fx,fy,nstp_max,nxmax,1,'@ptheta(s) [eV]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=babs_ob(0:nxmax-1,nobt)
       CALL GRD1D(3,fx,fy,nstp_max,nxmax,1,'@babs(s) [T]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=phi_ob(0:nxmax-1,nobt)
       CALL GRD1D(4,fx,fy,nstp_max,nxmax,1,'@phi(s) [V]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       CALL PAGES

       fy(1:nxmax,1)=vpara_ob(0:nxmax-1,nobt)
       CALL GRD1D(1,fx,fy,nstp_max,nxmax,1,'@vpara(s) [m/s]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=vperp_ob(0:nxmax-1,nobt)
       CALL GRD1D(2,fx,fy,nstp_max,nxmax,1,'@vperp(s) [m/s]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=psit_ob(0:nxmax-1,nobt)
       CALL GRD1D(3,fx,fy,nstp_max,nxmax,1,'@psit(s) [Bm^2]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=zeta_ob(0:nxmax-1,nobt)
       CALL GRD1D(4,fx,fy,nstp_max,nxmax,1,'@zeta(s) [deg]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       CALL PAGES

       fy(1:nxmax,1)=rr_ob(0:nxmax-1,nobt)
       CALL GRD1D(1,fx,fy,nstp_max,nxmax,1,'@rr(s) [m]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=zz_ob(0:nxmax-1,nobt)
       CALL GRD1D(2,fx,fy,nstp_max,nxmax,1,'@zz(s) [m]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=rs_ob(0:nxmax-1,nobt)
       CALL GRD1D(3,fx,fy,nstp_max,nxmax,1,'@rs(s) [m]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       fy(1:nxmax,1)=theta_ob(0:nxmax-1,nobt)*180.D0/Pi
       CALL GRD1D(4,fx,fy,nstp_max,nxmax,1,'@theta(s) [deg]@', &
                  NLMAX=1, &
                  LINE_RGB=line_rgb,LINE_MARK_SIZE=line_mark_size, &
                  LINE_MARK=line_mark,LINE_MARK_STEP=line_mark_step)

       CALL PAGEE
       
    END DO
    RETURN
  END SUBROUTINE ob_grf2

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
