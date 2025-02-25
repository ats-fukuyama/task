! wmeign.f90

MODULE wmeign

  PRIVATE
  PUBLIC wm_am0d,wm_am1d,wm_am2d,wm_scan,wm_eign,wm_wout

CONTAINS

!     ***** SINGLE CALCULATION *****

  SUBROUTINE wm_am0d(KID,LINE)

    USE wmcomm
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE wmemfp,ONLY: wm_bfield
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80):: LINE
    CHARACTER(LEN=1):: KID
    INTEGER:: MODE
    REAL(rkind):: AMPL
    EXTERNAL:: GUFLSH

    MODE=0
1   CONTINUE
    IF(NRANK.EQ.0) THEN
2      CONTINUE
       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4)') &
            ' ## FR,FI=',FRINI,FIINI
       WRITE(6,'(A)') ' ## INPUT : KID or FR,FI'
       CALL TASK_KLIN(LINE,KID,MODE,wm_parm)
       IF(MODE.EQ.0) &
            READ(LINE,*,ERR=2,END=2) FRINI,FIINI
    ENDIF
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE.EQ.2) THEN
       CALL wm_broadcast
       GOTO 1
    ENDIF
    IF(MODE.EQ.3) GOTO 1

    CALL mtx_broadcast1_real8(FRINI)
    CALL mtx_broadcast1_real8(FIINI)

    CALL DIAMIN(FRINI,FIINI,AMPL)
         
    CALL wm_bfield

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_am0d

!     ***** 1D PLOT OF AMPLITUDE *****

  SUBROUTINE wm_am1d(KID,LINE)

    USE wmcomm
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE libkio
    USE libmpi
    IMPLICIT NONE
    REAL,ALLOCATABLE::  GX(:),GF(:,:)
    CHARACTER(LEN=80)::  LINE
    CHARACTER(LEN=1):: KID
    INTEGER:: MODE,NGF
    REAL(rkind):: DELTFR,AMPMIN,FR,FI,AMPL
    EXTERNAL GUFLSH

    MODE=0
1   CONTINUE
    IF(NRANK.EQ.0) THEN
2      CONTINUE
       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4,I5,1PE12.4)') &
            ' ## FRMIN,FRMAX,NGFMAX,FI0=',FRMIN,FRMAX,NGFMAX,FI0
       WRITE(6,'(A)') ' ## INPUT : KID or FRMIN,FRMAX,NGFMAX,FI0'
       CALL TASK_KLIN(LINE,KID,MODE,wm_parm)
       IF(MODE.EQ.0) &
            READ(LINE,*,ERR=2,END=2) FRMIN,FRMAX,NGFMAX,FI0
    ENDIF
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE.EQ.2) THEN
       CALL wm_broadcast
       GOTO 1
    ENDIF
    IF(MODE.EQ.3) GOTO 1

    CALL mtx_broadcast1_integer(NGFMAX)
    CALL mtx_broadcast1_real8(FRMIN)
    CALL mtx_broadcast1_real8(FRMAX)
    CALL mtx_broadcast1_real8(FI0)

    IF(ALLOCATED(GX)) DEALLOCATE(GX,GF)
    ALLOCATE(GX(NGFMAX),GF(NGFMAX,1))

    DELTFR=(FRMAX-FRMIN)/(NGFMAX-1)
    AMPMIN=1.D30
    DO NGF=1,NGFMAX
       FR=FRMIN+(NGF-1)*DELTFR
       FI=FI0
       CALL DIAMIN(FR,FI,AMPL)
       IF(NRANK.EQ.0) THEN
          IF(LISTEG.GE.1) THEN
             WRITE(6,'(A,1P3E12.4)') &
                  '      FR,FI,AMPL = ',FR,FI,AMPL
             CALl GUFLSH
          ENDIF
       ENDIF
       GX(NGF)=GUCLIP(FR)
       GF(NGF,1)=GUCLIP(LOG10(AMPL))
       IF(AMPL.LE.AMPMIN) THEN
          FRINI=FR
          FIINI=FI
          AMPMIN=AMPL
       ENDIF
    ENDDO
    IF(NRANK.EQ.0) THEN
       WRITE(6,'(A,1P3E12.4)') '      FR,FI,AMPMIN=', &
            FRINI,FIINI,AMPMIN
       CALL GUFLSH
    ENDIF

    CALL WMEG1D(GX,GF,NGFMAX,NGFMAX,1,GUCLIP(FI0))

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_am1d

  SUBROUTINE WMEG1D(GX,GZ,NGXM,NGXMAX,NGYMAX,GFI0)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NGXM,NGXMAX,NGYMAX
    REAL,INTENT(IN):: GX(NGXMAX),GZ(NGXM,NGYMAX),GFI0
    INTEGER:: NGY
    REAL:: GXMIN1,GXMAX1,GZMIN1,GZMAX1
    REAL:: GXMIN,GXMAX,GXSCL,GZMIN,GZMAX,GZSCL
    INTERFACE
       FUNCTION NGULEN(X)
         USE bpsd_kinds
         REAL:: X
         INTEGER:: NGULEN
       END FUNCTION NGULEN
    END INTERFACE
    EXTERNAL:: GMNMX2,GQSCAL,PAGES,SETCHS,SETLIN,MOVE,TEXT,NUMBR,PAGEE
    EXTERNAL:: GDEFIN,GFRAME,GSCALE,GVALUE,GSCALL,GVALUL,GPLOTP


    GXMIN1=GX(1)
    GXMAX1=GX(NGXMAX)
    CALL GMNMX2(GZ,NGXMAX,1,NGXMAX,1,1,NGYMAX,1,GZMIN1,GZMAX1)
    IF(GZMIN1.LT.-10.0) GZMIN1=-10.0
    CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSCL)
    CALL GQSCAL(GZMIN1,GZMAX1,GZMIN,GZMAX,GZSCL)

    CALL PAGES
    CALL SETCHS(0.3,0.0)
    CALL SETLIN(0,0,4)
    CALL GDEFIN( 4.0,24.0, 3.0,16.0,GXMIN,GXMAX,GZMIN,GZMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,  GXSCL,0.0,0.0,0.2,9)
    GXMIN=NINT(GXMIN/(2*GXSCL)+1)*2*GXSCL
    CALL GVALUE(GXMIN,2*GXSCL,0.0,0.0,NGULEN(2*GXSCL))  
    CALL GSCALL(0.0,0,  GZMIN,3,0.2,9)
    GZMIN=NINT(GZMIN/(2*GZSCL)+1)*2*GZSCL
    CALL GVALUL(0.0,0,  GZMIN,1,1)
    DO NGY=1,NGYMAX
       CALL SETLIN(0,0,7-MOD(NGY-1,5))
       CALL GPLOTP(GX,GZ(1,NGY),1,NGXMAX,1,0,0,0) 
    ENDDO

    CALL MOVE(10.0,0.5)
    CALL TEXT('Re(f) [MHz]',11)
    CALL MOVE(4.0,17.0)
    CALL TEXT('AMP1D',5)
    CALL MOVE(7.0,17.0)
    CALL TEXT('Im(f)=',6)
    CALL NUMBR(GFI0,'(1PE12.4)',12)
    CALL MOVE(13.0,17.0)
    CALL TEXT('min=',4)
    CALL NUMBR(GZMIN,'(1PE12.4)',12)
    CALL TEXT('   max=',7)
    CALL NUMBR(GZMAX,'(1PE12.4)',12)

    CALL PAGEE
    RETURN
  END SUBROUTINE WMEG1D

!     ***** 2D PLOT OF AMPLITUDE *****

  SUBROUTINE wm_am2d(KID,LINE)

    USE wmcomm
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE libkio
    USE libmpi
    IMPLICIT NONE
    REAL,ALLOCATABLE:: GX(:),GY(:),GZ(:,:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    INTEGER,PARAMETER:: NINFMAX=5
    REAL:: GFINF(3,NINFMAX)
    CHARACTER(LEN=80):: LINE
    CHARACTER(LEN=1):: KID
    INTEGER:: MODE,NGX,NGY,NINF,NINF1
    INTEGER:: NGXN,NGXP,NGYN,NGYP
    REAL(rkind):: DELTFR,DELTFI,FR,FI,AMPL,AMPMIN
    REAL:: GAMPL
    EXTERNAL:: GUFLSH

    MODE=0
1   CONTINUE
    IF(NRANK.EQ.0) THEN
2      CONTINUE
       IF(MODE.EQ.0) THEN
          WRITE(6,'(A)') &
               ' ## FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX='
          WRITE(6,'(3X,1P2E12.4,I5,1P2E12.4,I5)') &
               FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX
       ENDIF
       WRITE(6,'(A,A)') ' ## INPUT : KID or ', &
            'FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX'
       CALL TASK_KLIN(LINE,KID,MODE,wm_parm)
       IF(MODE.EQ.0) &
            READ(LINE,*,ERR=2,END=2) &
            FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX
    ENDIF
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE.EQ.2) THEN
       CALL wm_broadcast
       GOTO 1
    ENDIF
    IF(MODE.EQ.3) GOTO 1

    CALL mtx_broadcast1_integer(NGXMAX)
    CALL mtx_broadcast1_real8(FRMIN)
    CALL mtx_broadcast1_real8(FRMAX)
    CALL mtx_broadcast1_integer(NGYMAX)
    CALL mtx_broadcast1_real8(FIMIN)
    CALL mtx_broadcast1_real8(FIMAX)

    IF(ALLOCATED(GX)) DEALLOCATE(GX,GY,GZ,KA)
    ALLOCATE(GX(NGXMAX),GY(NGYMAX),GZ(NGXMAX,NGYMAX),KA(8,NGXMAX,NGYMAX))

    DELTFR=(FRMAX-FRMIN)/(NGXMAX-1)
    DELTFI=(FIMAX-FIMIN)/(NGYMAX-1)
    DO NGX=1,NGXMAX
       FR=FRMIN+(NGX-1)*DELTFR
       GX(NGX)=GUCLIP(FR)
       DO NGY=1,NGYMAX
          FI=FIMIN+(NGY-1)*DELTFI
          CALL DIAMIN(FR,FI,AMPL)
          IF(NRANK.EQ.0) THEN
             IF(LISTEG.GE.1) THEN
                WRITE(6,'(A,1P3E12.4)') '      FR,FI,AMPL = ',FR,FI,AMPL
                CALL GUFLSH
             ENDIF
          ENDIF
          GY(NGY)=GUCLIP(FI)
          GZ(NGX,NGY)=GUCLIP(AMPL)
       ENDDO
    ENDDO

    DO NINF=1,NINFMAX
       GFINF(1,NINF)=0.0
       GFINF(2,NINF)=0.0
       GFINF(3,NINF)=1.E30
    ENDDO
    DO NGX=1,NGXMAX
       DO NGY=1,NGYMAX
          GAMPL=GZ(NGX,NGY)
          NGXN=MAX(1,NGX-1)
          NGXP=MIN(NGXMAX,NGX+1)
          NGYN=MAX(1,NGY-1)
          NGYP=MIN(NGYMAX,NGY+1)
          IF(GZ(NGXN,NGY).GE.GAMPL.AND. &
             GZ(NGXP,NGY).GE.GAMPL.AND. &
             GZ(NGX,NGYN).GE.GAMPL.AND. &
             GZ(NGX,NGYP).GE.GAMPL) THEN
             DO NINF=1,NINFMAX
                IF(GFINF(3,NINF).GT.GAMPL) THEN
                   DO NINF1=NINFMAX,NINF+1,-1
                      GFINF(1,NINF1)=GFINF(1,NINF1-1)
                      GFINF(2,NINF1)=GFINF(2,NINF1-1)
                      GFINF(3,NINF1)=GFINF(3,NINF1-1)
                   ENDDO
                   GFINF(1,NINF)=GX(NGX)
                   GFINF(2,NINF)=GY(NGY)
                   GFINF(3,NINF)=GAMPL
                   GOTO 1000
                ENDIF
             ENDDO
1000         CONTINUE
          ENDIF
       ENDDO
    ENDDO
    IF(GFINF(3,1).LT.1.E30) THEN
       FRINI=DBLE(GFINF(1,1))
       FIINI=DBLE(GFINF(2,1))
       AMPMIN=DBLE(GFINF(3,1))
    ENDIF
    IF(NRANK.EQ.0) THEN
       WRITE(6,'(A,1P3E12.4)') &
            '      FR,FI,AMPMIN=',FRINI,FIINI,AMPMIN
       CALL GUFLSH
    ENDIF

    CALL WMEG2D(GX,GY,GZ,NGXMAX,NGXMAX,NGYMAX,KA,GFINF,NINFMAX)

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_am2d

  SUBROUTINE WMEG2D(GX,GY,GZ,NGXM,NGXMAX,NGYMAX,KA,GFINF,NINFMAX)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NGXM,NGXMAX,NGYMAX,NINFMAX
    REAL,INTENT(IN):: GX(NGXMAX),GY(NGYMAX),GZ(NGXM,NGYMAX)
    REAL,INTENT(IN):: GFINF(3,NINFMAX)
    INTEGER,INTENT(IN):: KA(8,NGXM,NGYMAX)
    REAL,PARAMETER:: RGBS(3,5)=RESHAPE([0.0,1.0,1.0, &
                                           0.0,0.0,1.0, &
                                           0.0,1.0,0.0, &
                                           1.0,0.5,0.0, &
                                           1.0,0.0,0.0],[3,5])
    INTEGER:: IGMAX,IGMIN,IPRD,I,NINF
    REAL:: GXMIN,GXMAX,GXMIN1,GXMAX1,GXSCL,GXORG
    REAL:: GYMIN,GYMAX,GYMIN1,GYMAX1,GYSCL,GYORG
    REAL:: GZMIN,GZMAX,GZLMIN,GZLMAX,GZORG,GZSTEP,GL
    INTERFACE
       FUNCTION NGULEN(X)
         USE bpsd_kinds
         REAL:: X
         INTEGER:: NGULEN
       END FUNCTION NGULEN
    END INTERFACE
    EXTERNAL:: GQSCAL,GMNMX2,PAGES,SETFNT,SETCHS,SETLIN,SETLNW
    EXTERNAL:: GDEFIN,GFRAME,GSCALE,GVALUE,SETRGB,CONTQ1
    EXTERNAL:: SETVEW,OFFVEW,MOVE,TEXT,NUMBR,NUMBI,PAGEE
    
    GXMIN=GX(1)
    GXMAX=GX(NGXMAX)
    GYMIN=GY(1)
    GYMAX=GY(NGYMAX)
    CALL GQSCAL(GXMIN,GXMAX,GXMIN1,GXMAX1,GXSCL)
    CALL GQSCAL(GYMIN,GYMAX,GYMIN1,GYMAX1,GYSCL)

    IF(GXMIN*GXMAX.LE.0.0) THEN
       GXORG=0.0
    ELSE
       GXORG=(INT(GXMIN/(2*GXSCL))+1)*2*GXSCL
    ENDIF
    IF(GYMIN*GYMAX.LE.0.0) THEN
       GYORG=0.0
    ELSE
       GYORG=(INT(GYMIN/(2*GYSCL))+1)*2*GYSCL
    ENDIF

    CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
    GZLMAX=LOG10(GZMAX)
    GZLMIN=LOG10(GZMIN)
    IF(GZLMAX.LT.0.0) THEN
       IGMAX=-INT(-GZLMAX)-1
    ELSE
       IGMAX=INT(GZLMAX)
    ENDIF
    IF(GZLMIN.LT.0.0) THEN
       IGMIN=-INT(-GZLMIN)-1
    ELSE
       IGMIN=INT(GZLMIN)
    ENDIF
    IF(IGMAX.GT.2) THEN
       IF(IGMIN.GT.-2) THEN
          IGMAX=IGMIN+4
       ELSE
          IGMAX=2
       ENDIF
    ENDIF

    CALL PAGES
    CALL SETFNT(32)
    CALL SETCHS(0.3,0.0)
    CALL SETLIN(0,0,7)
    CALL SETLNW(0.035)
    CALL GDEFIN( 4.0,24.0, 3.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    CALL GSCALE(GXORG,  GXSCL,0.0,0.0,0.2,9)
    CALL GVALUE(GXORG,2*GXSCL,0.0,0.0,NGULEN(2*GXSCL))
    CALL GSCALE(0.0,0.0,GYORG,  GYSCL,0.2,9)
    CALL GVALUE(0.0,0.0,GYORG,2*GYSCL,NGULEN(2*GYSCL))
    CALL GSCALE(0.0,0.0,0.0,GYMAX-GYMIN,0.0,0)
    CALL SETFNT(0)

    IPRD=0
    GZSTEP=GZMAX-GZMIN

!      WRITE(6,*) '2D ',GZMIN,GZMAX,IGMAX,10.0**IGMAX,10.0**(IGMAX-4)

    DO I=1,5
       CALL SETRGB(RGBS(1,I),RGBS(2,I),RGBS(3,I))
       GL=10.0**(IGMAX-I+1)
       IF(10*GL.GT.GZMIN.AND.GL.LT.GZMAX) THEN
          GZORG=GL
          CALL SETLNW(0.035)
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,0,KA)
          GZORG=1.5*GL
          CALL SETLNW(0.018)
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,2,KA)
          GZORG=2.0*GL
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,2,KA)
          GZORG=3.0*GL
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,2,KA)
          GZORG=5.0*GL
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,2,KA)
          GZORG=7.0*GL
          CALL CONTQ1(GZ,NGXMAX,NGXMAX,NGYMAX,GZORG,GZSTEP,1,IPRD,2,KA)
       ENDIF
    ENDDO

    CALL SETRGB(0.0,1.0,0.0)
    CALL SETVEW( 4.0,24.0, 3.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
    DO NINF=1,NINFMAX
       IF(GFINF(3,NINF).LT.1.E30) THEN
          CALL MOVE(GFINF(1,NINF),GFINF(2,NINF))
          CALL NUMBI(NINF,'(I1)',1)
       ENDIF
    ENDDO
    CALL OFFVEW

    CALL SETFNT(32)
    CALL SETLIN(0,0,7)
    CALL MOVE(12.0,0.5)
    CALL TEXT('Re(f) [MHz]',11)
    CALL MOVE(4.0,17.0)
    CALL TEXT('AMP2D',5)
    CALL MOVE(7.0,17.0)
    CALL TEXT('   min=',7)
    CALL NUMBR(GZMIN,'(1PE12.4)',12)
    CALL TEXT('   max=',7)
    CALL NUMBR(GZMAX,'(1PE12.4)',12)

    CALL PAGEE
    RETURN
  END SUBROUTINE WMEG2D

!     ***** FIND EIGEN VALUE *****

  SUBROUTINE wm_eign(KID,LINE)

    USE wmcomm
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE wmemfp,ONLY: wm_bfield
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80),INTENT(OUT):: LINE
    CHARACTER(LEN=1),INTENT(OUT):: KID
    INTEGER:: MODE,IERR
    REAL(rkind):: X,Y,XX,YY

    MODE=0
1   CONTINUE
    IF(NRANK.EQ.0) THEN
2      CONTINUE
       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4)') &
            ' ## FRINI,FIINI=',FRINI,FIINI
       WRITE(6,'(A)') ' ## INPUT : KID or FRINI,FINII'
       CALL TASK_KLIN(LINE,KID,MODE,wm_parm)
       IF(MODE.EQ.0) READ(LINE,*,ERR=2,END=2) FRINI,FIINI
    ENDIF
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE.EQ.2) THEN
       CALL wm_broadcast
       GOTO 1
    ENDIF
    IF(MODE.EQ.3) GOTO 1

    CALL mtx_broadcast1_real8(FRINI)
    CALL mtx_broadcast1_real8(FIINI)

    X=FRINI
    Y=FIINI
    IF (MODENW.EQ.0) THEN
       CALL NEWTN0(DIAMIN,X,Y,XX,YY,DLTNW,EPSNW,LMAXNW,LISTNW,NRANK,IERR)
    ELSEIF(MODENW.EQ.1) THEN
       CALL NEWTN1(DIAMIN,X,Y,XX,YY,DLTNW,EPSNW,LMAXNW,LISTNW,NRANK,IERR)
    ELSE
       XX=X
       YY=Y
    ENDIF

    FRINI=XX
    FIINI=YY
    CALL wm_bfield

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_eign

!     ***** PARAMETER SCAN *****

  SUBROUTINE wm_scan(KID,LINE)

    USE wmcomm
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE wmemfp,ONLY: wm_bfield
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80),INTENT(OUT):: LINE
    CHARACTER(LEN=1),INTENT(OUT):: KID

    REAL(rkind),ALLOCATABLE:: &
         PN_SAVE(:),PTPR_SAVE(:),PTPP_SAVE(:),PU_SAVE(:), &
         PNITB_SAVE(:),PTITB_SAVE(:),PUITB_SAVE(:),RHOITB_SAVE(:)
    REAL,ALLOCATABLE:: GX(:),GY(:,:)
    INTEGER:: MODE,ISCAN,NPH0_SAVE,MODE1,NSC,NS,IERR,I
    CHARACTER(LEN=4):: KV
    REAL(rkind):: DSC,X,Y,SC,XX,YY,AMPL
    REAL(rkind):: RR_SAVE,QA_SAVE,BB_SAVE,PNA_SAVE,PNAL_SAVE
    REAL(rkind):: QMIN_SAVE,RHOMIN_SAVE
    EXTERNAL:: GUFLSH,GMNMX1,GQSCAL,PAGES,SETCHS,SETLIN,PAGEE
    EXTERNAL:: GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,MOVE,TEXT

    ALLOCATE(PN_SAVE(nsmax),PTPR_SAVE(nsmax),PTPP_SAVE(nsmax),PU_SAVE(nsmax))
    ALLOCATE(PNITB_SAVE(nsmax),PTITB_SAVE(nsmax),PUITB_SAVE(nsmax))
    ALLOCATE(RHOITB_SAVE(nsmax))
    
    MODE=0
    ISCAN=0
1   CONTINUE
    IF(NRANK.EQ.0) THEN
2      CONTINUE
       IF(MODE.EQ.0) THEN
          WRITE(6,'(A)') &
               ' ## ISCAN: 1:NPH0 2:RR 3:QA 4:BB 5:PNA 6:PNAL 7:PT 8:PN'
          WRITE(6,'(A)') &
               ' ##        9:qmin 10:rhomin 11:itb 12:rhoitb 13:PU'
       ENDIF
       WRITE(6,'(A)') ' ## INPUT : KID or ISCAN'
       CALL TASK_KLIN(LINE,KID,MODE,wm_parm)
       IF(MODE.EQ.0) &
            READ(LINE,*,ERR=2,END=2) ISCAN
    ENDIF
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE.EQ.2) THEN
       CALL wm_broadcast
       GOTO 1
    ENDIF
    IF(MODE.EQ.3) GOTO 1

    CALL mtx_broadcast1_integer(ISCAN)
    IF(ISCAN.EQ.0.OR.ISCAN.EQ.9) GOTO 1

    SELECT CASE(ISCAN)
    CASE(1)
       KV='NPH0'
       NPH0_SAVE=NPH0
       SCMIN=DBLE(NPH0)
    CASE(2)
       KV='RR  '
       RR_SAVE=RR
       SCMIN=RR
    CASE(3)
       KV='QA  '
       QA_SAVE=QA
       SCMIN=QA
    CASE(4)
       KV='BB  '
       BB_SAVE=BB
       SCMIN=BB
    CASE(5)
       KV='PNA '
       PNA_SAVE=PNA
       SCMIN=PNA
    CASE(6)
       KV='PNAL'
       PNAL_SAVE=PNAL
       SCMIN=PNAL
    CASE(7)
       KV='PT  '
       DO NS=1,NSMAX
          PTPR_SAVE(NS)=PTPR(NS)
          PTPP_SAVE(NS)=PTPP(NS)
       ENDDO
       SCMIN=PTPR(1)
    CASE(8)
       KV='PN  '
       DO NS=1,NSMAX
          PN_SAVE(NS)=PN(NS)
       ENDDO
       SCMIN=PN(1)
    CASE(9)
       KV='QMIN'
       QMIN_SAVE=QMIN
       SCMIN=QMIN
    CASE(10)
       KV='RHOM'
       RHOMIN_SAVE=RHOMIN
       SCMIN=RHOMIN
    CASE(11)
       KV='PNTB'
       DO NS=1,NSMAX
          PNITB_SAVE(NS)=PNITB(NS)
          PTITB_SAVE(NS)=PTITB(NS)
          PUITB_SAVE(NS)=PUITB(NS)
       ENDDO
       SCMIN=0.D0
       SCMAX=1.D0
    CASE(12)
       KV='RHOB'
       DO NS=1,NSMAX
          RHOITB_SAVE(NS)=RHOITB(NS)
       END DO
       SCMIN=0.D0
       SCMAX=1.D0
    CASE(13)
       KV='PU  '
       DO NS=1,NSMAX
          PU_SAVE(NS)=PU(NS)
       ENDDO
       SCMIN=PU(1)
    END SELECT

    MODE1=0
11  CONTINUE
    IF(NRANK.EQ.0) THEN
12     CONTINUE
       IF(MODE1.EQ.0) WRITE(6,'(A/1P2E12.4,I5,1P2E12.4)') &
            ' ## SCMIN,SCMAX,NSCMAX,FRINI,FIINI', &
                 SCMIN,SCMAX,NSCMAX,FRINI,FIINI
       WRITE(6,'(A)') &
            ' ## INPUT : KID or SCMIN,SCMAX,NSCMAX,FRINI,FIINI'
       CALL TASK_KLIN(LINE,KID,MODE1,wm_parm)
       IF(MODE1.EQ.0) &
            READ(LINE,*,ERR=12,END=12) SCMIN,SCMAX,NSCMAX,FRINI,FIINI
    ENDIF
    CALL mtx_broadcast1_integer(MODE1)
    IF(MODE1.EQ.1) THEN
       CALL mtx_broadcast1_character(KID)
       GOTO 9000
    ENDIF
    IF(MODE1.EQ.2) THEN
       CALL wm_broadcast
       GOTO 11
    ENDIF
    IF(MODE1.EQ.3) GOTO 11

    CALL mtx_broadcast1_real8(SCMIN)
    CALL mtx_broadcast1_real8(SCMAX)
    CALL mtx_broadcast1_integer(NSCMAX)
    CALL mtx_broadcast1_real8(FRINI)
    CALL mtx_broadcast1_real8(FIINI)

    IF(ALLOCATED(GX)) DEALLOCATE(GX,GY)
    ALLOCATE(GX(NSCMAX),GY(NSCMAX,3))

    DSC=(SCMAX-SCMIN)/(NSCMAX-1)
    X=FRINI
    Y=FIINI

    DO NSC=1,NSCMAX
       SC=SCMIN+DSC*(NSC-1)
       GX(NSC)=GUCLIP(SC)
    ENDDO

    DO NSC=1,NSCMAX

       SC=SCMIN+DSC*(NSC-1)
       SELECT CASE(ISCAN)
       CASE(1)
          NPH0=NINT(SC)
       CASE(2)
          RR=SC
       CASE(3)
          QA=SC
       CASE(4)
          BB=SC
       CASE(5)
          PNA=SC
       CASE(6)
          PNAL=SC
       CASE(7)
          DO NS=1,NSMAX
             PTPR(NS)=SC*PTPR_SAVE(NS)/PTPR_SAVE(1)
             PTPP(NS)=SC*PTPP_SAVE(NS)/PTPR_SAVE(1)
          ENDDO
       CASE(8)
          DO NS=1,NSMAX
             PN(NS)=PN_SAVE(NS)*SC/PN_SAVE(1)
          ENDDO
       CASE(9)
          QMIN=SC
       CASE(10)
          RHOMIN=SC
       CASE(11)
          DO NS=1,NSMAX
             PNITB(NS)=PNITB_SAVE(NS)*SC
             PTITB(NS)=PTITB_SAVE(NS)*SC
             PUITB(NS)=PUITB_SAVE(NS)*SC
          ENDDO
       CASE(12)
          DO NS=1,NSMAX
             RHOITB(NS)=RHOITB_SAVE(NS)*SC
          END DO
       CASE(13)
          IF(PU_SAVE(1).EQ.0.D0) THEN
             DO NS=1,NSMAX
                PU(NS)=SC
             ENDDO
          ELSE
             DO NS=1,NSMAX
                PU(NS)=PU_SAVE(NS)*SC/PU_SAVE(1)
             ENDDO
          ENDIF
       END SELECT

       IF(MODENW.EQ.0) THEN
          CALL NEWTN0(DIAMIN,X,Y,XX,YY, &
                      DLTNW,EPSNW,LMAXNW,LISTNW,NRANK,IERR)
       ELSEIF(MODENW.EQ.1) THEN
          CALL NEWTN1(DIAMIN,X,Y,XX,YY, &
                      DLTNW,EPSNW,LMAXNW,LISTNW,NRANK,IERR)
       ELSE
          XX=X
          YY=Y
       ENDIF
       IF(IERR.NE.0) GOTO 3000

       CALL DIAMIN(XX,YY,AMPL)
       GY(NSC,1)=GUCLIP(XX)
       GY(NSC,2)=GUCLIP(YY)
       GY(NSC,3)=GUCLIP(AMPL)
       X=XX
       Y=YY
    ENDDO
3000 CONTINUE

    IF(IERR.EQ.0) THEN
       FRINI=XX
       FIINI=YY
       CALL wm_bfield
    ENDIF

    SELECT CASE(ISCAN)
    CASE(1)
       NPH0=NPH0_SAVE
    CASE(2)
       RR=RR_SAVE
    CASE(3)
       QA=QA_SAVE
    CASE(4)
       BB=BB_SAVE
    CASE(5)
       PNA=PNA_SAVE
    CASE(6)
       PNAL=PNAL_SAVE
    CASE(7)
       DO NS=1,NSMAX
          PTPR(NS)=PTPR_SAVE(NS)
          PTPP(NS)=PTPP_SAVE(NS)
       ENDDO
    CASE(8)
       DO NS=1,NSMAX
          PN(NS)=PN_SAVE(NS)
       ENDDO
    CASE(9)
       QMIN=QMIN_SAVE
    CASE(10)
       RHOMIN=RHOMIN_SAVE
    CASE(11)
       DO NS=1,NSMAX
          PNITB(NS)=PNITB_SAVE(NS)
          PTITB(NS)=PTITB_SAVE(NS)
          PUITB(NS)=PUITB_SAVE(NS)
       ENDDO
    CASE(12)
       DO NS=1,NSMAX
          RHOITB(NS)=RHOITB_SAVE(NS)
       END DO
    CASE(13)
       DO NS=1,NSMAX
          PU(NS)=PU_SAVE(NS)
       ENDDO
    END SELECT

    CALL WMSC1G(GX,GY(1,1),NSCMAX,ISCAN,KV)
    CALL WMSC1G(GX,GY(1,2),NSCMAX,ISCAN,KV)

    IF(NRANK.EQ.0) THEN
       WRITE(6,*) '      '//KV
       WRITE(6,'(I5,1P4E12.4)') &
            (I,GX(I),GY(I,1),GY(I,2),GY(I,3),I=1,NSCMAX) 
       CALL GUFLSH
    ENDIF

    GOTO 11

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_scan

! --- wm_scan internal routine

  SUBROUTINE WMSC1G(GX,GY,NSCMAX,ISCAN,KV)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSCMAX,ISCAN
    REAL,INTENT(IN):: GX(NSCMAX),GY(NSCMAX)
    CHARACTER(LEN=4),INTENT(IN)::  KV
    INTEGER:: NKETA
    REAL:: GXMIN1,GXMAX1,GXMIN,GXMAX,GXSCL
    REAL:: GYMIN1,GYMAX1,GYMIN,GYMAX,GYSCL
    INTERFACE
       FUNCTION NGULEN(X)
         USE bpsd_kinds
         REAL:: X
         INTEGER:: NGULEN
       END FUNCTION NGULEN
    END INTERFACE
    EXTERNAL:: GMNMX1,GQSCAL,PAGES,SETCHS,SETLIN,GDEFIN,GFRAME,GSCALE,GVALUE
    EXTERNAL:: GPLOTP,MOVE,TEXT,PAGEE

    SELECT CASE(ISCAN)
    CASE(1)
       NKETA=1
    CASE(2)
       NKETA=3
    CASE(5,6)
       NKETA=4
    CASE DEFAULT
       NKETA=2
    END SELECT

!     ****MAP****

    CALL GMNMX1(GX,1,NSCMAX,1,GXMIN1,GXMAX1)
    CALL GMNMX1(GY,1,NSCMAX,1,GYMIN1,GYMAX1)
    CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSCL)
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSCL)

    CALL PAGES
    CALL SETCHS(0.35,0.0)
    CALL SETLIN(0,0,4)
    CALL GDEFIN(4.5,24.5,3.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,  GXSCL,0.0,0.0,0.2, 9)
    CALL GVALUE(GXMIN,2*GXSCL,0.0,0.0,NKETA)  
    CALL GSCALE(0.0,0.0,GYMIN,  GYSCL,0.2, 9)
    CALL GVALUE(0.0,0.0,GYMIN,2*GYSCL,NGULEN(2*GYSCL))  
    CALL GPLOTP(GX,GY,1,NSCMAX,1,-1,1,0) 

    CALL SETCHS(0.3,0.0)
    CALL MOVE(4.0,17.0)
    CALL TEXT(KV,4)

    CALL PAGEE
    RETURN

  END SUBROUTINE WMSC1G

!     ***** CALCULATE Inverse of wave amplitude to find its minimum *****

  SUBROUTINE DIAMIN(X,Y,F)

    USE wmcomm
    USE wmsetg
    USE wmsolv
    USE wmemfp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y
    REAL(rkind),INTENT(OUT):: F
    INTEGER:: ierr,NR,NHH,NTH,ND,NDX,MD,MDX
    REAL(rkind):: ESUM,EABS_MAX,EABS,ACE_MAX,EABS2
    COMPLEX(rkind):: CE_MAX,CFACT
    
    RF=X
    RFI=Y

    MODEEG=1

    CALL wm_setg(ierr)
    CALL wm_solv
    CALL wm_efield

    ESUM=0.D0
    EABS_MAX=0.D0
    CE_MAX=(0.D0,0.D0)
    DO NR=1,NRMAX+1
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             EABS=ABS(CEFLD(1,NTH,NHH,NR))**2 &
                 +ABS(CEFLD(2,NTH,NHH,NR))**2 &
                 +ABS(CEFLD(3,NTH,NHH,NR))**2
             ESUM=ESUM+EABS
             EABS=ABS(CEFLD(2,NTH,NHH,NR))**2
             IF(EABS.GT.EABS_MAX) THEN
                EABS_MAX=EABS
                CE_MAX=CEFLD(2,NTH,NHH,NR)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF(ABS(CE_MAX).EQ.0.D0) THEN
       CFACT=1.D0
    ELSE
       CFACT=1.D0/CE_MAX
    ENDIF
    DO NR=1,NRMAX+1
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             CEFLD(1,NTH,NHH,NR)=CFACT*CEFLD(1,NTH,NHH,NR)
             CEFLD(2,NTH,NHH,NR)=CFACT*CEFLD(2,NTH,NHH,NR)
             CEFLD(3,NTH,NHH,NR)=CFACT*CEFLD(3,NTH,NHH,NR)
          ENDDO
       ENDDO
    ENDDO

    EABS_MAX=0.D0
    CE_MAX=0.D0
    ACE_MAX=0.D0
    DO NR=1,NRMAX+1
       EABS2=0.D0
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             EABS2=EABS2+DREAL(CEFLDK(2,MDX,NDX,NR))**2 &
                        +DIMAG(CEFLDK(2,MDX,NDX,NR))**2
             IF(ABS(CEFLDK(2,MDX,NDX,NR)).GT.ACE_MAX) THEN
                CE_MAX=CEFLDK(2,MDX,NDX,NR)
                ACE_MAX=ABS(CE_MAX)
             ENDIF
          ENDDO
       ENDDO
       IF(EABS2.GT.EABS_MAX) EABS_MAX=EABS2
    ENDDO

    IF(EABS_MAX.EQ.0.D0) THEN
       CFACT=1.D0
    ELSE
       CFACT=ABS(CE_MAX)/(CE_MAX*SQRT(EABS_MAX))
    ENDIF
    DO NR=1,NRMAX+1
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             CEFLDK(1,MDX,NDX,NR)=CFACT*CEFLDK(1,MDX,NDX,NR)
             CEFLDK(2,MDX,NDX,NR)=CFACT*CEFLDK(2,MDX,NDX,NR)
             CEFLDK(3,MDX,NDX,NR)=CFACT*CEFLDK(3,MDX,NDX,NR)
          ENDDO
       ENDDO
    ENDDO

    !      WRITE(6,*) NRMAX,NTHMAX,NHHMAX,ESUM,RF
    
    F=NRMAX*NTHMAX*NHHMAX/(ESUM*RF**4)
    IF(F.LE.1.D-15) F=1.D-15
    AMP_EIGEN=F

    RETURN
  END SUBROUTINE DIAMIN

!     ****** TWO-DIMENSIONAL NEWTON METHOD ******

  SUBROUTINE NEWTN0(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,NRANK,IERR)

    USE wmcomm,ONLY: rkind
    USE libmpi
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: X,Y
    REAL(rkind),INTENT(IN):: DELT,EPS
    INTEGER,INTENT(IN):: ILMAX,LIST,NRANK
    REAL(rkind),INTENT(OUT):: XX,YY
    INTEGER,INTENT(OUT):: IERR
    LOGICAL:: LDUMP,LPRINT
    INTEGER:: ITER
    REAL(rkind):: HX,HY,FX,FY,FXX,FXY,FYY,DF,DET,FXN,FYN,DFN
    REAL(rkind):: F00,FM0,FP0,F0M,F0P,FMM,FMP,FPM,FPP
    REAL(rkind):: H11,H12,H21,H22,DX,DY,V,DV,TT
    INTERFACE
       SUBROUTINE SUB(X,Y,F)
         USE wmcomm,ONLY: rkind
         REAL(rkind),INTENT(IN):: X,Y
         REAL(rkind),INTENT(OUT):: F
       END SUBROUTINE SUB
    END INTERFACE
    EXTERNAL:: GUFLSH
    
    LDUMP=NRANK.EQ.0.AND.LIST.GE.2
    LPRINT=NRANK.EQ.0.AND.LIST.EQ.1

    IERR=0
    ITER=0
    IF(ABS(X).GT.1.D0) THEN
       HX=DELT*X
    ELSE
       HX=DELT
    ENDIF
    IF(ABS(Y).GT.1.D0) THEN
       HY=DELT*Y
    ELSE
       HY=DELT
    ENDIF
    IF(LDUMP) WRITE(6,'(A,1P2E12.4)') 'X,Y   =',X,Y
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X,   Y,   F00)
    IF(LPRINT) WRITE(6,600) X,Y,0.D0,0.D0,F00
    IF(LPRINT) CALL GUFLSH
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y,F00
    IF(LDUMP) CALL GUFLSH

    IF(F00.LE.2.D-15) GOTO 9001

    CALL SUB(X-HX,Y,   FM0)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X-HX,Y,FM0
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X+HX,Y,   FP0)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X+HX,Y,FP0
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X,   Y-HY,F0M)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y-HY,F0M
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X,   Y+HY,F0P)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y+HY,F0P
    IF(LDUMP) CALL GUFLSH

    FX=(FP0-FM0)/(2*HX)
    FY=(F0P-F0M)/(2*HY)

1   CONTINUE
    CALL SUB(X-HX,Y-HY,FMM)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X-HX,Y-HY,FMM
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X-HX,Y+HY,FMP)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X-HX,Y+HY,FMP
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X+HX,Y-HY,FPM)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X+HX,Y-HY,FPM
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X+HX,Y+HY,FPP)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X+HX,Y+HY,FPP
    IF(LDUMP) CALL GUFLSH

    FXX=(FP0-2*F00+FM0)/(HX*HX)
    FYY=(F0P-2*F00+F0M)/(HY*HY)
    FXY=(FPP-FPM-FMP+FMM)/(4*HX*HY)
    IF(LDUMP) WRITE(6,'(A,1P2E12.4)') 'FX,FY =',FX,FY
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'FXX,FXY,FXY =',FXX,FXY,FYY
    IF(LDUMP) CALL GUFLSH

    DF=SQRT(FX*FX+FY*FY)
    DET=FXX*FYY-FXY*FXY
    IF(LDUMP) WRITE(6,'(A,1P2E12.4)') 'DF,DET:0 =',DF,DET
    IF(LDUMP) CALL GUFLSH

    IF(DF.LE.2.D-15) GO TO 9000
    CALL mtx_barrier

    H11= FYY/DET
    H12=-FXY/DET
    H21=-FXY/DET
    H22= FXX/DET
    DX=-(H11*FX+H12*FY)
    DY=-(H21*FX+H22*FY)
    V=SQRT(X*X+Y*Y+EPS)
    DV=SQRT(DX*DX+DY*DY)
    IF(DV/V.LE.EPS) GO TO 9000

    TT=1.D0
2   CONTINUE
    X=X+TT*DX
    Y=Y+TT*DY

    CALL SUB(X,   Y,   F00)
    IF(LPRINT) WRITE(6,600) X,Y,TT*DX,TT*DY,F00
    IF(LPRINT) CALL GUFLSH
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y,F00
    IF(LDUMP) CALL GUFLSH

    IF(F00.LE.2.D-15) GOTO 9001
    CALL mtx_barrier

    CALL SUB(X-HX,Y,   FM0)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X-HX,Y,FM0
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X+HX,Y,   FP0)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X+HX,Y,FP0
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X,   Y-HY,F0M)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y-HY,F0M
    IF(LDUMP) CALL GUFLSH

    CALL SUB(X,   Y+HY,F0P)
    IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y+HY,F0P
    IF(LDUMP) CALL GUFLSH

    FXN=(FP0-FM0)/(2*HX)
    FYN=(F0P-F0M)/(2*HY)

    DFN=SQRT(FXN*FXN+FYN*FYN)
    IF(DFN.GT.DF) THEN
       X=X-TT*DX
       Y=Y-TT*DY
       TT=0.5D0*TT
       ITER=ITER+1
       IF(TT.LE.1.D-3) GOTO 8000
       IF(ITER.LE.ILMAX) GOTO 2
    ELSE
       FX=FXN
       FY=FYN
       DF=DFN
       ITER=ITER+1
       IF(ITER.LE.ILMAX) GO TO 1
    ENDIF

    IERR=2
    IF(LPRINT)  WRITE(6,*) 'XX NEWTN0: LOOP COUNT EXCEEDS UPPER BOUND.'
    GOTO 9001

8000 CONTINUE
    IERR=1
    IF(LPRINT) WRITE(6,*) 'XX NEWTN0: DOES NOT CONVERGE.'
    GOTO 9001

9000 CONTINUE
    CALL SUB(X,Y,F00)
9001 CONTINUE
    XX=X
    YY=Y
    RETURN
600 FORMAT(' ',3X,'X,Y,DX,DY,F = ',1P2E15.7,1P3E10.2)
  END SUBROUTINE NEWTN0

!     ****** TWO-DIMENSIONAL NEWTON METHOD (SIMPLER) ******

  SUBROUTINE NEWTN1(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,NRANK,IERR)

    USE wmcomm,ONLY: rkind
    USE libmpi
    IMPLICIT NONE
    INTERFACE
       SUBROUTINE SUB(X,Y,F)
         USE wmcomm,ONLY: rkind
         REAL(rkind),INTENT(IN):: X,Y
         REAL(rkind),INTENT(OUT):: F
       END SUBROUTINE SUB
    END INTERFACE
    REAL(rkind),INTENT(INOUT):: X,Y
    REAL(rkind),INTENT(IN):: DELT,EPS
    INTEGER,INTENT(IN):: ILMAX,LIST,NRANK
    REAL(rkind),INTENT(OUT):: XX,YY
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ITER
    REAL(rkind):: HX,HY,F00,FM0,FP0,F0M,F0P,FX,FY,FPP
    REAL(rkind):: FXX,FYY,FXY,DF,DET,H11,H12,H21,H22,DX,DY,TT
    REAL(rkind):: FXN,FYN,DFN,V,DV
    EXTERNAL:: GUFLSH

    IERR=0
    ITER=0
    IF(ABS(X).GT.1.D0) THEN
       HX=DELT*X
    ELSE
       HX=DELT
    ENDIF
    IF(ABS(Y).GT.1.D0) THEN
       HY=DELT*Y
    ELSE
       HY=DELT
    ENDIF
    CALL SUB(X,   Y,   F00)
    CALL SUB(X-HX,Y,   FM0)
    CALL SUB(X+HX,Y,   FP0)
    CALL SUB(X,   Y-HY,F0M)
    CALL SUB(X,   Y+HY,F0P)
    FX=(FP0-FM0)/(2*HX)
    FY=(F0P-F0M)/(2*HY)

1   CONTINUE
    CALL SUB(X+HX,Y+HY,FPP)

    FXX=(FP0-2*F00+FM0)/(HX*HX)
    FYY=(F0P-2*F00+F0M)/(HY*HY)
    FXY=(FPP-FP0-F0P+F00)/(HX*HY)
    DF=SQRT(FX*FX+FY*FY)
    DET=FXX*FYY-FXY*FXY
    H11= FYY/DET
    H12=-FXY/DET
    H21=-FXY/DET
    H22= FXX/DET

    DX=-(H11*FX+H12*FY)
    DY=-(H21*FX+H22*FY)
    TT=1.D0
2   CONTINUE
    X=X+TT*DX
    Y=Y+TT*DY

    CALL SUB(X,   Y,   F00)
    CALL SUB(X-HX,Y,   FM0)
    CALL SUB(X+HX,Y,   FP0)
    CALL SUB(X,   Y-HY,F0M)
    CALL SUB(X,   Y+HY,F0P)
    FXN=(FP0-FM0)/(2*HX)
    FYN=(F0P-F0M)/(2*HY)
    IF(NRANK.EQ.0) THEN
       IF(LIST.GT.0) WRITE(6,600) X,Y,DX,DY,F00
       IF(LIST.GT.0) CALL GUFLSH
    ENDIF
    DFN=SQRT(FXN*FXN+FYN*FYN)
    IF(DFN.GT.DF) THEN
       X=X-TT*DX
       Y=Y-TT*DY
       TT=0.5D0*TT
       ITER=ITER+1
       IF(TT.LE.1.D-3) GOTO 8000
       IF(ITER.LE.ILMAX) GOTO 2
    ELSE
       FX=FXN
       FY=FYN
       DF=DFN
       ITER=ITER+1
       V=SQRT(X*X+Y*Y+EPS)
       DV=SQRT(DX*DX+DY*DY)
       IF(DV/V.LE.EPS) GO TO 9000
       IF(ITER.LE.ILMAX) GO TO 1
    ENDIF

    IERR=2
    IF(NRANK.EQ.0) THEN
       IF(LIST.GT.0) WRITE(6,*) 'XX NEWTN1: LOOP COUNT EXCEEDS UPPER BOUND.'
    ENDIF
    GOTO 9000

8000 CONTINUE
    IERR=1
    IF(NRANK.EQ.0) THEN
       IF(LIST.GT.0) WRITE(6,*) 'XX NEWTN1: DOES NOT CONVERGE.'
    ENDIF
    GOTO 9000

9000 CONTINUE
    XX=X
    YY=Y
    RETURN
600 FORMAT(' ',3X,'X,Y,DX,DY,F = ',1P2E15.7,1P3E10.2)
  END SUBROUTINE NEWTN1

!     ***** WRITE WAVE DATA *****

  SUBROUTINE wm_wout

    USE wmcomm
    IMPLICIT NONE
    INTEGER:: NF,NR,ND,NDX,MD,MDX
    
    NF=26
    WRITE(NF,*) 'FR [MHz] = ',FRINI
    WRITE(NF,*) 'FI [MHz] = ',FIINI
    WRITE(NF,*) 'AMP_EIGEN= ',AMP_EIGEN
    WRITE(NF,*) 'NR COUNT = ',NRMAX+1
    WRITE(NF,*) 'MD COUNT = ',MDMAX-MDMIN+1
    WRITE(NF,*) 'ND COUNT = ',NDMAX-NDMIN+1
    WRITE(NF,*)

    WRITE(NF,*) 'r [m]:'
    WRITE(NF,'(1P5E14.6)') (XR(NR),NR=1,NRMAX+1)
    WRITE(NF,*)

    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1
          WRITE(NF,*) 'MD           = ',NTH0+MD
          WRITE(NF,*) 'ND           = ',NPH0+NHC*ND
          WRITE(NF,*) 'Re Er:'
          WRITE(NF,'(1P5E14.6)') &
               (DREAL(CEFLDK(1,MDX,NDX,NR)),NR=1,NRMAX+1)
          WRITE(NF,*) 'Im Er:'
          WRITE(NF,'(1P5E14.6)') &
               (DIMAG(CEFLDK(1,MDX,NDX,NR)),NR=1,NRMAX+1)
          WRITE(NF,*) 'Re Etheta:'
          WRITE(NF,'(1P5E14.6)') &
               (DREAL(CEFLDK(2,MDX,NDX,NR)),NR=1,NRMAX+1)
          WRITE(NF,*) 'Im Etheta:'
          WRITE(NF,'(1P5E14.6)') &
               (DIMAG(CEFLDK(2,MDX,NDX,NR)),NR=1,NRMAX+1)
          WRITE(NF,*) 'Re Ephi:'
          WRITE(NF,'(1P5E14.6)') &
               (DREAL(CEFLDK(3,MDX,NDX,NR)),NR=1,NRMAX+1)
          WRITE(NF,*) 'Im Ephi:'
          WRITE(NF,'(1P5E14.6)') &
               (DIMAG(CEFLDK(3,MDX,NDX,NR)),NR=1,NRMAX+1)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE wm_wout
END MODULE wmeign
