C     $Id$
C
C     ***** SINGLE CALCULATION *****
C
      SUBROUTINE WMAM0D(KID,LINE)
C
      INCLUDE 'wmcomm.inc'
      EXTERNAL WMPARM
      CHARACTER LINE*80,KID*1
C
      MODE=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
    2       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4)')
     &      ' ## FR,FI=',FRINI,FIINI
            WRITE(6,'(A)') ' ## INPUT : KID or FR,FI'
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
            IF(MODE.EQ.0)
     &      READ(LINE,*,ERR=2,END=2) FRINI,FIINI
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
         CALL MPBCDA(FRINI)
         CALL MPBCDA(FIINI)
C
         CALL DIAMIN(FRINI,FIINI,AMPL)
         
         CALL WMBFLD

      GOTO 1
C
 9000 CONTINUE
      RETURN
      END
C
C     ***** 1D PLOT OF AMPLITUDE *****
C
      SUBROUTINE WMAM1D(KID,LINE)
C
      INCLUDE 'wmcomm.inc'
      DIMENSION GX(NGZM),GZ(NGZM,1)
      CHARACTER LINE*80,KID*1
      EXTERNAL WMPARM
C
      MODE=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
    2       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4,I5,1PE12.4)')
     &      ' ## FRMIN,FRMAX,NGFMAX,FI0=',FRMIN,FRMAX,NGFMAX,FI0
            WRITE(6,'(A)') ' ## INPUT : KID or FRMIN,FRMAX,NGFMAX,FI0'
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
            IF(MODE.EQ.0)
     &      READ(LINE,*,ERR=2,END=2) FRMIN,FRMAX,NGFMAX,FI0
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
         CALL MPBCIA(NGFMAX)
         CALL MPBCDA(FRMIN)
         CALL MPBCDA(FRMAX)
         CALL MPBCDA(FI0)
C
         DELTFR=(FRMAX-FRMIN)/(NGFMAX-1)
         AMPMIN=1.D30
         DO NGF=1,NGFMAX
            FR=FRMIN+(NGF-1)*DELTFR
            FI=FI0
            CALL DIAMIN(FR,FI,AMPL)
            IF(MYRANK.EQ.0) THEN
               IF(LISTEG.GE.1) THEN
                  WRITE(6,'(A,1P3E12.4)') 
     &              '      FR,FI,AMPL = ',FR,FI,AMPL
                  CALl GUFLSH
               ENDIF
            ENDIF
            GX(NGF)=GUCLIP(FR)
            GZ(NGF,1)=GUCLIP(LOG10(AMPL))
            IF(AMPL.LE.AMPMIN) THEN
               FRINI=FR
               FIINI=FI
               AMPMIN=AMPL
            ENDIF
         ENDDO
         IF(MYRANK.EQ.0) THEN
            WRITE(6,'(A,1P3E12.4)') '      FR,FI,AMPMIN=',
     &                                     FRINI,FIINI,AMPMIN
            CALL GUFLSH
         ENDIF
C
         CALL WMEG1D(GX,GZ,NGZM,NGFMAX,1,GUCLIP(FI0))
C
      GOTO 1
C
 9000 CONTINUE
      RETURN
      END
C
      SUBROUTINE WMEG1D(GX,GZ,NGZM,NGXMAX,NGYMAX,GFI0)
C
      DIMENSION GX(NGZM),GZ(NGZM,NGYMAX)
C
      GXMIN1=GX(1)
      GXMAX1=GX(NGXMAX)
      CALL GMNMX2(GZ,NGZM,1,NGXMAX,1,1,NGYMAX,1,GZMIN1,GZMAX1)
      IF(GZMIN1.LT.-10.0) GZMIN1=-10.0
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSCL)
      CALL GQSCAL(GZMIN1,GZMAX1,GZMIN,GZMAX,GZSCL)
C
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
C
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
C
      CALL PAGEE
      RETURN
      END
C
C     ***** 2D PLOT OF AMPLITUDE *****
C
      SUBROUTINE WMAM2D(KID,LINE)
C
      INCLUDE 'wmcomm.inc'
      DIMENSION GX(NGZM),GY(NGZM)
      DIMENSION GZ(NGZM,NGZM)
      PARAMETER (NINFM=5)
      DIMENSION GFINF(3,NINFM),KA(8,NGZM,NGZM)
      CHARACTER LINE*80,KID*1
      EXTERNAL WMPARM
C
      MODE=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
    2       IF(MODE.EQ.0) THEN
               WRITE(6,'(A)') 
     &            ' ## FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX='
               WRITE(6,'(3X,1P2E12.4,I5,1P2E12.4,I5)')
     &                 FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX
            ENDIF
            WRITE(6,'(A,A)') ' ## INPUT : KID or ',
     &               'FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX'
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
            IF(MODE.EQ.0)
     &      READ(LINE,*,ERR=2,END=2)
     &                         FRMIN,FRMAX,NGXMAX,FIMIN,FIMAX,NGYMAX
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
         CALL MPBCIA(NGXMAX)
         CALL MPBCDA(FRMIN)
         CALL MPBCDA(FRMAX)
         CALL MPBCIA(NGYMAX)
         CALL MPBCDA(FIMIN)
         CALL MPBCDA(FIMAX)
C
         DELTFR=(FRMAX-FRMIN)/(NGXMAX-1)
         DELTFI=(FIMAX-FIMIN)/(NGYMAX-1)
         DO NGX=1,NGXMAX
            FR=FRMIN+(NGX-1)*DELTFR
            GX(NGX)=GUCLIP(FR)
            DO NGY=1,NGYMAX
               FI=FIMIN+(NGY-1)*DELTFI
               CALL DIAMIN(FR,FI,AMPL)
               IF(MYRANK.EQ.0) THEN
               IF(LISTEG.GE.1) THEN
               WRITE(6,'(A,1P3E12.4)') '      FR,FI,AMPL = ',FR,FI,AMPL
               CALL GUFLSH
               ENDIF
               ENDIF
               GY(NGY)=GUCLIP(FI)
               GZ(NGX,NGY)=GUCLIP(AMPL)
            ENDDO
         ENDDO
C
         DO NINF=1,NINFM
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
               IF(GZ(NGXN,NGY).GE.GAMPL.AND.
     &            GZ(NGXP,NGY).GE.GAMPL.AND.
     &            GZ(NGX,NGYN).GE.GAMPL.AND.
     &            GZ(NGX,NGYP).GE.GAMPL) THEN
                  DO NINF=1,NINFM
                     IF(GFINF(3,NINF).GT.GAMPL) THEN
                        DO NINF1=NINFM,NINF+1,-1
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
 1000             CONTINUE
               ENDIF
            ENDDO
         ENDDO
         IF(GFINF(3,1).LT.1.E30) THEN
            FRINI=DBLE(GFINF(1,1))
            FIINI=DBLE(GFINF(2,1))
            AMPMIN=DBLE(GFINF(3,1))
         ENDIF
         IF(MYRANK.EQ.0) THEN
            WRITE(6,'(A,1P3E12.4)') '      FR,FI,AMPMIN=',
     &                              FRINI,FIINI,AMPMIN
            CALL GUFLSH
         ENDIF
C
         CALL WMEG2D(GX,GY,GZ,NGZM,NGXMAX,NGYMAX,
     &               KA,GFINF,NINFM)
C
      GOTO 1
C
 9000 CONTINUE
      RETURN
      END
C
      SUBROUTINE WMEG2D(GX,GY,GZ,NGZM,NGXMAX,NGYMAX,
     &                  KA,GFINF,NINFM)
C
      COMMON /WMEGN7/ NCONT,ILN1,IBL1,ICL1,ILN2,IBL2,ICL2
      DIMENSION GX(NGZM),GY(NGZM),GZ(NGZM,NGZM)
      DIMENSION GFINF(3,NINFM),KA(8,NGZM,NGZM)
      DIMENSION RGBS(3,5)
      DATA RGBS/0.0,1.0,1.0,
     &          0.0,0.0,1.0,
     &          0.0,1.0,0.0,
     &          1.0,0.5,0.0,
     &          1.0,0.0,0.0/
C
      GXMIN=GX(1)
      GXMAX=GX(NGXMAX)
      GYMIN=GY(1)
      GYMAX=GY(NGYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXMIN1,GXMAX1,GXSCL)
      CALL GQSCAL(GYMIN,GYMAX,GYMIN1,GYMAX1,GYSCL)
C
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
C
      CALL GMNMX2(GZ,NGZM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
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
C
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
C
      IPRD=0
      GZW=GZMAX-GZMIN
C
C      WRITE(6,*) '2D ',GZMIN,GZMAX,IGMAX,10.0**IGMAX,10.0**(IGMAX-4)
      DO I=1,5
         CALL SETRGB(RGBS(1,I),RGBS(2,I),RGBS(3,I))
         GL=10.0**(IGMAX-I+1)
         IF(10*GL.GT.GZMIN.AND.GL.LT.GZMAX) THEN
            GZORG=GL
            CALL SETLNW(0.035)
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,0,KA)
            GZORG=1.5*GL
            CALL SETLNW(0.018)
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,2,KA)
            GZORG=2.0*GL
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,2,KA)
            GZORG=3.0*GL
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,2,KA)
            GZORG=5.0*GL
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,2,KA)
            GZORG=7.0*GL
            CALL CONTQ1(GZ,NGZM,NGXMAX,NGYMAX,
     &                  GZORG,GZW,1,IPRD,2,KA)
         ENDIF
      ENDDO
C
      CALL SETRGB(0.0,1.0,0.0)
      CALL SETVEW( 4.0,24.0, 3.0,16.0,GXMIN,GXMAX,GYMIN,GYMAX)
      DO NINF=1,NINFM
         IF(GFINF(3,NINF).LT.1.E30) THEN
            CALL MOVE(GFINF(1,NINF),GFINF(2,NINF))
            CALL NUMBI(NINF,'(I1)',1)
         ENDIF
      ENDDO
      CALL OFFVEW
C
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
C
      CALL PAGEE
      RETURN
      END
C
C     ***** FIND EIGEN VALUE *****
C
      SUBROUTINE WMEIGN(KID,LINE)
C
      INCLUDE 'wmcomm.inc'
      CHARACTER LINE*80,KID*1
      DIMENSION XA(2),WORK(2,2)
      EXTERNAL DIAMIN,WMPARM
C
      MODE=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
    2       IF(MODE.EQ.0) WRITE(6,'(A,1P2E12.4)')
     &      ' ## FRINI,FIINI=',FRINI,FIINI
            WRITE(6,'(A)') ' ## INPUT : KID or FRINI,FINII'
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
            IF(MODE.EQ.0)
     &      READ(LINE,*,ERR=2,END=2) FRINI,FIINI
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
         CALL MPBCDA(FRINI)
         CALL MPBCDA(FIINI)
C
         X=FRINI
         Y=FIINI
         IF (MODENW.EQ.0) THEN
            CALL NEWTN0(DIAMIN,X,Y,XX,YY,
     &                  DLTNW,EPSNW,LMAXNW,LISTNW,MYRANK,IERR)
         ELSEIF(MODENW.EQ.1) THEN
            CALL NEWTN1(DIAMIN,X,Y,XX,YY,
     &                  DLTNW,EPSNW,LMAXNW,LISTNW,MYRANK,IERR)
         ELSE
            XX=X
            YY=Y
         ENDIF
C
         FRINI=XX
         FIINI=YY
         CALL WMBFLD

      GOTO 1
C
 9000 CONTINUE
      RETURN
      END
C
C     ***** PARAMETER SCAN *****
C
      SUBROUTINE WMSCAN(KID,LINE)
C
      INCLUDE 'wmcomm.inc'
      DIMENSION PNSAVE(NSM),PTPRSAVE(NSM),PTPPSAVE(NSM)
      DIMENSION PNITBSAVE(NSM),PTITBSAVE(NSM),PUITBSAVE(NSM)
      DIMENSION PUSAVE(NSM)
      DIMENSION GX(NGZM),GY(NGZM,3)
      DIMENSION XA(2),WORK(2,2)
      CHARACTER LINE*80,KID*1,KV*4
      EXTERNAL DIAMIN,WMPARM
C
      MODE=0
      ISCAN=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
    2       IF(MODE.EQ.0) THEN
               WRITE(6,'(A)') 
     &      ' ## ISCAN: 1:NPH0 2:RR 3:QA 4:BB 5:PNA 6:PNAL 7:PT 8:PN'
               WRITE(6,'(A)') 
     &      ' ##        9:qmin 10:rhomin 11:itb 12:rhoitb 13:PU'
            ENDIF
            WRITE(6,'(A)') ' ## INPUT : KID or ISCAN'
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
            IF(MODE.EQ.0)
     &      READ(LINE,*,ERR=2,END=2) ISCAN
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
         CALL MPBCIA(ISCAN)
         IF(ISCAN.EQ.0.OR.ISCAN.EQ.9) GOTO 1
C
         IF(ISCAN.EQ.1) THEN
            KV='NPH0'
            NPH0SAVE=NPH0
            SCMIN=DBLE(NPH0)
         ELSE IF(ISCAN.EQ.2) THEN
            KV='RR  '
            RRSAVE=RR
            SCMIN=RR
         ELSE IF(ISCAN.EQ.3) THEN
            KV='QA  '
            QASAVE=QA
            SCMIN=QA
         ELSE IF(ISCAN.EQ.4) THEN
            KV='BB  '
            BBSAVE=BB
            SCMIN=BB
         ELSE IF(ISCAN.EQ.5) THEN
            KV='PNA '
            PNASAVE=PNA
            SCMIN=PNA
         ELSE IF(ISCAN.EQ.6) THEN
            KV='PNAL'
            PNALSAVE=PNAL
            SCMIN=PNAL
         ELSE IF(ISCAN.EQ.7) THEN
            KV='PT  '
            DO NS=1,NSMAX
               PTPRSAVE(NS)=PTPR(NS)
               PTPPSAVE(NS)=PTPP(NS)
            ENDDO
            SCMIN=PTPR(1)
         ELSE IF(ISCAN.EQ.8) THEN
            KV='PN  '
            DO NS=1,NSMAX
               PNSAVE(NS)=PN(NS)
            ENDDO
            SCMIN=PN(1)
         ELSE IF(ISCAN.EQ.9) THEN
            KV='QMIN'
            QMINSAVE=QMIN
            SCMIN=QMIN
         ELSE IF(ISCAN.EQ.10) THEN
            KV='RHOM'
            RHOMINSAVE=RHOMIN
            SCMIN=RHOMIN
         ELSE IF(ISCAN.EQ.11) THEN
            KV='PNTB'
            DO NS=1,NSMAX
               PNITBSAVE(NS)=PNITB(NS)
               PTITBSAVE(NS)=PTITB(NS)
               PUITBSAVE(NS)=PUITB(NS)
            ENDDO
            SCMIN=0.D0
            SCMAX=1.D0
         ELSE IF(ISCAN.EQ.12) THEN
            KV='RHOB'
            RHOITBSAVE=RHOITB
            SCMIN=RHOITB
         ELSE IF(ISCAN.EQ.13) THEN
            KV='PU  '
            DO NS=1,NSMAX
               PUSAVE(NS)=PU(NS)
            ENDDO
            SCMIN=PU(1)
         END IF
C
      MODE1=0
   11 CONTINUE
         IF(MYRANK.EQ.0) THEN
   12       IF(MODE1.EQ.0) WRITE(6,'(A/1P2E12.4,I5,1P2E12.4)')
     &      ' ## SCMIN,SCMAX,NSCMAX,FRINI,FIINI',
     &           SCMIN,SCMAX,NSCMAX,FRINI,FIINI
            WRITE(6,'(A)') 
     &      ' ## INPUT : KID or SCMIN,SCMAX,NSCMAX,FRINI,FIINI'
            CALL TASK_KLIN(LINE,KID,MODE1,WMPARM)
            IF(MODE1.EQ.0)
     &       READ(LINE,*,ERR=12,END=12) SCMIN,SCMAX,NSCMAX,FRINI,FIINI
         ENDIF
         CALL MPBCIA(MODE1)
         IF(MODE1.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 9000
         ENDIF
         IF(MODE1.EQ.2) THEN
            CALL WMPRBC
            GOTO 11
         ENDIF
         IF(MODE1.EQ.3) GOTO 11
C
         CALL MPBCDA(SCMIN)
         CALL MPBCDA(SCMAX)
         CALL MPBCIA(NSCMAX)
         CALL MPBCDA(FRINI)
         CALL MPBCDA(FIINI)
C
         DSC=(SCMAX-SCMIN)/(NSCMAX-1)
         X=FRINI
         Y=FIINI
C
         DO NSC=1,NSCMAX
            SC=SCMIN+DSC*(NSC-1)
            GX(NSC)=GUCLIP(SC)
         ENDDO
C
         DO NSC=1,NSCMAX
C
            SC=SCMIN+DSC*(NSC-1)
            IF(ISCAN.EQ.1) THEN
               NPH0=NINT(SC)
            ELSE IF(ISCAN.EQ.2) THEN
               RR=SC
            ELSE IF(ISCAN.EQ.3) THEN
               QA=SC
            ELSE IF(ISCAN.EQ.4) THEN
               BB=SC
            ELSE IF(ISCAN.EQ.5) THEN
               PNA=SC
            ELSE IF(ISCAN.EQ.6) THEN
               PNAL=SC
            ELSE IF(ISCAN.EQ.7) THEN
               DO NS=1,NSMAX
                  PTPR(NS)=SC*PTPRSAVE(NS)/PTPRSAVE(1)
                  PTPP(NS)=SC*PTPPSAVE(NS)/PTPRSAVE(1)
               ENDDO
            ELSE IF(ISCAN.EQ.8) THEN
               PN(1)=SC
               DO NS=2,NSMAX
                  PN(NS)=PNSAVE(NS)*SC/PNSAVE(1)
               ENDDO
            ELSE IF(ISCAN.EQ.9) THEN
               QMIN=SC
            ELSE IF(ISCAN.EQ.10) THEN
               RHOMIN=SC
            ELSE IF(ISCAN.EQ.11) THEN
               DO NS=1,NSMAX
                  PNITB(NS)=PNITBSAVE(NS)*SC
                  PTITB(NS)=PTITBSAVE(NS)*SC
                  PUITB(NS)=PUITBSAVE(NS)*SC
               ENDDO
            ELSE IF(ISCAN.EQ.12) THEN
               RHOITB=SC
            ELSE IF(ISCAN.EQ.13) THEN
               IF(PUSAVE(1).EQ.0.D0) THEN
                  DO NS=1,NSMAX
                     PU(NS)=SC
                  ENDDO
               ELSE
                  PU(1)=SC
                  DO NS=2,NSMAX
                     PU(NS)=PUSAVE(NS)*SC/PUSAVE(1)
                  ENDDO
               ENDIF
            END IF
C
            IF(MODENW.EQ.0) THEN
               CALL NEWTN0(DIAMIN,X,Y,XX,YY,
     &                     DLTNW,EPSNW,LMAXNW,LISTNW,MYRANK,IERR)
            ELSEIF(MODENW.EQ.1) THEN
               CALL NEWTN1(DIAMIN,X,Y,XX,YY,
     &                     DLTNW,EPSNW,LMAXNW,LISTNW,MYRANK,IERR)
            ELSE
               XX=X
               YY=Y
            ENDIF
            IF(IERR.NE.0) GOTO 3000
C
            CALL DIAMIN(XX,YY,AMPL)
            GY(NSC,1)=GUCLIP(XX)
            GY(NSC,2)=GUCLIP(YY)
            GY(NSC,3)=GUCLIP(AMPL)
            X=XX
            Y=YY
         ENDDO
 3000    CONTINUE
C
         IF(IERR.EQ.0) THEN
            FRINI=XX
            FIINI=YY
            CALL WMBFLD
         ENDIF
C
         IF(ISCAN.EQ.1) THEN
            NPH0=NPH0SAVE
         ELSE IF(ISCAN.EQ.2) THEN
            RR=RRSAVE
         ELSE IF(ISCAN.EQ.3) THEN
            QA=QASAVE
         ELSE IF(ISCAN.EQ.4) THEN
            BB=BBSAVE
         ELSE IF(ISCAN.EQ.5) THEN
            PNA=PNASAVE
         ELSE IF(ISCAN.EQ.6) THEN
            PNAL=PNALSAVE
         ELSE IF(ISCAN.EQ.7) THEN
            DO NS=1,NSMAX
               PTPR(NS)=PTPRSAVE(NS)
               PTPP(NS)=PTPPSAVE(NS)
            ENDDO
         ELSE IF(ISCAN.EQ.8) THEN
            DO NS=1,NSMAX
               PN(NS)=PNSAVE(NS)
            ENDDO
         ELSE IF(ISCAN.EQ.9) THEN
            QMIN=QMINSAVE
         ELSE IF(ISCAN.EQ.10) THEN
            RHOMIN=RHOMINSAVE
         ELSE IF(ISCAN.EQ.11) THEN
            DO NS=1,NSMAX
               PNITB(NS)=PNITBSAVE(NS)
               PTITB(NS)=PTITBSAVE(NS)
               PUITB(NS)=PUITBSAVE(NS)
            ENDDO
         ELSE IF(ISCAN.EQ.12) THEN
            RHOITB=RHOITBSAVE
         ELSE IF(ISCAN.EQ.13) THEN
            DO NS=1,NSMAX
               PU(NS)=PUSAVE(NS)
            ENDDO
         ENDIF
C
         CALL WMSC1G(GX,GY(1,1),NSCMAX,ISCAN,KV)
         CALL WMSC1G(GX,GY(1,2),NSCMAX,ISCAN,KV)
C
         IF(MYRANK.EQ.0) THEN
            WRITE(6,*) '      '//KV
            WRITE(6,'(I5,1P4E12.4)') 
     &          (I,GX(I),GY(I,1),GY(I,2),GY(I,3),I=1,NSCMAX) 
            CALL GUFLSH
          ENDIF
C
      GOTO 11
C
 9000 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE WMSC1G(GX,GY,NSCMAX,ISCAN,KV)
C
      CHARACTER KV*4
C
      DIMENSION GX(NSCMAX),GY(NSCMAX)
C
      IF(ISCAN.EQ.1) THEN
         NKETA=1
      ELSE IF(ISCAN.EQ.2) THEN
         NKETA=3
      ELSE IF(ISCAN.EQ.3) THEN
         NKETA=2
      ELSE IF(ISCAN.EQ.4) THEN
         NKETA=2
      ELSE IF(ISCAN.EQ.5) THEN
         NKETA=4
      ELSE IF(ISCAN.EQ.6) THEN
         NKETA=4
      ELSE IF(ISCAN.EQ.7) THEN
         NKETA=2
      ELSE
         NKETA=2
      END IF
C
C     ****MAP****
C
      CALL GMNMX1(GX,1,NSCMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX1(GY,1,NSCMAX,1,GYMIN1,GYMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSCL)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSCL)
C
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
C
      CALL SETCHS(0.3,0.0)
      CALL MOVE(4.0,17.0)
      CALL TEXT(KV,4)
C
      CALL PAGEE
      RETURN
C
      END
C
C     ***** CALCULATE AMPLITUDE *****
C
      SUBROUTINE DIAMIN(X,Y,F)
C
      INCLUDE 'wmcomm.inc'
C
      RF=X
      RFI=Y
      CRF=DCMPLX(RF,RFI)
C
      MODEEG=1

      CALL WMSETG(IERR)
      CALL WMSOLV
      CALL WMEFLD
C
      ESUM=0.D0
      EABSMAX=0.D0
      CEMAX=0.D0
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         EABS=ABS(CEFLD(1,NTH,NPH,NR))**2
     &       +ABS(CEFLD(2,NTH,NPH,NR))**2
     &       +ABS(CEFLD(3,NTH,NPH,NR))**2
         ESUM=ESUM+EABS
         EABS=ABS(CEFLD(2,NTH,NPH,NR))**2
         IF(EABS.GT.EABSMAX) THEN
            EABSMAX=EABS
            CEMAX=CEFLD(2,NTH,NPH,NR)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      IF(ABS(CEMAX).EQ.0.D0) THEN
         CFACT=1.D0
      ELSE
         CFACT=1.D0/CEMAX
      ENDIF
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CEFLD(1,NTH,NPH,NR)=CFACT*CEFLD(1,NTH,NPH,NR)
         CEFLD(2,NTH,NPH,NR)=CFACT*CEFLD(2,NTH,NPH,NR)
         CEFLD(3,NTH,NPH,NR)=CFACT*CEFLD(3,NTH,NPH,NR)
      ENDDO
      ENDDO
      ENDDO
C
      EABSMAX=0.D0
      CEMAX=0.D0
      ACEMAX=0.D0
      DO NR=1,NRMAX+1
         EABS2=0.D0
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         EABS2=EABS2+DREAL(CEFLDK(2,MDX,NDX,NR))**2
     &              +DIMAG(CEFLDK(2,MDX,NDX,NR))**2
         IF(ABS(CEFLDK(2,MDX,NDX,NR)).GT.ACEMAX) THEN
            CEMAX=CEFLDK(2,MDX,NDX,NR)
            ACEMAX=ABS(CEMAX)
         ENDIF
      ENDDO
      ENDDO
         IF(EABS2.GT.EABSMAX) EABSMAX=EABS2
      ENDDO
C
      IF(EABSMAX.EQ.0.D0) THEN
         CFACT=1.D0
      ELSE
         CFACT=ABS(CEMAX)/(CEMAX*SQRT(EABSMAX))
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
C
C      WRITE(6,*) NRMAX,NTHMAX,NPHMAX,ESUM,RF
      F=NRMAX*NTHMAX*NPHMAX/(ESUM*RF**4)
      IF(F.LE.1.D-15) F=1.D-15
      AMPEIGEN=F
C
      RETURN
      END
C
C     ****** TWO-DIMENSIONAL NEWTON METHOD ******
C
      SUBROUTINE NEWTN0(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,MYRANK,IER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LPRINT,LDUMP
      EXTERNAL SUB
C
      LDUMP=MYRANK.EQ.0.AND.LIST.GE.2
      LPRINT=MYRANK.EQ.0.AND.LIST.EQ.1
C
      IER=0
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
C
      CALL SUB(X,   Y,   F00)
      IF(LPRINT) WRITE(6,600) X,Y,0.D0,0.D0,F00
      IF(LPRINT) CALL GUFLSH
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y,F00
      IF(LDUMP) CALL GUFLSH
      IF(F00.LE.2.D-15) GOTO 9001
C
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
C
      FX=(FP0-FM0)/(2*HX)
      FY=(F0P-F0M)/(2*HY)
C
    1 CONTINUE
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
C
      FXX=(FP0-2*F00+FM0)/(HX*HX)
      FYY=(F0P-2*F00+F0M)/(HY*HY)
      FXY=(FPP-FPM-FMP+FMM)/(4*HX*HY)
      IF(LDUMP) WRITE(6,'(A,1P2E12.4)') 'FX,FY =',FX,FY
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'FXX,FXY,FXY =',
     &                                         FXX,FXY,FYY
      IF(LDUMP) CALL GUFLSH
C
      DF=SQRT(FX*FX+FY*FY)
      DET=FXX*FYY-FXY*FXY
      IF(LDUMP) WRITE(6,'(A,1P2E12.4)') 'DF,DET:0 =',DF,DET
      IF(LDUMP) CALL GUFLSH
      IF(DF.LE.2.D-15) GO TO 9000
      CALL MPSYNC
C
      H11= FYY/DET
      H12=-FXY/DET
      H21=-FXY/DET
      H22= FXX/DET
      DX=-(H11*FX+H12*FY)
      DY=-(H21*FX+H22*FY)
      V=SQRT(X*X+Y*Y+EPS)
      DV=SQRT(DX*DX+DY*DY)
      IF(DV/V.LE.EPS) GO TO 9000
C
      TT=1.D0
    2 X=X+TT*DX
      Y=Y+TT*DY
C
      CALL SUB(X,   Y,   F00)
      IF(LPRINT) WRITE(6,600) X,Y,TT*DX,TT*DY,F00
      IF(LPRINT) CALL GUFLSH
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y,F00
      IF(LDUMP) CALL GUFLSH
      IF(F00.LE.2.D-15) GOTO 9001
      CALL MPSYNC
C
      CALL SUB(X-HX,Y,   FM0)
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X-HX,Y,FM0
      IF(LDUMP) CALL GUFLSH
      CALL MPSYNC
      CALL SUB(X+HX,Y,   FP0)
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X+HX,Y,FP0
      IF(LDUMP) CALL GUFLSH
      CALL MPSYNC
      CALL SUB(X,   Y-HY,F0M)
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y-HY,F0M
      IF(LDUMP) CALL GUFLSH
      CALL SUB(X,   Y+HY,F0P)
      IF(LDUMP) WRITE(6,'(A,1P3E12.4)') 'X,Y,F =',X,Y+HY,F0P
      IF(LDUMP) CALL GUFLSH
C
      FXN=(FP0-FM0)/(2*HX)
      FYN=(F0P-F0M)/(2*HY)
C
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
C
      IER=2
      IF(LPRINT) 
     &   WRITE(6,*) 'XX NEWTN0: LOOP COUNT EXCEEDS UPPER BOUND.'
      GOTO 9001
C
 8000 IER=1
      IF(LPRINT)
     &   WRITE(6,*) 'XX NEWTN0: DOES NOT CONVERGE.'
      GOTO 9001
C
 9000 CALL SUB(X,Y,F00)
 9001 XX=X
      YY=Y
      RETURN
  600 FORMAT(' ',3X,'X,Y,DX,DY,F = ',1P2E15.7,1P3E10.2)
      END
C
C     ****** TWO-DIMENSIONAL NEWTON METHOD (SIMPLER) ******
C
      SUBROUTINE NEWTN1(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,MYRANK,IER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SUB
C
      IER=0
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
C
    1 CONTINUE
      CALL SUB(X+HX,Y+HY,FPP)
C
      FXX=(FP0-2*F00+FM0)/(HX*HX)
      FYY=(F0P-2*F00+F0M)/(HY*HY)
      FXY=(FPP-FP0-F0P+F00)/(HX*HY)
      DF=SQRT(FX*FX+FY*FY)
      DET=FXX*FYY-FXY*FXY
      H11= FYY/DET
      H12=-FXY/DET
      H21=-FXY/DET
      H22= FXX/DET
C
      DX=-(H11*FX+H12*FY)
      DY=-(H21*FX+H22*FY)
      TT=1.D0
    2 X=X+TT*DX
      Y=Y+TT*DY
C
      CALL SUB(X,   Y,   F00)
      CALL SUB(X-HX,Y,   FM0)
      CALL SUB(X+HX,Y,   FP0)
      CALL SUB(X,   Y-HY,F0M)
      CALL SUB(X,   Y+HY,F0P)
      FXN=(FP0-FM0)/(2*HX)
      FYN=(F0P-F0M)/(2*HY)
      IF(MYRANK.EQ.0) THEN
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
C         IF(DF.LE.EPS) GO TO 9000
         IF(ITER.LE.ILMAX) GO TO 1
      ENDIF
C
      IER=2
      IF(MYRANK.EQ.0) THEN
      IF(LIST.GT.0)
     &   WRITE(6,*) 'XX NEWTN2: LOOP COUNT EXCEEDS UPPER BOUND.'
      ENDIF
      GOTO 9000
C
 8000 IER=1
      IF(MYRANK.EQ.0) THEN
      IF(LIST.GT.0)
     &   WRITE(6,*) 'XX NEWTN2: DOES NOT CONVERGE.'
      ENDIF
      GOTO 9000
C
 9000 XX=X
      YY=Y
      RETURN
  600 FORMAT(' ',3X,'X,Y,DX,DY,F = ',1P2E15.7,1P3E10.2)
      END
C
      SUBROUTINE PRINTFG(LPRINT, ITER, F, X, G, N)
       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       DIMENSION X(N), G(N)
      WRITE(6, '(1X,A,I4,A,1PE15.7)') 'ITER=', ITER, '  F=',F 
      IF(LPRINT .GE. 2) WRITE(6, '(3X,A,1P5E15.7/(5X,5E15.7))')
     &                      'X=',(X(I),I=1,N)
      IF(LPRINT .GE. 3) WRITE(6, '(3X,A,1P5E15.7/(5X,5E15.7))')
     &                      'G=',(G(I),I=1,N)
      RETURN
      END
C
C     ***** WRITE WAVE DATA *****
C
      SUBROUTINE WMWOUT
C
      INCLUDE 'wmcomm.inc'
C
      NF=26
      WRITE(NF,*) 'FR [MHz] = ',FRINI
      WRITE(NF,*) 'FI [MHz] = ',FIINI
      WRITE(NF,*) 'AMPN     = ',AMPEIGEN
      WRITE(NF,*) 'NR COUNT = ',NRMAX+1
      WRITE(NF,*) 'MD COUNT = ',MDMAX-MDMIN+1
      WRITE(NF,*) 'ND COUNT = ',NDMAX-NDMIN+1
      WRITE(NF,*)
C
      WRITE(NF,*) 'r [m]:'
      WRITE(NF,'(1P5E14.6)') (XR(NR),NR=1,NRMAX+1)
      WRITE(NF,*)
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            WRITE(NF,*) 'MD           = ',MD+NTH0
            WRITE(NF,*) 'ND           = ',ND+NPH0
            WRITE(NF,*) 'Re Er:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DREAL(CEFLDK(1,MDX,NDX,NR)),NR=1,NRMAX+1)
            WRITE(NF,*) 'Im Er:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DIMAG(CEFLDK(1,MDX,NDX,NR)),NR=1,NRMAX+1)
            WRITE(NF,*) 'Re Etheta:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DREAL(CEFLDK(2,MDX,NDX,NR)),NR=1,NRMAX+1)
            WRITE(NF,*) 'Im Etheta:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DIMAG(CEFLDK(2,MDX,NDX,NR)),NR=1,NRMAX+1)
            WRITE(NF,*) 'Re Ephi:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DREAL(CEFLDK(3,MDX,NDX,NR)),NR=1,NRMAX+1)
            WRITE(NF,*) 'Im Ephi:'
            WRITE(NF,'(1P5E14.6)') 
     &           (DIMAG(CEFLDK(3,MDX,NDX,NR)),NR=1,NRMAX+1)
         ENDDO
      ENDDO
C
      RETURN
      END
