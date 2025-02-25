C
C     ***********************************************
C     ****** GSAF APLLICATION V3.5 : PLOT 2    ******
C     ***********************************************
C
C
C     ****** POINT PLOT : XY, VARIABLE STEP, PATTERN ******
C
      SUBROUTINE PLOTG2(Z,X,Y,NFMIN,NFMAX,NFSTEP,
     &                  VL,RGBL,ILN,WLN,NLMAX,
     &                  VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX,IPRD)
C
      IMPLICIT LOGICAL(L)
C
      EXTERNAL PLOTV1
      DIMENSION Z(NFMAX),X(NFMAX),Y(NFMAX)
      DIMENSION VL(NLMAX),RGBL(3,NLMAX),ILN(NLMAX),WLN(NLMAX)
      DIMENSION VM(NMMAX),RGBM(3,NMMAX),SIZEM(NMMAX),WLM(NMMAX)
      INTEGER MARKM(NMMAX),STEPM(NMMAX)
C
      CALL PLOTG0(Z,X,Y,NFMIN,NFMAX,NFSTEP,
     &            VL,RGBL,ILN,WLN,NLMAX,
     &            VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX,IPRD,3,PLOTV1)
C
      RETURN
      END
C
C     ****** CONTOUR PLOT : R-THETA, VARIABLE R, PATTERN ******
C
      SUBROUTINE PLOTG3(Z,R,TH,NFMIN,NFMAX,NFSTEP,
     &                  VL,RGBL,ILN,WLN,NLMAX,
     &                  VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSCTR5/ RMAX,RT,TT,XT,YT,NTHMAX
C
      EXTERNAL PLOTV2
      DIMENSION Z(NFMAX),R(NFMAX),TH(NFMAX),T(2)
      DIMENSION VL(NLMAX),RGBL(3,NLMAX),ILN(NLMAX),WLN(NLMAX)
      DIMENSION VM(NMMAX),RGBM(3,NMMAX),SIZEM(NMMAX),WLM(NMMAX)
      INTEGER MARKM(NMMAX),STEPM(NMMAX)
C
      RMAX=R(1)
      DO NF=2,NFMAX
         RMAX=MAX(RMAX,R(NF))
      END DO
      NTHMAX=120
      T(1)=2*3.1415926/NTHMAX
      T(2)=0.0
C
      CALL PLOTG0(Z,R,TH,NFMIN,NFMAX,NFSTEP,
     &            VL,RGBL,ILN,WLN,NLMAX,
     &            VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX,2,1,PLOTV2)
C
      RETURN
      END
C
C     ****** CONTOUR PLOT : R-THETA, VARIABLE STEP, PATTERN ******
C
      SUBROUTINE PLOTG4(Z,R,TH,NFMIN,NFMAX,NFSTEP,
     &                  VL,RGBL,ILN,WLN,NLMAX,
     &                  VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX,IPRD)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSCTR5/ RMAX,RT,TT,XT,YT,NTHMAX
C
      EXTERNAL PLOTV2
      DIMENSION Z(NFMAX),R(NFMAX),TH(NFMAX)
      DIMENSION VL(NLMAX),RGBL(3,NLMAX),ILN(NLMAX),WLN(NLMAX)
      DIMENSION VM(NMMAX),RGBM(3,NMMAX),SIZEM(NMMAX),WLM(NMMAX)
      INTEGER MARKM(NMMAX),STEPM(NMMAX)
C
      RMAX=R(NXMAX)
      NTHMAX=120
C
      CALL PLOTG0(Z,R,TH,NFMIN,NFMAX,NFSTEP,
     &                  VL,RGBL,ILN,WLN,NLMAX,
     &                  VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX,2,1,PLOTV2)
      RETURN
      END
C
C     ****** CONTOUR PLOT : COMMON SUB ******
C
      SUBROUTINE PLOTG0(Z,X,Y,NFMIN,NFMAX,NFSTEP,
     &                  VL,RGBL,ILN,WLN,NLMAX,
     &                  VM,RGBM,MARKM,STEPM,SIZEM,WLM,NMMAX, 
     &                  IPRD,INDX,SUBV)
C
      IMPLICIT LOGICAL(L)
      EXTERNAL SUBV
      DIMENSION Z(NFMAX),X(NFMAX),Y(NFMAX)
      DIMENSION VL(NLMAX),RGBL(3,NLMAX),ILN(NLMAX),WLN(NLMAX)
      DIMENSION VM(NMMAX),RGBM(3,NMMAX),SIZEM(NMMAX),WLM(NMMAX)
      INTEGER MARKM(NMMAX),STEPM(NMMAX)
C
      CALL INQMRK(IMARKS,HMRKS,WMRKS,ANGLS,TILTS)
      CALL INQRGB(RS,GS,BS)
      CALL INQLNW(WS)
C
      IF(NLMAX.GT.0) THEN
         DO NF=NFMIN,NFMAX,NFSTEP
            CALL SUBV(X(NF),Y(NF),XP,YP)
            IF(NF.EQ.NFMIN) THEN
               IF(NLMAX.EQ.1) THEN
                  NL=1
               ELSE
                  CALL SETIDX_PLOT(Z(NF),VL,NLMAX,NL)
               END IF
               CALL SETRGB(RGBM(1,NL),RGBM(2,NL),RGBM(3,NL))
               CALL SETLNW(WLN(NL))
               CALL MOVEPT(XP,IP,ILN(NL))
            ELSE
               IF(NLMAX.EQ.1) THEN
                  NL=1
               ELSE
                  CALL SETIDX_PLOT(Z(NF),VL,NLMAX,NL)
                  CALL SETRGB(RGBM(1,NL),RGBM(2,NL),RGBM(3,NL))
                  CALL SETLNW(WLN(NL))
               END IF
               CALL DRAWPT(XP,IP,ILN(NL))
            END IF
         END DO
      END IF
C
      MARKM(1)=1
      IF(NMMAX.GT.0) THEN
         IF(NMMAX <= 1) THEN
            NM=1
            CALL SETMRK(MARKM(NM),SIZEM(NM),SIZEM(NM),0.0,0.0)
            CALL SETRGB(RGBM(1,NM),RGBM(2,NM),RGBM(3,NM))
            CALL SETLNW(WLM(NM))
         END IF
C
         DO NF=NFMIN,NFMAX,NFSTEP
            CALL SUBV(X(NF),Y(NF),XP,YP)
            IF(NMMAX.GT.1) THEN
               CALL SETIDX_PLOT(Z(NF),VM,NMMAX,NM)
               CALL SETMRK(MARKM(NM),SIZEM(NM),SIZEM(NM),0.0,0.0)
               CALL SETRGB(RGBM(1,NM),RGBM(2,NM),RGBM(3,NM))
               CALL SETLNW(WLM(NM))
            END IF
            CALL MARK(XP,YP)
         END DO
      END IF

      CALL SETMRK(IMARKS,HMRKS,WMRKS,ANGLS,TILTS)
      CALL SETRGB(RS,GS,BS)
      CALL SETLNW(WS)

      RETURN
      END
C
C     ****** CONTOUR PLOT SLAVE ROUTINE : CONV XY ******
C
      SUBROUTINE PLOTV1(XA,YA,XB,YB)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
C
      XB=DX*(XA-GXS)+PXS
      YB=DY*(YA-GYS)+PYS
      RETURN
      END
C
C     ****** CONTOUR PLOT SLAVE ROUTINE : CONV RT ******
C
      SUBROUTINE PLOTV2(RA,TA,XB,YB)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
      COMMON /GSCTR5/ RMAX,RT,TT,XT,YT,NTHMAX

      XB=DX*(RA(I)*COS(TA)-GXS)+PXS
      YB=DY*(RA(I)*SIN(TA)-GYS)+PYS
      RETURN
      END
C
C     ****** CONTOUR PLOT SLAVE ROUTINE : CONV RT ******
C
      SUBROUTINE SETIDX_PLOT(V,VL,NLMAX,NL)
C
      DIMENSION VL(NLMAX)
C
      NL=1
      DO N=2,NLMAX
         VB=0.5*(VL(N-1)+VL(N))
         IF(V.LT.VB) GOTO 1000
         NL=N
      END DO
 1000 CONTINUE
      RETURN
      END
