MODULE dpgsub

CONTAINS

!--------------------------------------------------

  SUBROUTINE DPFPGRA(PMAX,NTHM,NPMAX,NTHMAX,P,TH,F,STRING)

      IMPLICIT NONE
      REAL(8),INTENT(IN):: PMAX
      INTEGER(4),INTENT(IN):: NPMAX,NTHMAX
      real(8),dimension(NPMAX),INTENT(IN):: P
      real(8),dimension(NTHMAX),INTENT(IN):: TH
      real(8),DIMENSION(NTHM,NPMAX),INTENT(IN):: F
      CHARACTER(LEN=*),INTENT(IN):: STRING

      real(4),DIMENSION(NPMAX,NTHMAX):: GF
      real(4),dimension(NPMAX):: GP
      real(4),dimension(NTHMAX):: GTH

      INTEGER,PARAMETER:: NGLM=100
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      real(4):: GPMAX, GPMIN1, GPMAX1, GPSTEP
      real(4):: GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      integer:: NP,NTH
      integer:: NTHM,NGLMAX,NGL
      real(4):: GUCLIP
      integer:: NGULEN

      DO NP=1,NPMAX
         GP(NP)=GUCLIP(P(NP))
      END DO
      DO NTH=1,NTHMAX
         GTH(NTH)=GUCLIP(TH(NTH))
      END DO
      DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            GF(NP,NTH)=GUCLIP(-LOG10(ABS(F(NTH,NP))))
         END DO
      END DO

      GPMAX=GUCLIP(PMAX)

      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)

      CALL GMNMX2(GF,NPMAX,1,NPMAX,1,1,NTHMAX,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)

      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))

      NGLMAX=INT(SQRT((GFMAX1-GFMIN1)/(0.01*GFSTEP)))
      DO NGL=1,NGLMAX
         ZL(NGL)=GFMIN1+(0.1*GFSTEP*(NGL-1))**2
         RGB(1,NGL)=1.D0
         RGB(2,NGL)=0.9*FLOAT(NGL-1)/FLOAT(NGLMAX-1)
         RGB(3,NGL)=0.9*FLOAT(NGL-1)/FLOAT(NGLMAX-1)
         ILN(NGL)=0
         WLN(NGL)=0.07
      ENDDO
      CALL CONTG4X(GF,GP,GTH,NPMAX,NPMAX,NTHMAX,ZL,RGB,ILN,WLN,NGLMAX,0)

      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)

      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      CALL MOVE(18.0,12.5)
      CALL TEXT('STEP =',6)
      CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)

      CALL PAGEE

      RETURN
  END SUBROUTINE DPFPGRA

!     ****** CONTOUR PLOT : R-THETA, VARIABLE STEP, PATTERN ******

  SUBROUTINE CONTG4X(Z,R,T,NXA,NXMAX,NYMAX, &
                     ZL,RGB,ILN,WLN,NSTEP,ISPL)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NXA,NXMAX,NYMAX,NSTEP,ISPL
    REAL(4),INTENT(IN):: Z(NXA,NYMAX),R(NXMAX),T(NYMAX)
    REAL(4),INTENT(IN):: ZL(NSTEP),RGB(3,NSTEP),WLN(NSTEP)
    INTEGER,INTENT(IN):: ILN(NSTEP)
    COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
    REAL(4):: RMAX,RT,TT,XT,YT

      RMAX=R(NXMAX)

      CALL CONTG0X(Z,R,T,NXA,NXMAX,NYMAX, &
                   ZL,RGB,ILN,WLN,NSTEP,ISPL,0,3,CONTV2X)

      RETURN
  END SUBROUTINE CONTG4X

!     ****** CONTOUR PLOT : COMMON SUB ******

  SUBROUTINE CONTG0X(Z,X,Y,NXA,NXMAX,NYMAX, &
                     ZL,RGB,ILN,WLN,NSTEP, &
                     ISPL,IPRD,INDX,SUBV)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NXA,NXMAX,NYMAX,NSTEP,ISPL,IPRD,INDX
    REAL(4),INTENT(IN):: Z(NXA,NYMAX),X(NXMAX),Y(NYMAX)
    REAL(4),INTENT(IN):: ZL(NSTEP),RGB(3,NSTEP),WLN(NSTEP)
    INTEGER,INTENT(IN):: ILN(NSTEP)
    INTEGER:: KA(2,NXMAX*NYMAX)
    EXTERNAL SUBV
    INTEGER,PARAMETER:: NH=101
    REAL(4):: ZLS(NH),RGBS(3,NH),WLNS(NH)
    INTEGER:: ILNS(NH)
!      PARAMETER (NFMAX=2000,NGMAX=4000)
    INTEGER,PARAMETER:: NFMAX=200,NGMAX=400
    REAL(4):: XF(NFMAX),YF(NFMAX)
    REAL(4):: XP(NFMAX),YP(NFMAX)
    REAL(4):: XG(NGMAX),YG(NGMAX)
    REAL(4):: RS,GS,BS,WS,ZORG,U0,U1,U2,U3,U4,UMAX,UMIN
    REAL(4):: XA,XB,YA,YB,RT
    INTEGER:: KMAX,I,NXM,NYM,NX,NY,IE,NXP,NYP,K,J,NF,NP,NN
    INTEGER:: IEL,NXL,NYL,NXL1,NYL1,IE1
    INTEGER:: N1X,N1Y,N2X,N2Y,N3X,N3Y,N4X,N4Y,N4,I2X,I2Y,I3X,I3Y,I4X,I4Y
    INTEGER:: NAX,NAY,NBX,NBY,UA,UB,MODE
    INTEGER:: NSAX,NSAY,NSBX,NSBY,USA,USB,MODES,MODEL
    LOGICAL:: LINV,LEND

      IF(ISPL.GE.0) THEN
         CALL INQRGB(RS,GS,BS)
         CALL INQLNW(WS)
      ENDIF

      KMAX=NSTEP
      IF(KMAX.GT.NH) KMAX=NH

      IF(ZL(1).LE.ZL(KMAX)) THEN
         ZORG=ZL(1)
         DO I=1,KMAX
            ZLS(I)=ZL(I)-ZORG
         END DO
      ELSE
         ZORG=ZL(KMAX)
         DO I=1,KMAX
            ZLS(I)=ZL(KMAX-I+1)-ZORG
         END DO
      ENDIF

      IF(ISPL.GE.0) THEN
         IF(ZL(1).LE.ZL(KMAX)) THEN
            DO I=1,KMAX
               RGBS(1,I)=RGB(1,I)
               RGBS(2,I)=RGB(2,I)
               RGBS(3,I)=RGB(3,I)
               ILNS(I)=ILN(I)
               WLNS(I)=WLN(I)
            END DO
         ELSE
            DO I=1,KMAX
               RGBS(1,I)=RGB(1,KMAX-I+1)
               RGBS(2,I)=RGB(2,KMAX-I+1)
               RGBS(3,I)=RGB(3,KMAX-I+1)
               ILNS(I)=ILN(KMAX-I+1)
               WLNS(I)=WLN(KMAX-I+1)
            END DO
         ENDIF
      ELSE
         DO I=1,KMAX
            ILNS(I)=ILN(1)
            WLNS(I)=WLN(1)
         END DO
      ENDIF

      NXM=NXMAX-1
      NYM=NYMAX-1
      IF(IPRD.EQ.1.OR.IPRD.EQ.3) NXM=NXMAX
      IF(IPRD.EQ.2.OR.IPRD.EQ.3) NYM=NYMAX

      DO NY=1,NYM
      DO NX=1,NXM
         IE=NXM*(NY-1)+NX
         NXP=NX+1
         NYP=NY+1
         IF(NXP.GT.NXMAX) NXP=1
         IF(NYP.GT.NYMAX) NYP=1
         U1=Z(NX ,NY )-ZORG
         U2=Z(NXP,NY )-ZORG
         U3=Z(NXP,NYP)-ZORG
         U4=Z(NX ,NYP)-ZORG
         UMAX=MAX(U1,U2,U3,U4)
         UMIN=MIN(U1,U2,U3,U4)
         KA(1,IE)=1
         DO K=1,KMAX
            IF(ZLS(K).LT.UMIN) KA(1,IE)=K+1
            IF(ZLS(K).LE.UMAX) KA(2,IE)=K
         END DO
      END DO
      END DO

      DO K=1,KMAX
         U0=ZLS(K)
         IF(ISPL.GE.0) THEN
            CALL SETRGB(RGBS(1,K),RGBS(2,K),RGBS(3,K))
            CALL SETLNW(WLNS(K))
         ENDIF
      DO NY=1,NYM
      DO NX=1,NXM
         IE=NXM*(NY-1)+NX
         IF(KA(1,IE).LE.K.AND.KA(2,IE).GE.K) THEN
            N1X=NX
            N1Y=NY
            N2X=NX+1
            N2Y=NY
            N3X=NX+1
            N3Y=NY+1
            N4X=NX
            N4Y=NY+1
!     
            U1=Z(N1X,N1Y)-ZORG
            I2X=N2X
            I2Y=N2Y
            IF(I2X.GT.NXMAX) I2X=1
            IF(I2Y.GT.NYMAX) I2Y=1
            U2=Z(I2X,I2Y)-ZORG
            I3X=N3X
            I3Y=N3Y
            IF(I3X.GT.NXMAX) I3X=1
            IF(I3Y.GT.NYMAX) I3Y=1
            U3=Z(I3X,I3Y)-ZORG
            I4X=N4X
            I4Y=N4Y
            IF(I4X.GT.NXMAX) I4X=1
            IF(I4Y.GT.NYMAX) I4Y=1
            U4=Z(I4X,I4Y)-ZORG

            IF((U1.GT.U0.AND.U2.LE.U0).OR. &
               (U1.LE.U0.AND.U2.GT.U0)) THEN
               NAX=N1X
               NAY=N1Y
               UA=U1
               NBX=N2X
               NBY=N2Y
               UB=U2
               MODE=1
               IF((U2.GT.U0.AND.U3.LE.U0).OR. &
                  (U2.LE.U0.AND.U3.GT.U0)) THEN
                  NSAX=N2X
                  NSAY=N2Y
                  USA=U2
                  NSBX=N3X
                  NSBY=N3Y
                  USB=U3
                  MODES=2
               ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR. &
                      (U3.LE.U0.AND.U4.GT.U0)) THEN
                  NSAX=N3X
                  NSAY=N3Y
                  USA=U3
                  NSBX=N4X
                  NSBY=N4Y
                  USB=U4
                  MODES=3
               ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR. &
                      (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR. &
                   (U2.LE.U0.AND.U3.GT.U0)) THEN
               NAX=N2X
               NAY=N2Y
               UA=U2
               NBX=N3X
               NBY=N3Y
               UB=U3
               MODE=2
               IF((U3.GT.U0.AND.U4.LE.U0).OR. &
                  (U3.LE.U0.AND.U4.GT.U0)) THEN
                  NSAX=N3X
                  NSAY=N3Y
                  USA=U3
                  NSBX=N4X
                  NSBY=N4Y
                  USB=U4
                  MODES=3
               ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR. &
                      (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR. &
                   (U3.LE.U0.AND.U4.GT.U0)) THEN
               NAX=N3X
               NAY=N3Y
               UA=U3
               NBX=N4X
               NBY=N4Y
               UB=U4
               MODE=3
               IF((U4.GT.U0.AND.U1.LE.U0).OR. &
                 (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF

            IF(INDX.EQ.0.OR.INDX.EQ.2) THEN
               XA=X(1)*(NAX-1)+X(2)
               XB=X(1)*(NBX-1)+X(2)
            ELSE
               XA=X(NAX)
               XB=X(NBX)
            ENDIF
            IF(INDX.EQ.0.OR.INDX.EQ.1) THEN
               YA=Y(1)*(NAY-1)+Y(2)
               YB=Y(1)*(NBY-1)+Y(2)
            ELSE
               YA=Y(NAY)
               YB=Y(NBY)
            ENDIF
            J=1
            RT=(U0-UA)/(UB-UA)
            XF(J)=(XB-XA)*RT+XA
            YF(J)=(YB-YA)*RT+YA

            IEL=IE
            NXL=NX
            NYL=NY
            MODEL=MODE
            LINV=.FALSE.
            LEND=.FALSE.

            NXL1=NX
            NYL1=NY
            IF(MODEL.EQ.1) NYL1=NYL-1
            IF(MODEL.EQ.2) NXL1=NXL+1
            IF(MODEL.EQ.3) NYL1=NYL+1
            IF(MODEL.EQ.4) NXL1=NXL-1

            IF(NXL1.GE.1.AND.NXL1.LE.NXM.AND. &
               NYL1.GE.1.AND.NYL1.LE.NYM) THEN
               IE1=NXM*(NYL1-1)+NXL1
            ELSE
               IE1=IE
            ENDIF

  200       IF(MODEL.EQ.1) NYL=NYL-1
            IF(MODEL.EQ.2) NXL=NXL+1
            IF(MODEL.EQ.3) NYL=NYL+1
            IF(MODEL.EQ.4) NXL=NXL-1

            IF(NXL.GE.1.AND.NXL.LE.NXM.AND. &
               NYL.GE.1.AND.NYL.LE.NYM) THEN
               IEL=NXM*(NYL-1)+NXL
               N1X=NXL
               N1Y=NYL
               N2X=NXL+1
               N2Y=NYL
               N3X=NXL+1
               N3Y=NYL+1
               N4X=NXL
               N4Y=NYL+1
               U1=Z(N1X,N1Y)-ZORG
               I2X=N2X
               I2Y=N2Y
               IF(I2X.GT.NXMAX) I2X=1
               IF(I2Y.GT.NYMAX) I2Y=1
               U2=Z(I2X,I2Y)-ZORG
               I3X=N3X
               I3Y=N3Y
               IF(I3X.GT.NXMAX) I3X=1
               IF(I3Y.GT.NYMAX) I3Y=1
               U3=Z(I3X,I3Y)-ZORG
               I4X=N4X
               I4Y=N4Y
               IF(I4X.GT.NXMAX) I4X=1
               IF(I4Y.GT.NYMAX) I4Y=1
               U4=Z(I4X,I4Y)-ZORG

               IF(MODEL.EQ.1) THEN
                  IF((U4.GT.U0.AND.U1.LE.U0).OR. &
                     (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSEIF((U1.GT.U0.AND.U2.LE.U0).OR. &
                         (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR. &
                         (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSE IF(MODEL.EQ.2) THEN
                  IF((U1.GT.U0.AND.U2.LE.U0).OR. &
                     (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR. &
                         (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR. &
                         (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSEIF(MODEL.EQ.3) THEN
                  IF((U2.GT.U0.AND.U3.LE.U0).OR. &
                     (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR. &
                         (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR. &
                         (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSEIF(MODEL.EQ.4) THEN
                  IF((U3.GT.U0.AND.U4.LE.U0).OR. &
                     (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR. &
                         (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSEIF((U1.GT.U0.AND.U2.LE.U0).OR. &
                         (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ENDIF
            ELSE
               IF(LINV) THEN
                  LEND=.TRUE.
               ELSE

                  DO NF=NFMAX,NFMAX-J+1,-1
                     XF(NF)=XF(NF-NFMAX+J)
                     YF(NF)=YF(NF-NFMAX+J)
                  END DO
                  J=NFMAX-J+1

                  IEL=IE
                  NAX=NSAX
                  NAY=NSAY
                  UA=USA
                  NBX=NSBX
                  NBY=NSBY
                  UB=USB
                  NXL=NX
                  NYL=NY
                  MODEL=MODES
                  LINV=.TRUE.
               ENDIF
            ENDIF

            IF(.NOT.LEND) THEN

               IF(INDX.EQ.0.OR.INDX.EQ.2) THEN
                  XA=X(1)*(NAX-1)+X(2)
                  XB=X(1)*(NBX-1)+X(2)
               ELSE
                  XA=X(NAX)
                  XB=X(NBX)
               ENDIF
               IF(INDX.EQ.0.OR.INDX.EQ.1) THEN
                  YA=Y(1)*(NAY-1)+Y(2)
                  YB=Y(1)*(NBY-1)+Y(2)
               ELSE
                  YA=Y(NAY)
                  YB=Y(NBY)
               ENDIF
               IF(.NOT.LINV) THEN
                  IF(J.EQ.NFMAX) THEN
                     CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     CALL LINEPTX(XG,YG,NN,ILNS(K))            
                     J=1
                     XF(1)=XF(NFMAX)
                     YF(1)=YF(NFMAX)
                  ENDIF
                  J=J+1
               ELSE
                  IF(J.EQ.1) THEN
                     CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1, &
      &                           XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     CALL LINEPTX(XG,YG,NN,ILNS(K))
                     J=NFMAX
                     XF(NFMAX)=XF(1)
                     YF(NFMAX)=YF(1)
                  ENDIF
                  J=J-1
               ENDIF
               RT=(U0-UA)/(UB-UA)
               XF(J)=(XB-XA)*RT+XA
               YF(J)=(YB-YA)*RT+YA

               KA(1,IEL)=KA(1,IEL)+1
               IF(IEL.EQ.IE) LEND=.TRUE.
            ENDIF

            IF(.NOT.LEND) GOTO 200

            IF(.NOT.LINV) THEN
               CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
            ELSE
               CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1,XP,YP,NFMAX,NP,ISPL)
            ENDIF
            CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
            CALL LINEPTX(XG,YG,NN,ILNS(K))            

         ENDIF
         END DO
      END DO
      END DO
 1000 CONTINUE
      IF(ISPL.GE.0) THEN
         CALL SETLNW(WS)
         CALL SETRGB(RS,GS,BS)
      ENDIF
      RETURN
  END SUBROUTINE CONTG0X

!     ****** CONTOUR PLOT SLAVE ROUTINE : CONV RT ******

  SUBROUTINE CONTV2X(RA,TA,N,XB,YB,M,NN)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: N,M
    INTEGER,INTENT(OUT):: NN
    REAL(4),INTENT(IN):: RA(N),TA(N)
    REAL(4),INTENT(OUT):: XB(M),YB(M)
    COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
    REAL(4):: DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE
    LOGICAL:: LGF
    COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
    REAL(4):: RMAX,RT,TT,XT,YT
    INTEGER:: J,I,IMAX,K
    REAL(4):: RS,TS,XS,YS,DELR,DELT,RI,TI,XI,YI,X0,Y0,S
    
      RT=RA(1)
      TT=TA(1)
      XT=DX*(RT*COS(TT)-GXS)+PXS
      YT=DY*(RT*SIN(TT)-GYS)+PYS
      J=1
      XB(J)=XT
      YB(J)=YT

      DO I=2,N
         RS=RA(I)
         TS=TA(I)
         XS=DX*(RS*COS(TS)-GXS)+PXS
         YS=DY*(RS*SIN(TS)-GYS)+PYS

         IMAX=INT(SQRT((XS-XT)**2+(YS-YT)**2)*8/RMAX)+1
         DELR=(RS-RT)/IMAX
         DELT= TS-TT
         IF(DELT.GT. 4.0) DELT=DELT-2*3.1415926
         IF(DELT.LT.-4.0) DELT=DELT+2*3.1415926
         DELT=DELT/IMAX
         DO K=1,IMAX
            RI=DELR*K+RT
            TI=DELT*K+TT
            XI=DX*(RI*COS(TI)-GXS)+PXS
            YI=DY*(RI*SIN(TI)-GYS)+PYS
            IF(J.GE.M) THEN
               WRITE(6,*) 'XX GSAF CONTV2 ERROR: BUFFER OVER'
               NN=J
               RETURN
            ENDIF
            J=J+1
            XB(J)=XI
            YB(J)=YI
         END DO
         RT=RS
         TT=TS
         XT=XS
         YT=YS
      END DO
      NN=J

!     ***** slightly expand the region *****

      X0=0.D0
      Y0=0.D0
      DO J=1,NN
         X0=X0+XB(J)
         Y0=Y0+YB(J)
      ENDDO
      X0=X0/NN
      Y0=Y0/NN
      DO J=1,NN
         S=SQRT((XB(J)-X0)**2+(YB(J)-Y0)**2)
         IF(ABS(S).LE.1.E-32) S=1.D0
         XB(J)=XB(J)+0.02*(XB(J)-X0)/S
         YB(J)=YB(J)+0.02*(YB(J)-Y0)/S
      ENDDO
      RETURN
  END SUBROUTINE CONTV2X

!     ****** DRAW LINES WITH PATTERN ******

  SUBROUTINE LINEPTX(XG,YG,N,IPAT)

      DIMENSION XG(N),YG(N)

      CALL MOVEPT(XG(1),YG(1),IPAT)
      DO I=2,N
         CALL DRAWPT(XG(I),YG(I))
      END DO
      RETURN
  END SUBROUTINE LINEPTX

!     ****** SPLINE INTERPOLATION OF 2D LINES ******

  SUBROUTINE GUSP2DX(XH,YH,N,XP,YP,NPM,NP,ISPL)

      PARAMETER(NPA=2001)
      PARAMETER(M=3)
      DIMENSION XH(N),YH(N),XP(NPM),YP(NPM)
      DIMENSION IKN(0-M:NPA+M)

      IF(ISPL.EQ.0) THEN
         DO I=1,N
            XP(I)=XH(I)
            YP(I)=YH(I)
         ENDDO
         NP=N
         RETURN
      ENDIF

      NQ=NPM
      IF(NQ.GT.NPA) NQ=NPA

      IF((ABS(XH(N)-XH(1)).LT.1.E-30).AND. &
         (ABS(YH(N)-YH(1)).LT.1.E-30)) THEN
         IOC=1
      ELSE
         IOC=0
      ENDIF

      IF((ISPL.LE.0).OR.((N.EQ.2).AND.(IOC.EQ.0))) THEN
         DO I=1,N
            XP(I)=XH(I)
            YP(I)=YH(I)
         END DO
         NP=N
      ELSEIF(N.GT.2) THEN
         NP=ISPL*N
         IF(NP.GT.NQ) NP=NQ

         IF(IOC.EQ.0) THEN
            DO I=0-M,0
               IKN(I)=0
            END DO

            DO I=1,N-1
               IKN(I)=I
            END DO

            DO I=N,N+M-1
               IKN(I)=N-1
            END DO
         ELSEIF(IOC.EQ.1) THEN
            DO I=0-M,N+M-1
               IKN(I)=I
            END DO
         ENDIF
         CALL GUCSPLX(N-1,XH,YH,IKN,IOC,NP,XP,YP)
      ENDIF
      RETURN
  END SUBROUTINE GUSP2DX

!     ****** SPLINE ******

  SUBROUTINE GUCSPLX(N,X,Y,IKN,IOC,NP,XP,YP)

      PARAMETER (M=3)
      DIMENSION IKN(0-M:N+M),B(0-M:M,0:M)
      DIMENSION X(0:N),Y(0:N),XP(NP),YP(NP)

      H=(IKN(N)-IKN(0))/REAL(NP-1)

      DO J=1,NP
         TP=IKN(0)+H*(J-1)
         CALL GUBSPLX(TP,ITM,N,IKN,M,B)
         XV=0
         YV=0
         DO I=0,0-M,-1
            IV=ITM+I+2
            IF(IOC.EQ.1) THEN
               IF(IV.LT.0) THEN
                  IV=IV+N
               ELSEIF(IV.GT.N) THEN
                  IV=IV-N
               ENDIF
            ENDIF
            XV=XV+X(IKN(IV))*B(I,M)
            YV=YV+Y(IKN(IV))*B(I,M)
         END DO
         XP(J)=XV
         YP(J)=YV
      END DO
      RETURN
  END SUBROUTINE GUCSPLX

!     ****** B SPLINE ******

  SUBROUTINE GUBSPLX(TP,ITM,N,IKN,M,B)

      DIMENSION IKN(0-M:N+M),B(0-M:M,0:M)

      DO JT=N-1,0,-1
         IF(TP.GE.REAL(IKN(JT))) THEN
            ITM=JT
            GOTO 11
         ENDIF
      END DO
   11 CONTINUE

      DO K=0,M-1
      DO I=-1-K,M-K
         B(I,K)=0.0
      END DO
      END DO

      B(0,0)=1.0
      DO K=1,M
      DO I=0-K,0
         IV=ITM+I
         IF(IKN(IV+K).GT.IKN(IV)) THEN
            B0=(TP-IKN(IV))/(IKN(IV+K)-IKN(IV)+1.E-30) &
                 *B(I,K-1)
         ELSE
            B0=0
         ENDIF

         IF(IKN(IV+K+1).GT.IKN(IV+1)) THEN
            B1=(IKN(IV+K+1)-TP)/(IKN(IV+K+1)-IKN(IV+1)+1.E-30) &
                 *B(I+1,K-1)
         ELSE
            B1=0
         ENDIF

         B(I,K)=B0+B1
      END DO
      END DO
      RETURN
  END SUBROUTINE GUBSPLX
END MODULE dpgsub
