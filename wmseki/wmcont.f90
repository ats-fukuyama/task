!     $Id: wmcont.f90,v 1.2 2012/08/29 13:45:17 fukuyama Exp $
!
! ***********************************************************
!
!                    SELECTION OF GRAPH
!
! ***********************************************************
!
      MODULE WMCONT

      contains
!---------------------------------------
!
!     ****** CONTOUR PLOT : R-THETA, VARIABLE STEP, PATTERN ******
!
      SUBROUTINE CONTG4X(Z,R,T,NXA,NXMAX,NYMAX,     &
                        ZL,RGB,ILN,WLN,NSTEP,ISPL)
!
!      IMPLICIT NONE
      IMPLICIT LOGICAL(L)
      COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
!
!      EXTERNAL CONTV2X
      real(4),DIMENSION(NXA,NYMAX):: Z
      real(4),dimension(NXMAX):: R
      real(4),dimension(NYMAX):: T
      INTEGER,PARAMETER:: NGLM=30  
      real(4),DIMENSION(NGLM):: ZL,WLN
      integer,dimension(NGLM):: ILN
      real(4),dimension(3,NGLM):: RGB
!
      RMAX=R(NXMAX)
!
      CALL CONTG0X(Z,R,T,NXA,NXMAX,NYMAX, &
                  ZL,RGB,ILN,WLN,NSTEP,ISPL,0,3,CONTV2X)
!
      RETURN
      END SUBROUTINE CONTG4X
!
!     ****** CONTOUR PLOT : COMMON SUB ******
!
      SUBROUTINE CONTG0X(Z,X,Y,NXA,NXMAX,NYMAX, &
                        ZL,RGB,ILN,WLN,NSTEP,   &
                        ISPL,IPRD,INDX,SUBV)
!
      IMPLICIT LOGICAL(L)
      EXTERNAL SUBV
      real(4),DIMENSION(NXA,NYMAX):: Z
      real(4),dimension(NXMAX):: X
      real(4),dimension(NYMAx):: Y
      integer(4),dimension(2,NXMAX*NYMAX):: KA
      INTEGER,PARAMETER:: NGLM=30  
      real(4),DIMENSION(NGLM):: ZL,WLN
      integer,dimension(NGLM):: ILN
      real(4),dimension(3,NGLM):: RGB
      PARAMETER(NH=101)
      real(4),DIMENSION(NH):: ZLS,ILNS,WLNS
      real(4),dimension(3,NH):: RGBS
!      PARAMETER (NFMAX=2000,NGMAX=4000)
      PARAMETER (NFMAX=2000,NGMAX=4000)
      real(4),DIMENSION(NFMAX):: XF,YF,XP,YP
      real(4),DIMENSION(NGMAX):: XG,YG
!
      IF(ISPL.GE.0) THEN
         CALL INQRGB(RS,GS,BS)
         CALL INQLNW(WS)
      ENDIF
!
      KMAX=NSTEP
      IF(KMAX.GT.NH) KMAX=NH
!
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
!
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
!
      NXM=NXMAX-1
      NYM=NYMAX-1
      IF(IPRD.EQ.1.OR.IPRD.EQ.3) NXM=NXMAX
      IF(IPRD.EQ.2.OR.IPRD.EQ.3) NYM=NYMAX
!
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
!
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
!
            IF((U1.GT.U0.AND.U2.LE.U0).OR.   &
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
               IF((U3.GT.U0.AND.U4.LE.U0).OR.  &
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
!
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
!
            IEL=IE
            NXL=NX
            NYL=NY
            MODEL=MODE
            LINV=.FALSE.
            LEND=.FALSE.
!
            NXL1=NX
            NYL1=NY
            IF(MODEL.EQ.1) NYL1=NYL-1
            IF(MODEL.EQ.2) NXL1=NXL+1
            IF(MODEL.EQ.3) NYL1=NYL+1
            IF(MODEL.EQ.4) NXL1=NXL-1
!
            IF(NXL1.GE.1.AND.NXL1.LE.NXM.AND.  &
               NYL1.GE.1.AND.NYL1.LE.NYM) THEN
               IE1=NXM*(NYL1-1)+NXL1
            ELSE
               IE1=IE
            ENDIF
!
  200       IF(MODEL.EQ.1) NYL=NYL-1
            IF(MODEL.EQ.2) NXL=NXL+1
            IF(MODEL.EQ.3) NYL=NYL+1
            IF(MODEL.EQ.4) NXL=NXL-1
!
            IF(NXL.GE.1.AND.NXL.LE.NXM.AND.  &
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
!
               IF(MODEL.EQ.1) THEN
                  IF((U4.GT.U0.AND.U1.LE.U0).OR.   &
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
                  IF((U1.GT.U0.AND.U2.LE.U0).OR.  &
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
!c$$$                     WRITE(6,*) 'TYPE A'
!c$$$                     DO I=1,J
!c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
!c$$$     &                       I,XF(I),YF(I)
!c$$$                     ENDDO
!c$$$                     CALL GUFLSH
!c$$$                     PAUSE
                     CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     IPAT=ILNS(K)
!                     CALL LINEPTX(XG,YG,NN,ILNS(K))
                     CALL LINEPTX(XG,YG,NN,IPAT)
                     J=1
                     XF(1)=XF(NFMAX)
                     YF(1)=YF(NFMAX)
!c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
!c$$$     &                    '-3-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)
                  ENDIF
                  J=J+1
               ELSE
                  IF(J.EQ.1) THEN
!c$$$                     WRITE(6,*) 'TYPE B'
!c$$$                     DO I=J,NFMAX
!c$$$                        WRITE(6,'(I5,1P2E12.4)')
!c$$$     &                       I,XF(I),YF(I)
!c$$$                     ENDDO
!c$$$                     CALL GUFLSH
!c$$$                     PAUSE
                     CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1,  &
                                 XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     IPAT=ILNS(K)
!                     CALL LINEPTX(XG,YG,NN,ILNS(K))
                     CALL LINEPTX(XG,YG,NN,IPAT)
                     J=NFMAX
                     XF(NFMAX)=XF(1)
                     YF(NFMAX)=YF(1)
!c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
!c$$$     &                    '-2-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)
                  ENDIF
                  J=J-1
               ENDIF
               RT=(U0-UA)/(UB-UA)
               XF(J)=(XB-XA)*RT+XA
               YF(J)=(YB-YA)*RT+YA
!c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
!c$$$     &                    '-1-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)

               KA(1,IEL)=KA(1,IEL)+1
               IF(IEL.EQ.IE) LEND=.TRUE.
!               IF(IEL.EQ.IE.AND..NOT.LINV) LEND=.TRUE.
!               IF(IEL.EQ.IE1.AND..NOT.LINV) LEND=.TRUE.
            ENDIF
!
            IF(.NOT.LEND) GOTO 200
!
            IF(.NOT.LINV) THEN
!c$$$                     WRITE(6,*) 'TYPE C'
!c$$$                     DO I=1,J
!c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
!c$$$     &                       I,XF(I),YF(I)
!c$$$                     ENDDO
!c$$$                     CALL GUFLSH
               CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
            ELSE
!c$$$                     WRITE(6,*) 'TYPE D'
!c$$$                     DO I=J,NFMAX
!c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
!c$$$     &                       I,XF(I),YF(I)
!c$$$                     ENDDO
!c$$$                     CALL GUFLSH
               CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1,XP,YP,NFMAX,NP,ISPL)
            ENDIF
            CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
            IPAT=ILNS(K)
!            CALL LINEPTX(XG,YG,NN,ILNS(K))            
            CALL LINEPTX(XG,YG,NN,IPAT)
!
!            IF(.NOT.LINV) THEN
!               WRITE(6,*) J,NP,NN
!            ELSE
!               WRITE(6,*) NFMAX-J+1,NP,NN
!            ENDIF
!            CALL GUFLSH
!
         ENDIF
 1000 CONTINUE
      END DO
      END DO
      END DO
      IF(ISPL.GE.0) THEN
         CALL SETLNW(WS)
         CALL SETRGB(RS,GS,BS)
      ENDIF
      RETURN
      END SUBROUTINE CONTG0X
!
!     ****** CONTOUR PLOT SLAVE ROUTINE : CONV RT ******
!
      SUBROUTINE CONTV2X(RA,TA,N,XB,YB,M,NN)
!
      IMPLICIT LOGICAL(L)
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
      COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
      real(4),DIMENSION(N):: RA,TA
      real(4),dimension(M):: XB,YB
!
      RT=RA(1)
      TT=TA(1)
      XT=DX*(RT*COS(TT)-GXS)+PXS
      YT=DY*(RT*SIN(TT)-GYS)+PYS
      J=1
      XB(J)=XT
      YB(J)=YT
!
      DO  I=2,N
         RS=RA(I)
         TS=TA(I)
         XS=DX*(RS*COS(TS)-GXS)+PXS
         YS=DY*(RS*SIN(TS)-GYS)+PYS
!
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
               WRITE(6,*) 'XX GSAF CONTV2X ERROR: BUFFER OVER'
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
!
!     ***** slightly expand the region *****
!
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
!
!     ****** DRAW LINES WITH PATTERN ******
!
      SUBROUTINE LINEPTX(XG,YG,N,IPAT)
!
      real(4),DIMENSION(N):: XG,YG
      integer:: IPAT
!
      CALL MOVEPT(XG(1),YG(1),IPAT)
      DO I=2,N
         CALL DRAWPT(XG(I),YG(I))
      END DO
      RETURN
      END SUBROUTINE LINEPTX
!
!     ****** SPLINE INTERPOLATION OF 2D LINES ******
!
      SUBROUTINE GUSP2DX(XH,YH,N,XP,YP,NPM,NP,ISPL)
!
      PARAMETER(NPA=2001)
      PARAMETER(M=3)
      real(4),dimension(N):: XH,YH
      real(4),dimension(NPM):: XP,YP
      integer(4),DIMENSION(0-M:NPA+M):: IKN
!
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

      IF((ABS(XH(N)-XH(1)).LT.1.E-30).AND.  &
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
!
         IF(IOC.EQ.0) THEN
            DO I=0-M,0
               IKN(I)=0
            END DO
!
            DO I=1,N-1
               IKN(I)=I
            END DO
!
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
!
!     ****** SPLINE ******
!
      SUBROUTINE GUCSPLX(N,X,Y,IKN,IOC,NP,XP,YP)
!
      PARAMETER (M=3)
      integer(4),DIMENSION(0-M:N+M):: IKN
      real(4),dimension(0-M:M,0:M):: B
      real(4),DIMENSION(0:N):: X,Y
      real(4),dimension(NP):: XP,YP
!
      H=(IKN(N)-IKN(0))/REAL(NP-1)
!
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
!
!     ****** B SPLINE ******
!
      SUBROUTINE GUBSPLX(TP,ITM,N,IKN,M,B)
!
      integer(4),DIMENSION(0-M:N+M):: IKN
      real(4),dimension(0-M:M,0:M):: B
!
      DO JT=N-1,0,-1
         IF(TP.GE.REAL(IKN(JT))) THEN
            ITM=JT
            GOTO 11
         ENDIF
      END DO
   11 CONTINUE
!
      DO K=0,M-1
      DO I=-1-K,M-K
         B(I,K)=0.0
      END DO
      END DO
!
      B(0,0)=1.0
      DO K=1,M
      DO I=0-K,0
         IV=ITM+I
         IF(IKN(IV+K).GT.IKN(IV)) THEN
            B0=(TP-IKN(IV))/(IKN(IV+K)-IKN(IV)+1.E-30)  &
                 * B(I,K-1)
         ELSE
            B0=0
         ENDIF
!
         IF(IKN(IV+K+1).GT.IKN(IV+1)) THEN
            B1=(IKN(IV+K+1)-TP)/(IKN(IV+K+1)-IKN(IV+1)+1.E-30) &
                 *B(I+1,K-1)
         ELSE
            B1=0
         ENDIF
!
         B(I,K)=B0+B1
      END DO
      END DO
      RETURN
      END SUBROUTINE GUBSPLX


      END MODULE WMCONT
