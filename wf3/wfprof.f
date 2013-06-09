C     $Id$
C
C     ****** PSI ******
C
      SUBROUTINE WFSPSI(X,Y,Z,PSI)
C
      INCLUDE 'wfcomm.inc'
C
      IF(MODELB.EQ.0.OR.
     &   MODELB.EQ.1.OR.
     &   MODELB.EQ.2) THEN
         PSI=(X*X+Y*Y)/(RA*RA)
      ENDIF
      RETURN
      END
C
C     ****** MAGNETIC FIELD PROFILE ******
C
      SUBROUTINE WFSMAG(IN,BABS,AL)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION BLO(3),AL(3)
C
      X=XND(IN)
      Y=YND(IN)
      Z=ZND(IN)
      IF(MODELB.EQ.0) THEN
         BLO(1)=BB
         BLO(2)=0.D0
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.1) THEN
         BLO(1)=0.D0
         BLO(2)=BB
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.2) THEN
         BLO(1)=0.D0
         BLO(2)=0.D0
         BLO(3)=BB
      ENDIF
C
      BABS=0.D0
      DO I=1,3
         BABS=BABS+BLO(I)*BLO(I)
      ENDDO
      BABS=SQRT(BABS)
      DO I=1,3
         AL(I)=BLO(I)/BABS
      ENDDO
C
      RETURN
      END
C
C     ****** DENSITY & TEMPERATURE PROFILE ******
C
      SUBROUTINE WFSDEN(IN,RN,RTPR,RTPP,RZCL)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
C
      IF(MODELP.EQ.0) THEN
         FACT=1.D0
      ELSE
         CALL WFSPSI(XND(IN),YND(IN),ZND(IN),PSI)
         IF(PSI.LT.1.D0) THEN
            IF(MODELP.EQ.1) THEN
               FACT=1.D0
            ELSEIF(MODELP.EQ.2) THEN
               FACT=1.D0-PSI
            ELSEIF(MODELP.EQ.3) THEN
               FACT=EXP(-(ZPMAX-ZND(IN))/0.10D0)
            ELSE
               WRITE(6,*) 'XX WDSDEN: UNKNOWN MODELP = ',MODELP
            ENDIF
         ELSE
            FACT=0.D0
         ENDIF
C
         Z=ZND(IN)
         IF(Z.LT.ZPMIN.OR.Z.GT.ZPMAX) THEN
            FACT=0.D0
         ENDIF
      ENDIF
C
      DO NS=1,NSMAX
         RN(NS)  =(PN(NS)  -PNS(NS))*FACT+PNS(NS)
         RTPR(NS)=PTPR(NS)
         RTPP(NS)=PTPP(NS)
      ENDDO
C
      IF(RN(1).GT.0.D0) THEN
C
         TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
         TI=(RTPR(2)+2.D0*RTPP(2))*1.D3/3.D0
         RLAMEE= 8.0D0+2.3D0*(LOG10(TE)-0.5D0*LOG10(RN(1)))
         RLAMEI= RLAMEE+0.3D0
         RLAMII=12.1D0+2.3D0*(LOG10(TI)-0.5D0*LOG10(RN(1)))
         SN=1.D-20
         PNN0=PPN0/(PTN0*AEE)
C
         DO NS=1,NSMAX
            IF(PZCL(NS).EQ.0) THEN
               IF(NS.EQ.1) THEN
                  TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
                  VTE=SQRT(2.D0*TE*AEE/AME)
                  RNUEE=RN(1)*RLAMEE
     &                 /(1.24D-4*SQRT(TE*1.D-3)**3)
                  RNUEI=0.D0
                  DO NSI=2,NSMAX
                     RNUEI=RNUEI+PZ(NSI)**2*RN(NSI)
                  ENDDO
                  RNUEI=RNUEI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
                  RNUEN=PNN0*SN*0.88D0*VTE
                  RNUE=RNUEE+RNUEI+RNUEN
                  RZCL(NS)=RNUE/(2.D6*PI*RF)
               ELSE
                  TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
                  VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
                  RNUIE=PZ(NS)**2*RN(1)*RLAMEI
     &                 /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
                  RNUII=0.D0
                  DO NSI=2,NSMAX
                     RNUII=RNUII+PZ(NSI)**2*RN(NSI)
                  ENDDO
                  RNUII=RNUII*RLAMII
     &                 /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
                  RNUIN=PNN0*SN*0.88D0*VTI
                  RNUI=RNUIE+RNUII+RNUIN
                  RZCL(NS)=RNUI/(2.D6*PI*RF)
               ENDIF
            ELSE
               RZCL(NS)=PZCL(NS)
            ENDIF
         ENDDO
C     
      ELSE
         DO NS=1,NSMAX
            RZCL(NS)=0.D0
         ENDDO
      ENDIF
C
C         WRITE(6,*) 'ZND= ',ZND(IN)
C         WRITE(6,*) 'RN = ',RN(1),RN(2)
C         WRITE(6,*) 'RT = ',RTPR(1),RTPR(2)
C         WRITE(6,*) 'RZ = ',RZCL(1),RZCL(2)
C         WRITE(6,*) 'E  = ',RNUE,RNUEE,RNUEI,RNUEN
C         WRITE(6,*) 'I  = ',RNUI,RNUIE,RNUII,RNUIN
C         STOP
C
      RETURN
      END
