C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED)
C
C ************************************************************
C
      SUBROUTINE FPCALW
C
      INCLUDE 'fpcomm.inc'
C
C =============  CALCULATION OF DWPP AND DWPT  ===============
C
      FACT=0.5D0
      DO 1000 NRDO=1,NRMAX
         NR=NRDO
C
         DO 101 NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
         DO 100 NP=1,NPMAX+1
            CALL FPSUMW(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP),NR,
     &                  DLHA,DFWA,DECA,DECB,DECC)
            DWPP(NTH,NP,NR)  =ABS(COSM(NTH))*(DLHA+DFWA
     &                                       +SINM(NTH)**2*DECA)
            DWLHPP(NTH,NP,NR)=ABS(COSM(NTH))*DLHA
            DWFWPP(NTH,NP,NR)=ABS(COSM(NTH))*DFWA
            DWECPP(NTH,NP,NR)=ABS(COSM(NTH))*SINM(NTH)**2*DECA
            IF(NTH.LE.NTHMAX/2) THEN
               DWPT(NTH,NP,NR)  =-SINM(NTH)*(DLHA+DFWA-DECB)
               DWLHPT(NTH,NP,NR)=-SINM(NTH)*DLHA
               DWFWPT(NTH,NP,NR)=-SINM(NTH)*DFWA
               DWECPT(NTH,NP,NR)= SINM(NTH)*DECB
            ELSE
               DWPT(NTH,NP,NR)  = SINM(NTH)*(DLHA+DFWA-DECB)
               DWLHPT(NTH,NP,NR)= SINM(NTH)*DLHA
               DWFWPT(NTH,NP,NR)= SINM(NTH)*DFWA
               DWECPT(NTH,NP,NR)=-SINM(NTH)*DECB
            ENDIF
  100    CONTINUE
  101    CONTINUE
C
         IF(MODELA.EQ.1) THEN
         DO 200 NP=1,NPMAX+1
C
         DO 190 NTH=ITL(NR)+1,NTHMAX/2
            DWPP(NTH,NP,NR)  =(DWPP(NTH,NP,NR)
     &                        +DWPP(NTHMAX-NTH+1,NP,NR))*FACT
            DWPT(NTH,NP,NR)  =(DWPT(NTH,NP,NR)
     &                        +DWPT(NTHMAX-NTH+1,NP,NR))*FACT
            DWLHPP(NTH,NP,NR)=(DWLHPP(NTH,NP,NR)
     &                        +DWLHPP(NTHMAX-NTH+1,NP,NR))*FACT
            DWLHPT(NTH,NP,NR)=(DWLHPT(NTH,NP,NR)
     &                        +DWLHPT(NTHMAX-NTH+1,NP,NR))*FACT
            DWFWPP(NTH,NP,NR)=(DWFWPP(NTH,NP,NR)
     &                        +DWFWPP(NTHMAX-NTH+1,NP,NR))*FACT
            DWFWPT(NTH,NP,NR)=(DWFWPT(NTH,NP,NR)
     &                        +DWFWPT(NTHMAX-NTH+1,NP,NR))*FACT
            DWECPP(NTH,NP,NR)=(DWECPP(NTH,NP,NR)
     &                        +DWECPP(NTHMAX-NTH+1,NP,NR))*FACT
            DWECPT(NTH,NP,NR)=(DWECPT(NTH,NP,NR)
     &                        +DWECPT(NTHMAX-NTH+1,NP,NR))*FACT
            DWPP(NTHMAX-NTH+1,NP,NR)  =DWPP(NTH,NP,NR)
            DWPT(NTHMAX-NTH+1,NP,NR)  =DWPT(NTH,NP,NR)
            DWLHPP(NTHMAX-NTH+1,NP,NR)=DWLHPP(NTH,NP,NR)
            DWLHPT(NTHMAX-NTH+1,NP,NR)=DWLHPT(NTH,NP,NR)
            DWFWPP(NTHMAX-NTH+1,NP,NR)=DWFWPP(NTH,NP,NR)
            DWFWPT(NTHMAX-NTH+1,NP,NR)=DWFWPT(NTH,NP,NR)
            DWECPP(NTHMAX-NTH+1,NP,NR)=DWECPP(NTH,NP,NR)
            DWECPT(NTHMAX-NTH+1,NP,NR)=DWECPT(NTH,NP,NR)
  190    CONTINUE
         DWPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWLHPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWLHPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWLHPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWLHPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWLHPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWFWPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWFWPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWFWPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWFWPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWFWPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWECPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWECPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWECPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWECPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWECPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
C
         DWPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWLHPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWLHPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWLHPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWLHPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWLHPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWFWPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWFWPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWFWPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWFWPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWFWPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWECPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWECPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWECPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWECPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWECPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWPP(ITU(NR),NP,NR)  =DWPP(ITL(NR),NP,NR)
         DWPT(ITU(NR),NP,NR)  =DWPT(ITL(NR),NP,NR)
         DWLHPP(ITU(NR),NP,NR)=DWLHPP(ITL(NR),NP,NR)
         DWLHPT(ITU(NR),NP,NR)=DWLHPT(ITL(NR),NP,NR)
         DWFWPP(ITU(NR),NP,NR)=DWFWPP(ITL(NR),NP,NR)
         DWFWPT(ITU(NR),NP,NR)=DWFWPT(ITL(NR),NP,NR)
         DWECPP(ITU(NR),NP,NR)=DWECPP(ITL(NR),NP,NR)
         DWECPT(ITU(NR),NP,NR)=DWECPT(ITL(NR),NP,NR)
  200    CONTINUE
         ENDIF
 1000 CONTINUE
C
C =============  CALCULATION OF DWTP AND DWTT  ===============
C
      DO 2000 NRDO=1,NRMAX
         NR=NRDO
C
         DO 1101 NTH=1,NTHMAX+1
            IF(NTH.EQ.NTHMAX/2+1) GOTO 1101
         DO 1100 NP=1,NPMAX
            CALL FPSUMW(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),NR,
     &                  DLHA,DFWA,DECA,DECB,DECC)
C
            IF(NTH.LE.NTHMAX/2) THEN
               DWTP(NTH,NP,NR)=-SING(NTH)*(DLHA+DFWA-DECB)
            ELSE
               DWTP(NTH,NP,NR)= SING(NTH)*(DLHA+DFWA-DECB)
            ENDIF
            DWTT(NTH,NP,NR)=(SING(NTH)**2*(DLHA+DFWA)+DECC)
     &                     /ABS(COSG(NTH))
 1100    CONTINUE
 1101    CONTINUE
C
         DO 1200 NP=1,NPMAX
            CALL FPWAVE(PM(NR),0.D0,NR,EPSR(NR)*RR,0.0D0,
     &                  DLHL,DFWL,DECL)
            DWTP(NTHMAX/2+1,NP,NR)=0.D0
            DWTT(NTHMAX/2+1,NP,NR)=DLHL+DFWL
 1200    CONTINUE
C
         IF(MODELA.EQ.1) THEN
         DO 1300 NTH=ITL(NR)+1,NTHMAX/2
         DO 1300 NP=1,NPMAX
            DWTP(NTH,NP,NR)=(DWTP(NTH,NP,NR)
     &                      +DWTP(NTHMAX-NTH+2,NP,NR))*FACT
            DWTT(NTH,NP,NR)=(DWTT(NTH,NP,NR)
     &                      +DWTT(NTHMAX-NTH+2,NP,NR))*FACT
            DWTP(NTHMAX-NTH+2,NP,NR)=DWTP(NTH,NP,NR)
            DWTT(NTHMAX-NTH+2,NP,NR)=DWTT(NTH,NP,NR)
 1300    CONTINUE
         ENDIF
 2000 CONTINUE
C
      RETURN
      END
C
C =======================================================
C
      SUBROUTINE FPSUMW(ETA,RSIN,RCOS,P,NR,SUM1,SUM2,SUM3,SUM4,SUM5)
C
      INCLUDE 'fpcomm.inc'
C
      DELH=2.D0*ETA/NAVMAX
C
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      SUM4=0.D0
      SUM5=0.D0
C
      DO 100 N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         X=EPSR(NR)*COS(ETAL)*RR
         Y=EPSR(NR)*SIN(ETAL)*RR
C
         IF(MODELA.EQ.0) THEN
            PSI=1.D0
            PSIN=RSIN
            PCOS=ABS(RCOS)
         ELSE
            PSI=(1.D0+EPSR(NR))/(1.D0+X/RR)
            PSIN=SQRT(PSI)*RSIN
            PCOS=SQRT(1.D0-PSI*RSIN**2)
         ENDIF
C
         PPERP=P*PSIN
         IF (RCOS.GT.0.0) THEN
            PPARA =  P*PCOS
         ELSE
            PPARA = -P*PCOS
         END IF
C
         CALL FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL)
C
         SUM1=SUM1+PCOS       *DLHL
         SUM2=SUM2+PCOS       *DFWL
         SUM3=SUM3+PSI/PCOS   *DECL
         SUM4=SUM4+PCOS       *DECL
         SUM5=SUM5+PCOS**3/PSI*DECL
  100 CONTINUE
         SUM1=SUM1*DELH/PI
         SUM2=SUM2*DELH/PI
         SUM3=SUM3*DELH/PI
         SUM4=SUM4*DELH/PI
         SUM5=SUM5*DELH/PI
      RETURN
      END
C
C     ******************************************
C        QUASI-LINEAR WAVE DIFFUSION COEF.
C     ******************************************
C
      SUBROUTINE FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL)
C
      INCLUDE 'fpcomm.inc'
C
      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      P2=PPARA**2+PPERP**2
      PVPARA=PPARA/SQRT(1.D0+P2*THETA0)
C
C        ***** LOWER HTYBRID WAVE *****
C
C        IF (PVPARA.GT.PLH1.AND.PVPARA.LT.PLH2.AND.
C    &       RM(NR).GE.RLH)THEN
C           DLHL=DLH*RNUD(NR,NSFP)
C        ELSE
C           DLHL=0.D0
C        END IF
C
         IF (DLH.GT.0.D0.AND.RM(NR).GE.RLH)THEN
            IF(ABS(PVPARA).GT.1.D-70) THEN
               ARG=((AMFP*VC/(PTFP0*PVPARA)-PLH1)/PLH2)**2
               IF(ARG.LT.20.D0) THEN
                  DLHL=DLH*RNUD(NR,NSFP)*EXP(-ARG)
               ELSE
                  DLHL=0.D0
               ENDIF
            ELSE
               DLHL=0.D0
            ENDIF
         ELSE
            DLHL=0.D0
         END IF
C
C        ***** FAST WAVE *****
C
C        IF (PVPARA.GT.PFW1.AND.PVPARA.LT.PFW2.AND.
C    &       RM(NR).GE.RFW)THEN
C           DFWL=DFW*RNUD(NR,NSFP)/(PPARA**2)
C    &         *(PPERP**2-2.D0*(TE(NR)+0.025D0*551.D0/TE0))**2
C        ELSE
C           DFWL=0.D0
C        END IF
C
         IF (DFW.GT.0.D0.AND.RM(NR).GE.RFW)THEN
            IF(ABS(PVPARA).GT.1.D-70) THEN
               AMI=2*AMP
               AEI=1*AEE
               WPI2=RNFP0*1.D20*RNFP(NR)*AEI*AEE/(AMI*EPS0)
               WFW2=(RFW*1.D6*2*PI)**2
               FACT=RTFP(NR)+WFW2*AMFP*VC*VC/(WPI2*RTFP0*1.D3)
               IF(ABS(PFW1-2.2D0).LT.1.D-30) THEN
                  ARG=((AMFP*VC/(PTFP0*PVPARA)-2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=DFW*RNUD(NR,NSFP)*EXP(-ARG)/(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP*VC/(PTFP0*PVPARA)+2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=DFW*RNUF(NR,NSFP)*EXP(-ARG)/(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSEIF(ABS(PFW1-1.6D0).LT.1.D-30) THEN
                  ARG=((AMFP*VC/(PTFP0*PVPARA)-1.6D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=1.5D0*DFW*RNUD(NR,NSFP)*EXP(-ARG)
     &                    /(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP*VC/(PTFP0*PVPARA)+2.8D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=0.5D0*DFW*RNUD(NR,NSFP)*EXP(-ARG)
     &                    /(PVPARA**2)
     &                    *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSE
                  ARG=((AMFP*VC/(PTFP0*PVPARA)-PFW1)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL=DFW*RNUD(NR,NSFP)*EXP(-ARG)
     &                    /(PVPARA**2)
     &                    *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL=0.D0
                  ENDIF
               ENDIF
            ELSE
               DFWL=0.D0
            ENDIF
         ELSE
            DFWL=0.D0
         END IF
C
C        ***** ELECTRON CYCLOTRON WAVE *****
C
         IF(DEC.GT.0.D0.AND.ABS(PPARA).GT.1.D-70) THEN
            ARG1=(Y/DELYEC)**2
            IF(ARG1.LT.20.0) THEN
               FACT1=EXP(-ARG1)
               W=RFEC/(1.D0+X/RR)
               PARAN=PPARA*SQRT(RTFP0*1.D3*AEE/(AMFP*VC*VC))
               FN=PEC1*PARAN-SQRT(1.D0+THETA0*P2)+W
               DELF=PEC2*PARAN
               ARG2=(FN/DELF)**2
               IF(ARG2.LT.20.0) THEN
                  FACT2=EXP(-ARG2)
                  DECL=DEC*RNUD(NR,NSFP)*FACT1*FACT2
               ELSE
                  DECL=0.D0
               ENDIF
            ELSE
               DECL=0.D0
            ENDIF
         ELSE
            DECL=0.D0
         ENDIF
C
      RETURN
      END
