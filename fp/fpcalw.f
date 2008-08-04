C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED)
C
C ************************************************************
C
      SUBROUTINE FPCALW(NSA)
C
      INCLUDE 'fpcomm.inc'
C
C =============  CALCULATION OF DWPP AND DWPT  ===============
C
      FACT=0.5D0
      DO NRDO=1,NRMAX
         NR=NRDO
         DO NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
            DO NP=1,NPMAX+1
               CALL FPSUMW(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP),NR,
     &                     DLHA,DFWA,DECA,DECB,DECC,NSA)
               DWPP(NTH,NP,NR,NSA)  =ABS(COSM(NTH))
     &                               *(DLHA+DFWA+SINM(NTH)**2*DECA)
               DWLHPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DLHA
               DWFWPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DFWA
               DWECPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*SINM(NTH)**2*DECA
               IF(NTH.LE.NTHMAX/2) THEN
                  DWPT(NTH,NP,NR,NSA)  =-SINM(NTH)*(DLHA+DFWA-DECB)
                  DWLHPT(NTH,NP,NR,NSA)=-SINM(NTH)*DLHA
                  DWFWPT(NTH,NP,NR,NSA)=-SINM(NTH)*DFWA
                  DWECPT(NTH,NP,NR,NSA)= SINM(NTH)*DECB
               ELSE
                  DWPT(NTH,NP,NR,NSA)  = SINM(NTH)*(DLHA+DFWA-DECB)
                  DWLHPT(NTH,NP,NR,NSA)= SINM(NTH)*DLHA
                  DWFWPT(NTH,NP,NR,NSA)= SINM(NTH)*DFWA
                  DWECPT(NTH,NP,NR,NSA)=-SINM(NTH)*DECB
               ENDIF
            ENDDO
  101       CONTINUE
         ENDDO
C
         IF(MODELA.EQ.1) THEN
            DO NP=1,NPMAX+1
               DO NTH=ITL(NR)+1,NTHMAX/2
                  DWPP(NTH,NP,NR,NSA)  =(DWPP(NTH,NP,NR,NSA)
     &                            +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA)
     &                            +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPP(NTH,NP,NR,NSA)=(DWLHPP(NTH,NP,NR,NSA)
     &                            +DWLHPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPT(NTH,NP,NR,NSA)=(DWLHPT(NTH,NP,NR,NSA)
     &                            +DWLHPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPP(NTH,NP,NR,NSA)=(DWFWPP(NTH,NP,NR,NSA)
     &                            +DWFWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPT(NTH,NP,NR,NSA)=(DWFWPT(NTH,NP,NR,NSA)
     &                            +DWFWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPP(NTH,NP,NR,NSA)=(DWECPP(NTH,NP,NR,NSA)
     &                            +DWECPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPT(NTH,NP,NR,NSA)=(DWECPT(NTH,NP,NR,NSA)
     &                            +DWECPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
                  DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
                  DWLHPP(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPP(NTH,NP,NR,NSA)
                  DWLHPT(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPT(NTH,NP,NR,NSA)
                  DWFWPP(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPP(NTH,NP,NR,NSA)
                  DWFWPT(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPT(NTH,NP,NR,NSA)
                  DWECPP(NTHMAX-NTH+1,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)
                  DWECPT(NTHMAX-NTH+1,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)
               ENDDO
               DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPP(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPP(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPP(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPP(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWLHPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWLHPP(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWLHPP(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWLHPP(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWLHPP(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWFWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWFWPP(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWFWPP(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWFWPP(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWFWPP(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWECPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWECPP(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWECPP(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWECPP(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWECPP(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
C
               DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPT(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPT(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPT(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPT(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWLHPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWLHPT(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWLHPT(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWLHPT(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWLHPT(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWFWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWFWPT(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWFWPT(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWFWPT(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWFWPT(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWECPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWECPT(ITL(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)-1,NR)
     &                    +DWECPT(ITL(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITL(NR)+1,NR)
     &                    +DWECPT(ITU(NR)-1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)-1,NR)
     &                    +DWECPT(ITU(NR)+1,NP,NR,NSA)
     &                                       /RLAMDA(ITU(NR)+1,NR))
               DWPP(ITU(NR),NP,NR,NSA)  =DWPP(ITL(NR),NP,NR,NSA)
               DWPT(ITU(NR),NP,NR,NSA)  =DWPT(ITL(NR),NP,NR,NSA)
               DWLHPP(ITU(NR),NP,NR,NSA)=DWLHPP(ITL(NR),NP,NR,NSA)
               DWLHPT(ITU(NR),NP,NR,NSA)=DWLHPT(ITL(NR),NP,NR,NSA)
               DWFWPP(ITU(NR),NP,NR,NSA)=DWFWPP(ITL(NR),NP,NR,NSA)
               DWFWPT(ITU(NR),NP,NR,NSA)=DWFWPT(ITL(NR),NP,NR,NSA)
               DWECPP(ITU(NR),NP,NR,NSA)=DWECPP(ITL(NR),NP,NR,NSA)
               DWECPT(ITU(NR),NP,NR,NSA)=DWECPT(ITL(NR),NP,NR,NSA)

            END DO
         ENDIF
      END DO
C
C =============  CALCULATION OF DWTP AND DWTT  ===============
C
      DO NRDO=1,NRMAX
         NR=NRDO
C
         DO NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
               DO NP=1,NPMAX
                  CALL FPSUMW(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
     &                        NR,DLHA,DFWA,DECA,DECB,DECC,NSA)
C
                  IF(NTH.LE.NTHMAX/2) THEN
                     DWTP(NTH,NP,NR,NSA)=-SING(NTH)*(DLHA+DFWA-DECB)
                  ELSE
                     DWTP(NTH,NP,NR,NSA)= SING(NTH)*(DLHA+DFWA-DECB)
                  ENDIF
                  DWTT(NTH,NP,NR,NSA)=(SING(NTH)**2*(DLHA+DFWA)+DECC)
     &                               /ABS(COSG(NTH))
               ENDDO
            ENDIF
         ENDDO
C
         DO NP=1,NPMAX
            CALL FPWAVE(PM(NR),0.D0,NR,EPSR(NR)*RR,0.0D0,
     &                  DLHL,DFWL,DECL,NSA)
            DWTP(NTHMAX/2+1,NP,NR,NSA)=0.D0
            DWTT(NTHMAX/2+1,NP,NR,NSA)=DLHL+DFWL
         ENDDO
C
         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
               DO NP=1,NPMAX
                  DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA)
     &                           +DWTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA)
     &                           +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTP(NTHMAX-NTH+2,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA)
                  DWTT(NTHMAX-NTH+2,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
C
      RETURN
      END
C
C =======================================================
C
      SUBROUTINE FPSUMW(ETA,RSIN,RCOS,P,NR,SUM1,SUM2,SUM3,SUM4,SUM5,NSA)
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
         CALL FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
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
      SUBROUTINE FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
C
      INCLUDE 'fpcomm.inc'
C
      P2=PPARA**2+PPERP**2
      PVPARA=PPARA/SQRT(1.D0+P2*THETA0(NSA))
      RNUDL=0.D0 
      DO NSB=1,NSBMAX
         IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
            RNUDL=RNUD(NR,NSB,NSA)
            RNUFL=RNUF(NR,NSB,NSA)
         ENDIF
      ENDDO
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
               ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-PLH1)/PLH2)**2
               IF(ARG.LT.20.D0) THEN
                  DLHL=DLH*RNUDL+EXP(-ARG)
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
               WPI2=RNFP0(NSA)*1.D20*RNFP(NR,NSA)*AEI*AEE/(AMI*EPS0)
               WFW2=(RFW*1.D6*2*PI)**2
               FACT=RTFP(NR,NSA)+WFW2*AMFP(NSA)*VC*VC
     &                             /(WPI2*RTFP0(NSA)*1.D3)
               IF(ABS(PFW1-2.2D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=DFW*RNUDL*EXP(-ARG)/(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=DFW*RNUFL*EXP(-ARG)/(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSEIF(ABS(PFW1-1.6D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-1.6D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=1.5D0*DFW*RNUDL*EXP(-ARG)
     &                    /(PVPARA**2)
     &                   *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.8D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=0.5D0*DFW*RNUDL*EXP(-ARG)
     &                    /(PVPARA**2)
     &                    *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSE
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-PFW1)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL=DFW*RNUDL*EXP(-ARG)
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
               PARAN=PPARA*SQRT(RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC))
               FN=PEC1*PARAN-SQRT(1.D0+THETA0(NSA)*P2)+W
               DELF=PEC2*PARAN
               ARG2=(FN/DELF)**2
               IF(ARG2.LT.20.0) THEN
                  FACT2=EXP(-ARG2)
                  DECL=DEC*RNUDL*FACT1*FACT2
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
