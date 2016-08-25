!     $Id: fpcalw.f90,v 1.13 2013/01/20 23:24:03 fukuyama Exp $
!
! ************************************************************
!
!            CALCULATION OF DW (BOUNCE AVERAGED)
!
! ************************************************************
!
      MODULE fpcalw

      USE fpcomm

      contains

!------------------------------------------

      SUBROUTINE FP_CALW(NSA)
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NSA, NRDO, NR, NP, NTH
      real(8):: DLHA, DFWA, DECA, DECB, DECC, DLHL, DFWL, DECL
      real(8):: FACT
!
! =============  CALCULATION OF DWPP AND DWPT  ===============
!
      FACT=0.5D0
      DO NRDO=NRSTART,NREND
         NR=NRDO
         DO NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               CALL FPSUMW(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP,NSA),NR, &
                           DLHA,DFWA,DECA,DECB,DECC,NSA)
               DWPP(NTH,NP,NR,NSA)  =ABS(COSM(NTH))  &
                                     *(DLHA+DFWA+SINM(NTH)**2*DECA)
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
!
         IF(MODELA.EQ.1) THEN
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               DO NTH=ITL(NR)+1,NTHMAX/2
                  DWPP(NTH,NP,NR,NSA)  =(DWPP(NTH,NP,NR,NSA) &
                                  +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA) &
                                  +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPP(NTH,NP,NR,NSA)=(DWLHPP(NTH,NP,NR,NSA) &
                                  +DWLHPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPT(NTH,NP,NR,NSA)=(DWLHPT(NTH,NP,NR,NSA) &
                                  +DWLHPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPP(NTH,NP,NR,NSA)=(DWFWPP(NTH,NP,NR,NSA) &
                                  +DWFWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPT(NTH,NP,NR,NSA)=(DWFWPT(NTH,NP,NR,NSA) &
                                  +DWFWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPP(NTH,NP,NR,NSA)=(DWECPP(NTH,NP,NR,NSA) &
                                  +DWECPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPT(NTH,NP,NR,NSA)=(DWECPT(NTH,NP,NR,NSA) &
                                  +DWECPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
                  DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
                  DWLHPP(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPP(NTH,NP,NR,NSA)
                  DWLHPT(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPT(NTH,NP,NR,NSA)
                  DWFWPP(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPP(NTH,NP,NR,NSA)
                  DWFWPT(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPT(NTH,NP,NR,NSA)
                  DWECPP(NTHMAX-NTH+1,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)
                  DWECPT(NTHMAX-NTH+1,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)
               ENDDO
               DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0     &
                        *( DWPP(ITL(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWPP(ITL(NR)+1,NP,NR,NSA)               &
                                            /RLAMDA(ITL(NR)+1,NR)  &
                          +DWPP(ITU(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWPP(ITU(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)+1,NR)) 
               DWLHPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWLHPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWLHPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWLHPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWLHPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWFWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWFWPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWFWPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWFWPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWFWPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWECPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWECPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWECPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWECPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWECPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
!
               DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0     &
                        *( DWPT(ITL(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWPT(ITL(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWPT(ITU(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWPT(ITU(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWLHPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWLHPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWLHPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWLHPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWLHPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWFWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWFWPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWFWPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWFWPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWFWPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWECPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWECPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWECPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWECPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWECPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
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
!
! =============  CALCULATION OF DWTP AND DWTT  ===============
!
      DO NRDO=NRSTART,NREND
         NR=NRDO
!
         DO NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  CALL FPSUMW(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP,NSA), &
                              NR,DLHA,DFWA,DECA,DECB,DECC,NSA)
!
                  IF(NTH.LE.NTHMAX/2) THEN
                     DWTP(NTH,NP,NR,NSA)=-SING(NTH)*(DLHA+DFWA-DECB)
                  ELSE
                     DWTP(NTH,NP,NR,NSA)= SING(NTH)*(DLHA+DFWA-DECB)
                  ENDIF
                  DWTT(NTH,NP,NR,NSA)=(SING(NTH)**2*(DLHA+DFWA)+DECC) &
                                     /ABS(COSG(NTH))
               ENDDO
            ENDIF
         ENDDO
!
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            CALL FPWAVE(PM(NP,NSA),0.D0,NR,RSRHON(RM(NR)),0.0D0, &
                        DLHL,DFWL,DECL,NSA)
            DWTP(NTHMAX/2+1,NP,NR,NSA)=0.D0
            DWTT(NTHMAX/2+1,NP,NR,NSA)=DLHL+DFWL
         ENDDO
!
         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA) &
                                 -DWTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA) &
                                 +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTP(NTHMAX-NTH+2,NP,NR,NSA)=-DWTP(NTH,NP,NR,NSA)
                  DWTT(NTHMAX-NTH+2,NP,NR,NSA)= DWTT(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE FP_CALW
!
! =======================================================
!
      SUBROUTINE FPSUMW(ETA,RSIN,RCOS,P,NR,SUM1,SUM2,SUM3,SUM4,SUM5,NSA)
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NR, NSA, N
      real(8):: ETA, RSIN, RCOS, P, SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: DELH, ETAL, X, Y, PSI, PSIN, PCOS, PPERP, PPARA
      real(8):: DLHL, DFWL, DECL, XM
!
      DELH=2.D0*ETA/NAVMAX
!
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      SUM4=0.D0
      SUM5=0.D0
!
      DO N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         X=RSRHON(RM(NR))*COS(ETAL)
         Y=RSRHON(RM(NR))*SIN(ETAL)
         XM=EPSRM(NR)*COS(ETAL)*RR
!
         IF(MODELA.EQ.0) THEN
            PSI=1.D0
            PSIN=RSIN
!            PCOS=ABS(RCOS)
            PCOS=RCOS
         ELSE
            PSI=(1.D0+EPSRM(NR))/(1.D0+XM/RR)
            PSIN=SQRT(PSI)*RSIN
!            PCOS=SQRT(1.D0-PSI*RSIN**2)
!            write(6,'(A,I5,1P5E12.4)') 'NR:',NR,&
!                 EPSRM(NR),X/RR,PSI,RSIN,1.D0-PSI*RSIN**2
            IF (RCOS.GT.0.0D0) THEN
               PCOS= SQRT(1.D0-PSI*RSIN**2)
            ELSE
               PCOS=-SQRT(1.D0-PSI*RSIN**2)
            END IF
         ENDIF
!
         PPERP=P*PSIN
!         IF (RCOS.GT.0.0) THEN
            PPARA =  P*PCOS
!         ELSE
!            PPARA = -P*PCOS
!         END IF
!
         CALL FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
!
         SUM1=SUM1+PCOS       *DLHL
         SUM2=SUM2+PCOS       *DFWL
!         SUM3=SUM3+PSI/PCOS   *DECL
!         SUM4=SUM4+PCOS       *DECL
!         SUM5=SUM5+PCOS**3/PSI*DECL
         SUM3=SUM3+DECL*RCOS/PCOS
         SUM4=SUM4+DECL/SQRT(PSI)
         SUM5=SUM5+DECL*PCOS/RCOS/PSI
!!! unverified for LH, FW
      END DO
      SUM1=SUM1*DELH
      SUM2=SUM2*DELH
      SUM3=SUM3*DELH
      SUM4=SUM4*DELH
      SUM5=SUM5*DELH

      RETURN
      END SUBROUTINE FPSUMW
!
!     ******************************************
!        QUASI-LINEAR WAVE DIFFUSION COEF.
!     ******************************************
!
      SUBROUTINE FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB
      real(8):: PPARA, PPERP, X, Y, DLHL, DFWL, DECL
      real(8):: P2, PVPARA, RNUDL, RNUFL, AMI, AEI, WPI2, FACT, FACT2
      real(8):: DFWL1, DFWL2, ARG, ARG1, FACT1, W, PARAN, FN, DELF, ARG2
      real(8):: WFW2, ARG3, FACT3
!
      P2=PPARA**2+PPERP**2
      PVPARA=PPARA/SQRT(1.D0+P2*THETA0(NSA))
      RNUDL=0.D0 
      DO NSB=1,NSBMAX
         IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
            RNUDL=RNUD(NR,NSB,NSA)
            RNUFL=RNUF(NR,NSB,NSA)
         ENDIF
      ENDDO
!
!        ***** LOWER HTYBRID WAVE *****
!
!        IF (PVPARA.GT.PLH1.AND.PVPARA.LT.PLH2.AND.
!    &       RM(NR).GE.RLH)THEN
!           DLHL=DLH*RNUD(NR,NSA)
!        ELSE
!           DLHL=0.D0
!        END IF
!
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
!
!        ***** FAST WAVE *****
!
!        IF (PVPARA.GT.PFW1.AND.PVPARA.LT.PFW2.AND.
!    &       RM(NR).GE.RFW)THEN
!           DFWL=DFW*RNUD(NR,NSA)/(PPARA**2)
!    &         *(PPERP**2-2.D0*(TE(NR)+0.025D0*551.D0/TE0))**2
!        ELSE
!           DFWL=0.D0
!        END IF
!
         IF (DFW.GT.0.D0.AND.RM(NR).GE.RFW)THEN
            IF(ABS(PVPARA).GT.1.D-70) THEN
               AMI=2*AMP
               AEI=1*AEE
               WPI2=RNFP0(NSA)*1.D20*RNFP(NR,NSA)*AEI*AEE/(AMI*EPS0)
               WFW2=(RFW*1.D6*2*PI)**2
               FACT=RTFP(NR,NSA)+WFW2*AMFP(NSA)*VC*VC &
                                   /(WPI2*RTFP0(NSA)*1.D3)
               IF(ABS(PFW1-2.2D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=DFW*RNUDL*EXP(-ARG)/(PVPARA**2) &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=DFW*RNUFL*EXP(-ARG)/(PVPARA**2) &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSEIF(ABS(PFW1-1.6D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-1.6D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=1.5D0*DFW*RNUDL*EXP(-ARG)   &
                          /(PVPARA**2)                 &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.8D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=0.5D0*DFW*RNUDL*EXP(-ARG)   &
                          /(PVPARA**2)                 &
                          *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSE
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-PFW1)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL=DFW*RNUDL*EXP(-ARG)          &
                          /(PVPARA**2)                 &
                          *(PPERP**2-2.D0*FACT)**2
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
!
!        ***** ELECTRON CYCLOTRON WAVE *****
!
         IF(DEC.GT.0.D0.AND.ABS(PPARA).GT.1.D-70) THEN
            ARG1=(Y/DELYEC)**2
            IF(ARG1.LT.20.0) THEN
               FACT1=EXP(-ARG1)
            ELSE
               FACT1=0.D0
            END IF

            W=RFEC/(1.D0+X/RR)
            PARAN=PPARA*SQRT(RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC))
            FN=PEC1*PARAN-SQRT(1.D0+THETA0(NSA)*P2)+W
            DELF=PEC2*PARAN
            ARG2=(FN/DELF)**2
            IF(ARG2.LT.20.0) THEN
               FACT2=EXP(-ARG2)
            ELSE
               FACT2=0.D0
            END IF

            ARG3=(X-PEC3)/PEC4
            IF(ARG3.GT.20.D0) THEN
               FACT3=1.D0
            ELSE IF(ARG3.LT.-20.D0) THEN
               FACT3=0.D0
            ELSE
               FACT3=0.5D0*(1.D0+TANH(ARG3))
            ENDIF
              
            DECL=DEC*RNUDL*FACT1*FACT2*FACT3
         ELSE
            DECL=0.D0
         ENDIF
!
      RETURN
      END SUBROUTINE FPWAVE

!
!***********************************************************************
!     Calculate PSIN, PCOS, PSI
!***********************************************************************
!
      SUBROUTINE FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI,NSA)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,NSA
      REAL(8),INTENT(IN):: ETAL,RSIN,RCOS
      REAL(8),INTENT(OUT):: PSIN,PCOS,PSI
      REAL(8):: ARG

      IF(MODELA.EQ.0) THEN
         PSI=1.D0
         PSIN=RSIN
         PCOS=RCOS
      ELSE
         PSI=(1.D0+EPSRM(NR))/(1.D0+EPSRM(NR)*COS(ETAL))
         PSIN=SQRT(PSI)*RSIN
         ARG=1.D0-PSI*RSIN**2
         IF(ARG.LT.0.D0) ARG=0.D0
         IF (RCOS.GT.0.0D0) THEN
            PCOS= SQRT(ARG)
         ELSE
            PCOS=-SQRT(ARG)
         END IF
      ENDIF
      RETURN
      END SUBROUTINE FPDWRP
!---------------------------------
      END MODULE fpcalw
      
