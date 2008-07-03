C     $Id$
C
C *************************
C     SAVE PROFILE DATA
C *************************
C
      SUBROUTINE FPSPRF
C
      INCLUDE 'fpcomm.inc'
C
      IF(NTG1.LT.NTG1M) NTG1=NTG1+1
C
      IF(ISAVE.EQ.0) CALL FPSSUB
C
      RTG(NTG1)=TIMEFP
C
      NTG4 = INT( (NTG1+NSFPMA-NSFPMI-1)/(NSFPMA-NSFPMI+1)) +1
      if(NSFPMA.eq.NSFPMI) NTG4 = NTG1
      DO 1000 NR=1,NRMAX
         RNT(NR,NTG1) = RNS(NR)
         RJT(NR,NTG1) = RJS(NR)
         RWT(NR,NTG1) = RWS(NR)
         RPCT(NR,NTG1)= RPCS(NR)
         RPWT(NR,NTG1)= RPWS(NR)
         RPET(NR,NTG1)= RPES(NR)
         RLHT(NR,NTG1)= RLHS(NR)
         RFWT(NR,NTG1)= RFWS(NR)
         RECT(NR,NTG1)= RECS(NR)
         Do NS=1,NSMAX
            RNT2(NR,NTG4,NS) = RNS2(NR,NS)
            RJT2(NR,NTG4,NS) = RJS2(NR,NS)
            RPCT2(NR,NTG4,NS)= RPCS2(NR,NS)
            RWT2(NR,NTG4,NS) = RWS2(NR,NS)
            RTT2(NR,NTG4,NS) = RWS2(NR,NS)*1.D6
     &           /(1.5D0*RNS2(NR,NS)*1.D20*AEE*1.D3)

         END DO
C
         RTT(NR,NTG1) = RWS(NR)*1.D6/(1.5D0*RNS(NR)*1.D20*AEE*1.D3)
         RET(NR,NTG1) = E1(NR)
         RS=RSRHON(RM(NR))
         RQT(NR,NTG1) = RS*BB*2.D0/(RR*(BP(NR)+BP(NR+1)))

 1000 CONTINUE
C
      RETURN
      END
C
C *************************
C     SAVE GLOBAL DATA
C *************************
C
      SUBROUTINE FPSGLB
C
      INCLUDE 'fpcomm.inc'
C
      IF(NTG2.LT.NTG2M) NTG2=NTG2+1
      NTG3 = INT( (NTG2+NSFPMA-NSFPMI-1)/(NSFPMA-NSFPMI+1)) +1
      if(NSFPMA.eq.NSFPMI) NTG3 = NTG2
C
      IF(ISAVE.EQ.0) CALL FPSSUB
C
      PTG(NTG2)=TIMEFP
C
      PNT(NTG2)=0.D0
      PIT(NTG2)=0.D0
      PWT(NTG2)=0.D0
      PPCT(NTG2)=0.D0
      PPWT(NTG2)=0.D0
      PPET(NTG2)=0.D0
      PLHT(NTG2)=0.D0
      PFWT(NTG2)=0.D0
      PECT(NTG2)=0.D0
      DO NS=1,NSMAX
         PPCT2(NTG3,NS)= 0.D0
      END DO
      IF(NTG2.eq.1)THEN
         DO NS=1,NSMAX
            PNT2(1,NS) = 0.D0
            PWT2(1,NS) = 0.D0
            PIT2(1,NS) = 0.D0
         END DO
      END IF

      PNT2(NTG3,NSFP) = 0.D0
      PWT2(NTG3,NSFP) = 0.D0
      PIT2(NTG3,NSFP) = 0.D0

C
      DO 1000 NR=1,NRMAX
         RHOL=RM(NR)
         RHOL1=RG(NR)
         RHOL2=RG(NR+1)
         FACT=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1))
         PNT(NTG2) =PNT(NTG2) +RNS(NR)*FACT
         PIT(NTG2) =PIT(NTG2) +RJS(NR)*FACT
         PWT(NTG2) =PWT(NTG2) +RWS(NR)*FACT
         PPCT(NTG2)=PPCT(NTG2)+RPCS(NR)*FACT
         PPWT(NTG2)=PPWT(NTG2)+RPWS(NR)*FACT
         PPET(NTG2)=PPET(NTG2)+RPES(NR)*FACT
         PLHT(NTG2)=PLHT(NTG2)+RLHS(NR)*FACT
         PFWT(NTG2)=PFWT(NTG2)+RFWS(NR)*FACT
         PECT(NTG2)=PECT(NTG2)+RECS(NR)*FACT
         DO NS=1,NSMAX
            PPCT2(NTG2,NS)=PPCT2(NTG2,NS)+RPCS2(NR,NS)*FACT
            IF(NTG2.eq.1)THEN
               PNT2(1,NS) =PNT2(1,NS) +RNS2(NR,NS)*FACT
               PWT2(1,NS) =PWT2(1,NS) +RWS2(NR,NS)*FACT
               PIT2(1,NS) =PIT2(1,NS) +RJS2(NR,NS)*FACT
            END IF
         END DO
         IF(NTG2.ne.1)THEN
            PNT2(NTG3,NSFP) =PNT2(NTG3,NSFP) +RNS2(NR,NSFP)*FACT
            PWT2(NTG3,NSFP) =PWT2(NTG3,NSFP) +RWS2(NR,NSFP)*FACT
            PIT2(NTG3,NSFP) =PIT2(NTG3,NSFP) +RJS2(NR,NSFP)*FACT
         END IF

 1000 CONTINUE
C
      PNT(NTG2) =PNT(NTG2) *2.0*PI*RR
      PIT(NTG2) =PIT(NTG2) 
      PWT(NTG2) =PWT(NTG2) *2.D0*PI*RR
      PPCT(NTG2)=PPCT(NTG2)*2.D0*PI*RR
      PPWT(NTG2)=PPWT(NTG2)*2.D0*PI*RR
      PPET(NTG2)=PPET(NTG2)*2.D0*PI*RR
      PLHT(NTG2)=PLHT(NTG2)*2.D0*PI*RR
      PFWT(NTG2)=PFWT(NTG2)*2.D0*PI*RR
      PECT(NTG2)=PECT(NTG2)*2.D0*PI*RR
      PTT(NTG2) =PWT(NTG2)*1.D6/(1.5D0*PNT(NTG2)*1.D20*AEE*1.D3)
      RS=RSRHON(RM(1))
      PQT(NTG2) =RS*BB*2.D0/(RR*(BP(1)+BP(2)))
      PET(NTG2) =E1(NRMAX)
      
      Do NS=1,NSMAX
         PPCT2(NTG2,NS)=PPCT2(NTG2,NS)*2.D0*PI*RR
         IF(NTG3.eq.1)THEN
            PNT2(1,NS) =PNT2(1,NS) *2.D0*PI*RR
            PWT2(1,NS) =PWT2(1,NS) *2.D0*PI*RR
            PTT2(1,NS) =
     &           PWT2(1,NS)*1.D6/(1.5D0*PNT2(1,NS)*1.D20*AEE*1.D3)
            PNT2(1,NS)=PNT2(1,NS)/(2.D0*PI*RR
     &           *2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)))
         END IF
      END DO

      IF(NTG2.ne.1)THEN
         PNT2(NTG3,NSFP) =PNT2(NTG3,NSFP) *2.D0*PI*RR
         PWT2(NTG3,NSFP) =PWT2(NTG3,NSFP) *2.D0*PI*RR
         PTT2(NTG3,NSFP) =
     &       PWT2(NTG3,NSFP)*1.D6/(1.5D0*PNT2(NTG3,NSFP)*1.D20*AEE*1.D3)
         PNT2(NTG3,NSFP)=PNT2(NTG3,NSFP)/(2.D0*PI*RR
     &       *2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)))
      END IF

C     density
      PNT(NTG2)=PNT(NTG2)/(2.D0*PI*RR
     &     *2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)))
c      PNT(NTG2)=RNS(1)
c      PNT2(NTG3,NSFP)=RNS2(1,NSFP)

      RETURN
      END
C
C *************************
C     SAVE DATA ROUTINE
C *************************
C
      SUBROUTINE FPSSUB
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION RSUM10(NSMAX),RSUM33(NSMAX),RSUM11(NSMAX)
     &     ,RSUM22(NSMAX)
      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
C
      DO 1000 NR=1,NRMAX
         RSUM1=0.D0
         RSUM2=0.D0
         RSUM3=0.D0
         RSUM4=0.D0
         RSUM5=0.D0
         RSUM6=0.D0
         RSUM7=0.D0
         RSUM8=0.D0
         RSUM9=0.D0
         Do NS=1,NSMAX
            RSUM10(NS)=0.D0
            RSUM11(NS)=0.D0
            RSUM22(NS)=0.D0
            RSUM33(NS)=0.D0
         END DO
C
         DO 100 NP=1,NPMAX
         DO 100 NTH=1,NTHMAX
            RSUM1 = RSUM1+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &       *F(NTH,NP,NR)
            DO NS=1,NSMAX
               RSUM11(NS) = 
     &          RSUM11(NS)+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &          *FNS2(NTH,NP,NR,NS)
            END DO
  100    CONTINUE
C
         IF(MODELR.EQ.0) THEN
            DO 200 NP=1,NPMAX
            DO 200 NTH=1,NTHMAX
               RSUM2 = RSUM2+VOL(NTH,NP)*F(NTH,NP,NR)*PM(NP)*COSM(NTH)
               RSUM3 = RSUM3+VOL(NTH,NP)*F(NTH,NP,NR)
c*RLAMDA(NTH,NR)
     &                       *0.5D0*PM(NP)**2
               DO NS=1,NSMAX
                  RSUM22(NS) = RSUM22(NS)
     &                 +VOL(NTH,NP)*FNS2(NTH,NP,NR,NS)*PM(NP)*COSM(NTH)
                  RSUM33(NS) = RSUM33(NS)+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &                 *FNS2(NTH,NP,NR,NS)*0.5D0*PM(NP)**2
               END DO
  200       CONTINUE
         ELSE
            DO 300 NP=1,NPMAX
               PV=SQRT(1.D0+THETA0*PM(NP)**2)
            DO 300 NTH=1,NTHMAX
               RSUM2 = RSUM2+VOL(NTH,NP)*F(NTH,NP,NR)*PM(NP)*COSM(NTH)
     &                       /PV
               RSUM3 = RSUM3+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &          *F(NTH,NP,NR)
     &                       *(PV-1.D0)/THETA0
               DO NS=1,NSMAX
                  IF(NTG2.eq.0)THEN
                     AMFD=PA(NS)*AMP
                     RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
                     PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)
                     THETB0=RTFD0*1.D3*AEE/(AMFD*VC*VC)
                     PV2=SQRT(1.D0+THETB0*PM(NP)**2)

                     RSUM22(NS) = RSUM22(NS) + + VOL(NTH,NP)
     &                    *FNS2(NTH,NP,NR,NS)*PM(NP)*COSM(NTH)/PV2
                     RSUM33(NS) = RSUM33(NS)+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &                    *FNS(NTH,NP,NR,NS)*(PV2-1.D0)/THETB0
                  ELSE
                     RSUM22(NS) = RSUM22(NS) + VOL(NTH,NP)
     &                    *FNS2(NTH,NP,NR,NS)*PM(NP)*COSM(NTH)/PV
                     RSUM33(NS) = RSUM33(NS)+VOL(NTH,NP)
c*RLAMDA(NTH,NR)
     &                    *F(NTH,NP,NR)*(PV-1.D0)/THETA0                     
                  END IF
               END DO
  300       CONTINUE
         ENDIF

C
c         NP=1
c         NTH=NTHMAX
c         write(*,645)NP, NTH, RSUM33(1)
c     &           ,PG(NP),F(NTH,NP,NR),FNS(NTH,NP,NR,NSFP)
         DO 400 NP=2,NPMAX
            PV=SQRT(1.D0+THETA0*PG(NP)**2)
         DO 400 NTH=1,NTHMAX
            WPL=WEIGHP(NTH  ,NP,NR)
            IF(NTH.EQ.1) THEN
               WPM=0.D0
            ELSE
               WPM=WEIGHP(NTH-1,NP,NR)
            ENDIF
            IF(NTH.EQ.NTHMAX) THEN
               WPP=0.D0
            ELSE
               WPP=WEIGHP(NTH+1,NP,NR)
            ENDIF
            DFP=    PG(NP)*RLAMDA(NTH,NR)/DELP
     &                   *(F(NTH,NP,NR)-F(NTH,NP-1,NR))
            IF(NTH.EQ.1) THEN
               DFT=    1.D0/(2.D0*DELTH)
     &                   *( RLAMDA(NTH+1,NR)
     &                      *((1.D0-WPP)*F(NTH+1,NP  ,NR)
     &                             +WPP *F(NTH+1,NP-1,NR)))
            ELSE IF(NTH.EQ.NTHMAX) THEN
               DFT=    1.D0/(2.D0*DELTH)
     &                   *(-RLAMDA(NTH-1,NR)
     &                      *((1.D0-WPM)*F(NTH-1,NP  ,NR)
     &                             +WPM *F(NTH-1,NP-1,NR)))
            ELSE
               DFT=    1.D0/(2.D0*DELTH)
     &                   *( RLAMDA(NTH+1,NR)
     &                      *((1.D0-WPP)*F(NTH+1,NP  ,NR)
     &                             +WPP *F(NTH+1,NP-1,NR))
     &                     -RLAMDA(NTH-1,NR)
     &                      *((1.D0-WPM)*F(NTH-1,NP  ,NR)
     &                             +WPM *F(NTH-1,NP-1,NR)))
            ENDIF
            FFP=    PG(NP)*RLAMDA(NTH,NR)
     &                   *((1.D0-WPL)*F(NTH  ,NP  ,NR)
     &                          +WPL *F(NTH  ,NP-1,NR))  
            RSUM4 = RSUM4+PG(NP)**2*SINM(NTH)/PV
     &              *(DCPP(NTH,NP,NR)*DFP
     &               +DCPT(NTH,NP,NR)*DFT
     &               -FCPP(NTH,NP,NR)*FFP)
            RSUM5 = RSUM5+PG(NP)**2*SINM(NTH)/PV
     &              *(DWPP(NTH,NP,NR)*DFP
     &               +DWPT(NTH,NP,NR)*DFT)
            RSUM6 = RSUM6-PG(NP)**2*SINM(NTH)/PV
     &              *(FEPP(NTH,NP,NR)*FFP)
            RSUM7 = RSUM7+PG(NP)**2*SINM(NTH)/PV
     &              *(DWLHPP(NTH,NP,NR)*DFP
     &               +DWLHPT(NTH,NP,NR)*DFT)
            RSUM8 = RSUM8+PG(NP)**2*SINM(NTH)/PV
     &              *(DWFWPP(NTH,NP,NR)*DFP
     &               +DWFWPT(NTH,NP,NR)*DFT)
            RSUM9 = RSUM9+PG(NP)**2*SINM(NTH)/PV
     &              *(DWECPP(NTH,NP,NR)*DFP
     &               +DWECPT(NTH,NP,NR)*DFT)
            DO NS=1,NSMAX
               RSUM10(NS)=RSUM10(NS)+PG(NP)**2*SINM(NTH)/PV
     &              *(DCPP2(NTH,NP,NR,NS)*DFP
     &               +DCPT2(NTH,NP,NR,NS)*DFT
     &               -FCPP2(NTH,NP,NR,NS)*FFP)
            END DO

 645     FORMAT(2I3,5E16.8)
  400    CONTINUE
         FACT=RNFP0*1.D20
         RNS(NR) = RSUM1*FACT                   *1.D-20
         RJS(NR) = RSUM2*FACT*AEFP*PTFP0/AMFP*1.D-6
         DO NS=1,NSMAX
            RNS2(NR,NS) = RSUM11(NS) * PN(NS)
            RJS2(NR,NS) = RSUM22(NS)*FACT*AEFP*PTFP0/AMFP*1.D-6
         END DO
         FACT=RNFP0*1.D20*PTFP0**2/AMFP 
         RWS(NR) = RSUM3*FACT               *1.D-6

         RPCS(NR)=-RSUM4*FACT*2.D0*PI*DELP*DELTH *1.D-6
         RPWS(NR)=-RSUM5*FACT*2.D0*PI*DELP*DELTH *1.D-6 
         RPES(NR)=-RSUM6*FACT*2.D0*PI*DELP*DELTH *1.D-6
         RLHS(NR)=-RSUM7*FACT*2.D0*PI*DELP*DELTH *1.D-6
         RFWS(NR)=-RSUM8*FACT*2.D0*PI*DELP*DELTH *1.D-6
         RECS(NR)=-RSUM9*FACT*2.D0*PI*DELP*DELTH

         Do NS=1,NSMAX
            AMFD=PA(NS)*AMP
            RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
            PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)

            RWS2(NR,NS) = RSUM33(NS)*PN(NS)*1.D20*PTFD0**2/AMFD*1.D-6
            RPCS2(NR,NS)=-RSUM10(NS)*FACT*2.D0*PI*DELP*DELTH *1.D-6
         END DO
C
C         IF(MODELA.EQ.1) THEN
C            RJS(NR)=RJS(NR)*RLAMDB(NR)
C         ENDIF
C
 1000 CONTINUE

      ISAVE=1
      RETURN
      END
C
C ***********************************************************
C
C                         RESULT
C
C ***********************************************************
C
      SUBROUTINE FPWRIT
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
C
      WRITE(6,101) TIMEFP*1000
      WRITE(6,102) PNT(NTG2),PIT(NTG2),PWT(NTG2),PTT(NTG2)

      IF(NTG2.eq.1)THEN
         DO NS=1,NSMAX
            write(6,1467)  NS, PTT2(NTG2,NS)
         END DO
         WRITE(6,*) "NTG2=",NTG2,"     NTG3=",NTG3,"  NTG1=" ,NTG1
      END IF

      if(NTG2.ne.1)then
c-----check of conductivity--------
      FACT1 = 
     &2.D0*PI*RSRHON(RM(NRMAX))*(RSRHON(RG(2))-RSRHON(RG(NRMAX)))
      FACT2 = 2.D0*PI*RR
      rntv = 1.6D0/1.D19*1.5D0*
     &FACT1*FACT2*(PTT(NTG2)-PTT(NTG2-1))*PNT(NTG2)*1.D3*1.D20/DELT

c      write(6,*) "Pw+Pc",(PPWT(NTG2)+PPCT(NTG2))*1.D6,"n*delta T*V",rntv
c      write(6,*) "c*epsilon*RabsE*S",RABSE**2*FACT1*8.8D0*VC/1.D12
c     & ,"Pw[W]",PPWT(NTG2)*1.D6
c      write(6,*) "Pc+PE+Pw",PPCT(NTG2)+PPET(NTG2)+PPWT(NTG2),
c     &     "(Pc+Pe+Pw)/Pc",(PPCT(NTG2)+PPET(NTG2)+PPWT(NTG2))/PPCT(NTG2)
c      write(6,*) "Pc+-Pc-/Pc",(PPCT(NTG2)-PPCT(NTG2-1))/PPCT(NTG2)
c     Spitzer
c      resist = 1.D0/(RNFP0*RNFP(1)*1.D20*AEFP**2/AMFP
c     &     /RNUD(1,NSFP)*RNFD(1,NSFP)/RNFP0)*SQRT(2.D0)
c      write(6,*) "J/E*eta*1.D6", PIT(NTG2)/E0*1.D6/FACT1*resist
c     & ,"THETA0", (PTFP0/(AMFP*VC))**2
c----end of conductivity check---------
c----check of beam power-------
      IF(NSBM.ne.0)THEN
c         su1=0.D0
c         if(NTG2.eq.1)then
c            sumt1=0.D0
c            sumt2=0.D0
c            sumt3=0.D0
c         end if
c         
c         Do NS=2,NSMAX
c            su1=su1+PN(NS)*PZ(NS)**2/PA(NS)
c         END DO
c         WBC=14.8D0*PTPP(1)*(1-R1)**2
c     &        *( PA(NSBM)**(3.D0/2.D0)*su1/PN(1) )**(2.D0/3.D0)
c         write(6,*)"Wbc[keV],Wb_init",WBC,PTT(1)
c         write(6,*)"Wb/Wbc, PC3=>1, PC3=>2"
c         sumt1=sumt1+PPCT2(NTG2,1)
c         sumt2=sumt2+PPCT2(NTG2,2)
c         sumt3=sumt3+PPCT(NTG2)
c         write(6,999)PTT(NTG2)/WBC
c     &        ,PPCT2(NTG2,1)/PPCT(NTG2),PPCT2(NTG2,2)/PPCT(NTG2)
c         write(6,*)" "
c         write(6,999)PTT(1)/WBC,sumt1/sumt3,sumt2/sumt3
      END IF
c-----end of beam check----
c      write(6,*)DCTT2(1,10,1,1),DCTT2(1,10,1,2)
c      write(6,*)FCPP2(1,10,1,1),FCPP2(1,10,1,2)
c      write(6,*)NTG2,PTPR(NSFP)*(1.D0-R1**2)
c      write(6,*)NTG2,DCPP2(1,10,1,1),DCPP2(1,10,1,2)
c-----slowing down-------------
      IF(NTG2.GE.NSMAX+1.and.NSFP.eq.NSFPMA+1)THEN
         rtaue=1.09D16*SQRT(PTT2(NTG3,1)**3)
     &        /PNT(NTG2)/1.d20/17.D0
         rtauei=rtaue*PA(2)/PA(1)*0.5D0
         rtaui =6.60D17*SQRT(PA(2)/PA(2)*PTT2(NTG3,2)**3)
     &        /PNT(NTG2)/17.D0/1.D20
      RTAUEI2=
     &     -(PTT2(NTG3-1,1)-PTT2(NTG3-1,2))/
     &     (PTT2(NTG3,1)-PTT2(NTG3-1,1))
     &*DELT

      RTAUIE2=
     &     -(PTT2(NTG3-1,2)-PTT2(NTG3-1,1))/
     &     (PTT2(NTG3,2)-PTT2(NTG3-1,2))
     &*DELT
      write(6,1582)rtaue,rtaui,rtauei
      write(6,*)"slowing e-i from diff",rtauei2
      write(6,*)"slowing i-e from diff",rtauie2
      write(6,*)"ratio", rtauei2/rtauie2
c      write(7,*)NTG3,rtauei2,rtauie2
      END IF
 1582 FORMAT("tau_e, tau_i, tau_ei [s]", 3E14.5)
c-----------------------------

      Do NS=1,NSMAX
         write(6,99) NSFP,NS,PPCT2(NTG2,NS)
      END DO


 99   FORMAT(1H ," PC[MW] ",I2," to ",I2," = ",E14.7)
 999  FORMAT(f14.6,2E14.6)
 1467 FORMAT(1H ,"T_",I1,"[keV]= ",E13.6)

      end if

C      WRITE(6,'(1PE12.5)') (RNS(NR),NR=1,NRMAX)
C
      IF(TIMEFP.GT.0.D0) THEN
         WRITE(6,103) PPCT(NTG2),PPWT(NTG2),PPET(NTG2)
         IF(ABS(PPWT(NTG2)).GT.0.D0) THEN
            WRITE(6,104) PLHT(NTG2),PFWT(NTG2),PECT(NTG2),
     &                   PIT(NTG2)/PPWT(NTG2)
            DO NRAY=1,NRAYMX
               EFFIP =PIT(NTG2)/PPWT(NTG2)
               AVENE =PN(1)
               RGAMMEF=AVENE*RR*EFFIP
            WRITE(6,201) RAYIN(1,NRAY),RAYIN(6,NRAY),
     &                   RAYIN(7,NRAY),RGAMMEF
            ENDDO
C
            DO NS=1,NSMAX
            DO 1000 NR=1,NRMAX
               RPN=RPWT(NR,NTG1)*AMFP*1.D6
     &               /(RNFP0*RNFP(NR)*1.D20
     &                *PTFP(NR)*PTFP(NR)*RNUD(NR,NS))
               RJN=RJT(NR,NTG1)*AMFP*1.D6
     &               /(RNFP0*RNFP(NR)*1.D20*AEFP*PTFP(NR))
               IF(ABS(RPN).LT.1.D-70) THEN
                  RJP=0.D0
               ELSE
                  RJP=RJN/RPN
               ENDIF
               WRITE(6,106) RM(NR),RPN,RJN,RJP
 1000       CONTINUE
            ENDDO
            DO 1100 NR=1,NRMAX
               WRITE(6,105) RM(NR),RPCT(NR,NTG1),RPWT(NR,NTG1),
     &                             RJT(NR,NTG1)
 1100       CONTINUE
         ENDIF
C
         IF (ABS(PPET(NTG2)).GT.0.D0) THEN
            DO NS=1,NSMAX
            DO 2000 NR=1,NRMAX
               REN=RET(NR,NTG1)*AEFP/(RNUD(NR,NS)*PTFP(NR))
               RJN=RJT(NR,NTG1)*AMFP*1.D6
     &               /(RNFP0*RNFP(NR)*1.D20*AEFP*PTFP(NR))
               IF(ABS(REN).LT.1.D-70) THEN
                  RJE=0.D0
               ELSE
                  RJE=RJN/REN
               ENDIF
               WRITE(6,107) RM(NR),REN,RJN,RJE
 2000       CONTINUE
            ENDDO
         ENDIF
      ENDIF

      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      WRITE(6,*)" "
      write(6,1333) THETA0, RTFP(1),PTT2(1,NSFP)/RTFP(1)
      write(6,1334) THETA0, RNFP(1),PNT2(NTG3,NSFP)/RNFP(1)
      write(6,*) RNFP(1), RN(1),RNT(1,1)
      WRITE(6,*)" "
 1333 FORMAT(1H ,"THETA0=",E15.7,"   T_given=",E15.7,"   ratio=",E15.7)
 1334 FORMAT(1H ,"THETA0=",E15.7,"   n_given=",E15.7,"   ratio=",E15.7)

      IF(NSFP.eq.NSFPMA)THEN
         rtotI=0.D0
         rtotW=0.D0
         DO NS=1,NSMAX
            write(6,1467)  NS, PTT2(NTG3,NS)
            rtotI = rtotI + PIT2(NTG3,NS)
            rtotW = rtotW + PWT2(NTG3,NS)
         END DO
         PITT(NTG3) = rtotI
         PWTT(NTG3) = rtotW
         WRITE(6,1600)  PITT(NTG3)
         WRITE(6,1590)  PWTT(NTG3)
         WRITE(6,*) "NTG2=",NTG2,"     NTG3=",NTG3,"  NTG1=" ,NTG1
      END IF

 1600 FORMAT(1H ,'TOTAL CURRENT [MA]      = ',E15.7)
 1590 FORMAT(1H ,'TOTAL STORED ENERGY [MJ]= ',E15.7)


      IF(mod(NTG2,NSFPMA-NSFPMI+1).eq.1)THEN
         write(6,*) "-------------------------------------------------"
      ELSE
         WRITE(6,*) "NTG2=",NTG2,"     NTG3=",NTG3,"  NTG1=" ,NTG1
         WRITE(6,*) "    "
      END IF


      RETURN
C
  101 FORMAT(1H ,' TIME=',F12.3,' MS')
  102 FORMAT(1H ,' N [20]=',1PE11.4,' I [MA]=',1PE11.4,
     &           ' W [MJ]=',1PE11.4,' T[keV]=',1PE12.5)
  103 FORMAT(1H ,' PC[MW]=',1PE11.4,' PW[MW]=',1PE11.4,
     &           ' PE[MW]=',1PE11.4)
  104 FORMAT(1H ,' PLH   =',1PE11.4,' PFW   =',1PE11.4,
     &           ' PEC   =',1PE11.4,' I/P   =',1PE11.4)
  201 FORMAT(1H ,' F[MHz]=',1PE11.4,' THP[D]=',1PE11.4,
     &           ' THT[D]=',1PE11.4,' NOREFF=',1PE11.4)  
  105 FORMAT(1H ,' RHO   =',F6.3,5X,' PC    =',1PE11.4,
     &           ' PW    =',1PE11.4,' J     =',1PE11.4)
  106 FORMAT(1H ,' RHO   =',F6.3,5X,' PN    =',1PE11.4,
     &           ' JN    =',1PE11.4,' J/P   =',1PE11.4)
  107 FORMAT(1H ,' RHO   =',F6.3,5X,' EN    =',1PE11.4,
     &           ' JN    =',1PE11.4,' J/E   =',1PE11.4)
      END

