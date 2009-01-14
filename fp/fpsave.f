C     $Id$
C
C *************************
C     SAVE DATA ROUTINE
C *************************
C
      SUBROUTINE FPSSUB
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION RSUM10(NSBMAX)

      DO NR=1,NRMAX
         DO NSA=1,NSAMAX
            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            RSUM5=0.D0
            RSUM6=0.D0
            RSUM7=0.D0
            RSUM8=0.D0
            RSUM9=0.D0
            RSUM123=0.D0
            DO NSB=1,NSBMAX
               RSUM10(NSB)=0.D0
            ENDDO

            NS=NS_NSA(NSA)
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
               END DO
            ENDDO
C
            IF(MODELR.EQ.0) THEN
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM2 = RSUM2
     &                     +VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                                 *PM(NP)*COSM(NTH)
                     RSUM3 = RSUM3
     &                     +VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                                 *0.5D0*PM(NP)**2
                  END DO
               ENDDO
            ELSE
               DO NP=1,NPMAX
                  PV=SQRT(1.D0+THETA0(NSA)*PM(NP)**2)
                  DO NTH=1,NTHMAX
                     RSUM2 = RSUM2
     &                     +VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                                 *PM(NP)*COSM(NTH)/PV
                     RSUM3 = RSUM3
     &                    +VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                                 *(PV-1.D0)/THETA0(NSA)
                     RSUM123 = RSUM123
     &                    +VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                                 *(PV-1.D0)/THETA0(NSA)
                  END DO
               END DO
            ENDIF

c      AMFDL=PA(NS)*AMP
c      AEFDL=PZ(NS)*AEE
c      RNFD0L=PN(NS)
c      RTFDL=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
c      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
c         THETAL=THETA0(NSA)*RTFDL/RTFD0L
c         Z=1.D0/THETAL
c               DKBSL1=BESKN(1,Z)
c               DKBSL2=BESKN(2,Z)
c               DKBSL1=1.D0+3.D0/8.D0/Z -15.D0/128.D0/z**2
c               DKBSL2=1.D0+15.D0/8.D0/z +105.D0/128.D0/z**2

c              write(*,*)Rsum123, (3.D0/z+DKBSL1/DKBSL2)/THETA0(NSA)
c     &              *RN(NS)/PN(NS)
c     &              ,Rsum123/(3.D0/z+DKBSL1/DKBSL2)*THETA0(NSA)
c     &              *PN(NS)/RN(NS)
c     &              ,1.D0/z
c     &              ,THETA0(NSA),RN(NS)
c     &             ,Rsum123/(3.D0+DKBSL1/DKBSL2/THETA0(NSA))/THETA0(NSA)

            IF(NTG2.EQ.1)      CALL FPWEIGHT(NSA,IERR)
            DO NP=2,NPMAX
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP)**2)
               DO NTH=1,NTHMAX
                  WPL=WEIGHP(NTH  ,NP,NR,NSA)
                  IF(NTH.EQ.1) THEN
                     WPM=0.D0
                  ELSE
                     WPM=WEIGHP(NTH-1,NP,NR,NSA)
                  ENDIF
                  IF(NTH.EQ.NTHMAX) THEN
                     WPP=0.D0
                  ELSE
                     WPP=WEIGHP(NTH+1,NP,NR,NSA)
                  ENDIF
                  DFP=    PG(NP)*RLAMDA(NTH,NR)
     &                 /DELP
     &                   *(FNS(NTH,NP,NR,NS)-FNS(NTH,NP-1,NR,NS))
                  IF(NTH.EQ.1) THEN
                     DFT=    1.D0/(2.D0*DELTH)
     &                   *( RLAMDA(NTH+1,NR)*
     &                      ((1.D0-WPP)*FNS(NTH+1,NP  ,NR,NS)
     &                             +WPP *FNS(NTH+1,NP-1,NR,NS)))
                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/(2.D0*DELTH)
     &                   *(- RLAMDA(NTH-1,NR)
     &                      *((1.D0-WPM)*FNS(NTH-1,NP  ,NR,NS)
     &                             +WPM *FNS(NTH-1,NP-1,NR,NS)))
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)
     &                   *( RLAMDA(NTH+1,NR)*
     &                      ((1.D0-WPP)*FNS(NTH+1,NP  ,NR,NS)
     &                             +WPP *FNS(NTH+1,NP-1,NR,NS))
     &                     - RLAMDA(NTH-1,NR)
     &                      *((1.D0-WPM)*FNS(NTH-1,NP  ,NR,NS)
     &                             +WPM *FNS(NTH-1,NP-1,NR,NS)))
                  ENDIF
                  FFP=    PG(NP)*RLAMDA(NTH,NR)
     &                   *((1.D0-WPL)*FNS(NTH  ,NP  ,NR,NS)
     &                          +WPL *FNS(NTH  ,NP-1,NR,NS))  
                  RSUM4 = RSUM4+PG(NP)**2*SINM(NTH)/PV
     &                   *(DCPP(NTH,NP,NR,NSA)*DFP
     &                    +DCPT(NTH,NP,NR,NSA)*DFT
     &                    -FCPP(NTH,NP,NR,NSA)*FFP
     &                 )
                  RSUM5 = RSUM5+PG(NP)**2*SINM(NTH)/PV
     &                   *(DWPP(NTH,NP,NR,NSA)*DFP
     &                    +DWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM6 = RSUM6-PG(NP)**2*SINM(NTH)/PV
     &                   *(FEPP(NTH,NP,NR,NSA)*FFP)
                  RSUM7 = RSUM7+PG(NP)**2*SINM(NTH)/PV
     &                   *(DWLHPP(NTH,NP,NR,NSA)*DFP
     &                    +DWLHPT(NTH,NP,NR,NSA)*DFT)
                  RSUM8 = RSUM8+PG(NP)**2*SINM(NTH)/PV
     &                   *(DWFWPP(NTH,NP,NR,NSA)*DFP
     &                    +DWFWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM9 = RSUM9+PG(NP)**2*SINM(NTH)/PV
     &                   *(DWECPP(NTH,NP,NR,NSA)*DFP
     &                    +DWECPT(NTH,NP,NR,NSA)*DFT)
                  DO NSB=1,NSBMAX
                     RSUM10(NSB)=RSUM10(NSB)+PG(NP)**2*SINM(NTH)/PV
     &                   *(DCPP2(NTH,NP,NR,NSB,NSA)*DFP
     &                    +DCPT2(NTH,NP,NR,NSB,NSA)*DFT
     &                    -FCPP2(NTH,NP,NR,NSB,NSA)*FFP)
                  END DO
c                  WRITE(*,*)NP,NTH,FFP,DCPT(NTH,NP,NR,NSA)*DFT
               ENDDO
            ENDDO

            FACT=RNFP0(NSA)*1.D20
            RNS(NR,NSA) = RSUM1*FACT                   *1.D-20
            RJS(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA)
     &                    /AMFP(NSA)*1.D-6

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RWS(NR,NSA) = RSUM3*FACT               *1.D-6
            RWS123(NR,NSA) = RSUM123*FACT                   *1.D-6
            RPCS(NR,NSA)=-RSUM4*FACT*2.D0*PI*DELP*DELTH *1.D-6
            RPWS(NR,NSA)=-RSUM5*FACT*2.D0*PI*DELP*DELTH *1.D-6 
            RPES(NR,NSA)=-RSUM6*FACT*2.D0*PI*DELP*DELTH *1.D-6
            RLHS(NR,NSA)=-RSUM7*FACT*2.D0*PI*DELP*DELTH *1.D-6
            RFWS(NR,NSA)=-RSUM8*FACT*2.D0*PI*DELP*DELTH *1.D-6
            RECS(NR,NSA)=-RSUM9*FACT*2.D0*PI*DELP*DELTH

            DO NSB=1,NSBMAX
               RPCS2(NR,NSB,NSA)=-RSUM10(NSB)
     &                          *FACT*2.D0*PI*DELP*DELTH *1.D-6
            END DO
         ENDDO
      ENDDO

      ISAVE=1
      RETURN
      END
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
      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RNT(NR,NSA,NTG1) = RNS(NR,NSA)
            RJT(NR,NSA,NTG1) = RJS(NR,NSA)
            RWT(NR,NSA,NTG1) = RWS(NR,NSA)
            RPCT(NR,NSA,NTG1)= RPCS(NR,NSA)
            RPWT(NR,NSA,NTG1)= RPWS(NR,NSA)
            RPET(NR,NSA,NTG1)= RPES(NR,NSA)
            RLHT(NR,NSA,NTG1)= RLHS(NR,NSA)
            RFWT(NR,NSA,NTG1)= RFWS(NR,NSA)
            RECT(NR,NSA,NTG1)= RECS(NR,NSA)
            DO NSB=1,NSBMAX
               RPCT2(NR,NSB,NSA,NTG1)= RPCS2(NR,NSB,NSA)
            END DO
C
            IF(RNS(NR,NSA).NE.0.D0) THEN
               RTT(NR,NSA,NTG1) = RWS(NR,NSA)*1.D6
     &                            /(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
            ELSE
               RTT(NR,NSA,NTG1) = 0.D0
            ENDIF
            RET(NR,NTG1) = E1(NR)
            RS=RSRHON(RM(NR))
            RQT(NR,NTG1) = RS*BB*2.D0/(RR*(BP(NR)+BP(NR+1)))
         ENDDO
      ENDDO

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
C
      IF(ISAVE.EQ.0) CALL FPSSUB
C
      PTG(NTG2)=TIMEFP
C
      DO NSA=1,NSAMAX
         PNT(NSA,NTG2)=0.D0
         PIT(NSA,NTG2)=0.D0
         PWT(NSA,NTG2)=0.D0
         PPCT(NSA,NTG2)=0.D0
         PPWT(NSA,NTG2)=0.D0
         PPET(NSA,NTG2)=0.D0
         PLHT(NSA,NTG2)=0.D0
         PFWT(NSA,NTG2)=0.D0
         PECT(NSA,NTG2)=0.D0
         PWT2(NSA,NTG2)=0.D0
         DO NSB=1,NSBMAX
            PPCT2(NSB,NSA,NTG2)= 0.D0
         END DO
      ENDDO
C
      TVOL=0.D0
      DO NR=1,NRMAX
         RHOL=RM(NR)
         RHOL1=RG(NR)
         RHOL2=RG(NR+1)
         TVOL=TVOL+2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1))
      ENDDO

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RHOL=RM(NR)
            RHOL1=RG(NR)
            RHOL2=RG(NR+1)
            FACT=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1))
            PNT(NSA,NTG2) =PNT(NSA,NTG2) +RNS(NR,NSA)*FACT
            PIT(NSA,NTG2) =PIT(NSA,NTG2) +RJS(NR,NSA)*FACT
            PWT(NSA,NTG2) =PWT(NSA,NTG2) +RWS(NR,NSA)*FACT
            PPCT(NSA,NTG2)=PPCT(NSA,NTG2)+RPCS(NR,NSA)*FACT
            PPWT(NSA,NTG2)=PPWT(NSA,NTG2)+RPWS(NR,NSA)*FACT
            PPET(NSA,NTG2)=PPET(NSA,NTG2)+RPES(NR,NSA)*FACT
            PLHT(NSA,NTG2)=PLHT(NSA,NTG2)+RLHS(NR,NSA)*FACT
            PFWT(NSA,NTG2)=PFWT(NSA,NTG2)+RFWS(NR,NSA)*FACT
            PECT(NSA,NTG2)=PECT(NSA,NTG2)+RECS(NR,NSA)*FACT
            IF(MODELR.eq.1) then
               CALL FPNEWTON(NR,NSA,rtemp)
            else
               EAVE=RWS123(NR,NSA)*AMFP(NSA)*THETA0(NSA)
     &              /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
               THETAL=2.d0*EAVE/3.d0
               rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
            endif
            PWT2(NSA,NTG2) =PWT2(NSA,NTG2) +rtemp*FACT/1.D6
     &                            *(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
            DO NSB=1,NSBMAX
               PPCT2(NSB,NSA,NTG2)=PPCT2(NSB,NSA,NTG2)
     &                            +RPCS2(NR,NSB,NSA)*FACT
            END DO
         ENDDO
         PIT(NSA,NTG2) =PIT(NSA,NTG2) 
         PPCT(NSA,NTG2)=PPCT(NSA,NTG2)*2.D0*PI*RR
         PPWT(NSA,NTG2)=PPWT(NSA,NTG2)*2.D0*PI*RR
         PPET(NSA,NTG2)=PPET(NSA,NTG2)*2.D0*PI*RR
         PLHT(NSA,NTG2)=PLHT(NSA,NTG2)*2.D0*PI*RR
         PFWT(NSA,NTG2)=PFWT(NSA,NTG2)*2.D0*PI*RR
         PECT(NSA,NTG2)=PECT(NSA,NTG2)*2.D0*PI*RR
         IF(PNT(NSA,NTG2).NE.0.d0) THEN
            PTT(NSA,NTG2) =PWT(NSA,NTG2)*1.D6
     &           /(1.5D0*PNT(NSA,NTG2)*1.D20*AEE*1.D3)
            PTT2(NSA,NTG2) =PWT2(NSA,NTG2)*1.D6
     &           /(1.5D0*PNT(NSA,NTG2)*1.D20*AEE*1.D3)
         ELSE
            PTT(NSA,NTG2)=0.D0
            PTT2(NSA,NTG2)=0.D0
         ENDIF
         PNT(NSA,NTG2) =PNT(NSA,NTG2)/TVOL
         PWT(NSA,NTG2) =PWT(NSA,NTG2) *2.D0*PI*RR
         PWT2(NSA,NTG2) =PWT2(NSA,NTG2) *2.D0*PI*RR
         DO NSB=1,NSBMAX
            PPCT2(NSB,NSA,NTG2)=PPCT2(NSB,NSA,NTG2)*2.D0*PI*RR
         END DO
      ENDDO
         
      RETURN
      END
C
C ***********************************************************
C
C                         RESULT
C
C ***********************************************************
C
      SUBROUTINE FPWRT2
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
C
      WRITE(6,*)"--------------------------------------------"
      WRITE(6,*)"-----Global data"
      WRITE(6,101) PTG(NTG2)*1000

      DO NSA=1,NSAMAX
         IF(MODELR.eq.0)THEN
            WRITE(6,102) NSA,NS_NSA(NSA),
     &        PNT(NSA,NTG2),PTT(NSA,NTG2),PWT(NSA,NTG2),PIT(NSA,NTG2)
         ELSE
            WRITE(6,102) NSA,NS_NSA(NSA),
     &        PNT(NSA,NTG2),PTT2(NSA,NTG2),PWT(NSA,NTG2),PIT(NSA,NTG2)
         END IF
c         WRITE(6,103) PPCT(NSA,NTG2),PPWT(NSA,NTG2),PPET(NSA,NTG2)
         IF(NSBMAX.GT.1) THEN
c            write(6,104)(PPCT2(NSB,NSA,NTG2),NSB=1,NSBMAX)
         ENDIF
      ENDDO

      rtotalPW=0.D0
      rtotalPC=0.D0
      DO NSA=1,NSAMAX
         WRITE(6,103) NSA,NS_NSA(NSA),
     &        PPCT(NSA,NTG2),PPWT(NSA,NTG2),PPET(NSA,NTG2)
         rtotalPW=rtotalPW + PPWT(NSA,NTG2)
         rtotalPC=rtotalPC + PPCT(NSA,NTG2)
      END DO
      DO NSA=1,NSAMAX
         IF(NSBMAX.GT.1) THEN
            write(6,104) NSA,NS_NSA(NSA),
     &           (PPCT2(NSB,NSA,NTG2),NSB=1,NSBMAX)
         ENDIF
      END DO
      write(*,105) rtotalpw
      write(*,107) rtotalPC

c      write(8,106) PTG(NTG2)*1000, PPCT(1,NTG2),PPCT(2,NTG2)
c     &  ,PPCT2(1,1,NTG2),PPCT2(2,1,NTG2),PPCT2(1,2,NTG2),PPCT2(2,2,NTG2)
c     &     ,PPET(1,NTG2),PPET(2,NTG2)

c      NP2=99
c      WRITE(*,*)ITL(1)-2,DWPP(35,NP2,1,3),DCPP(35,NP2,1,3)
c      WRITE(*,*)ITL(1)-1,DWPP(36,NP2,1,3),DCPP(36,NP2,1,3)
c      WRITE(*,*)ITL(1),DWPP(37,NP2,1,3),DCPP(37,NP2,1,3)
c      WRITE(*,*)ITL(1)+1,DWPP(38,NP2,1,3),DCPP(38,NP2,1,3)
c      WRITE(*,*)ITL(1)+2,DWPP(39,NP2,1,3),DCPP(39,NP2,1,3)

c      WRITE(*,*)FNS(36,NP2,1,3),RLAMDA(36,1),RLAMDC(36,1)
c      WRITE(*,*)FNS(37,NP2,1,3),RLAMDA(37,1),RLAMDC(37,1)
c      WRITE(*,*)FNS(38,NP2,1,3),RLAMDA(38,1),RLAMDC(38,1)

c      DO NTH=1,NTHMAX
c         WRITE(*,1001)NTH,
c     &        DWPT(NTH,50,1,3),DWTP(NTH,50,1,3),DWTT(NTH,50,1,3)
c     &        ,DWPP(NTH,50,1,3),DWTT(NTH,50,1,3)
c     &        ,DWPT(NTH,50,1,3),FNS(NTH,50,1,3)
c      DO NP=1,NPMAX
c         WRITE(*,1001)NP,
c     &  FNS(ITL(1)-2,NP,1,3),FNS(ITL(1)-1,NP,1,3),FNS(ITL(1),NP,1,3),
c     &        FNS(ITL(1)+1,NP,1,3),FNS(ITL(1)+2,NP,1,3)
c     &   DTP(ITL(1)-2,NP,1,3),DTP(ITL(1)-1,NP,1,3),DTP(ITL(1),NP,1,3)
c     &        ,DTP(ITL(1)+1,NP,1,3),DTP(ITL(1)+2,NP,1,3)
c      END DO
 1001 FORMAT(I4,9E14.6)
      RETURN
  101 FORMAT(' TIME=',F12.3,' ms')
  102 FORMAT(' NSA,NS=',2I2,' n,T,W,I=',1P4E12.4)
c  103 FORMAT('             PC,PW,PE=',11X,1P3E12.4)
c  104 FORMAT('             PCAB    =',11X,1P3E12.4)
  103 FORMAT('        ',2I2,' PC,PW,PE=',11X,1P3E12.4)
  104 FORMAT('        ',2I2,' PCAB    =',11X,1P3E12.4)

 105  FORMAT('total absorption [MW]', E12.4)
 106  FORMAT(F12.4, 8E12.4)
 107  FORMAT('total collision power [MW]', E12.4)
      END

C ***********************************************************

      SUBROUTINE FPWRT1
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
C
      WRITE(6,*)"-----Radial profile data"
      WRITE(6,101) TIMEFP*1000

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
      RTFDL=RTFP(NR,NSA)
      RTFD0L=(PTPR(NSA)+2.D0*PTPP(NSA))/3.D0
         THETAL=THETA0(NSA)*RTFDL/RTFD0L

c            WRITE(6,102) NSA,NS_NSA(NSA),
c     &                   RM(NR),RNT(NR,NSA,NTG1),RTT(NR,NSA,NTG1),
c     &                          RJT(NR,NSA,NTG1),RPCT(NR,NSA,NTG1),
c     &                          RPET(NR,NSA,NTG1)
c            IF(RPWT(NR,NSA,NTG1).GT.0.D0) THEN
c            WRITE(6,103)        RPWT(NR,NSA,NTG1),RLHT(NR,NSA,NTG1),
c     &                          RFWT(NR,NSA,NTG1),RECT(NR,NSA,NTG1)
c            ENDIF
            IF(MODELR.eq.1)THEN
               CALL FPNEWTON(NR,NSA,rtemp)
               RTT(NR,NSA,NTG1)=rtemp
            END IF

            IF(RPWT(NR,NSA,NTG1).GT.0.D0) THEN
               WRITE(6,104) NSA,NS_NSA(NSA),
     &              RM(NR),RNT(NR,NSA,NTG1),RTT(NR,NSA,NTG1),
     &              RJT(NR,NSA,NTG1),RPCT(NR,NSA,NTG1),
     &              RPET(NR,NSA,NTG1),RPWT(NR,NSA,NTG1),
     &              RLHT(NR,NSA,NTG1),
     &              RFWT(NR,NSA,NTG1),RECT(NR,NSA,NTG1)
c               write(*,103)THETAL,rtemp, RTFP(NR,NSA)
            ELSE
               WRITE(6,102) NSA,NS_NSA(NSA),
     &              RM(NR),RNT(NR,NSA,NTG1),RTT(NR,NSA,NTG1),
     &              RJT(NR,NSA,NTG1),RPCT(NR,NSA,NTG1),
     &              RPET(NR,NSA,NTG1)
c               write(*,103)THETAL,rtemp, RTFP(NR,NSA)
            ENDIF
c               RTT(NR,NSA,NTG1)=rtemp
         ENDDO
      ENDDO
      RETURN
  101 FORMAT(' TIME=',F12.3,' ms'/
     &     'NSA/NS',5X,'RM',10X,' n',8X,' T//PW',6X,
     &     ' j//PLH',5X,'PC//PIC',5X,'PE//PEC')
  102 FORMAT(2I3,1P6E12.4)
  103 FORMAT(30X,1P4E14.6)
  104 FORMAT(2I3,1P10E12.4) 
      END
C ***********************************************************
C
      SUBROUTINE FPWRT3
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
C
      WRITE(6,101) TIMEFP*1000

c$$$      if(NTG2.ne.1)then
c$$$c-----check of conductivity--------
      IF(NTG2.ne.1)then
         FACT1 = 
     &2.D0*PI*RSRHON(RM(NRMAX))*(RSRHON(RG(2))-RSRHON(RG(NRMAX)))
c         FACT2 = 2.D0*PI*RR
c         rntv = 1.6D0/1.D19*1.5D0*
c     &FACT1*FACT2*(PTT(NTG2)-PTT(NTG2-1))*PNT(NTG2)*1.D3*1.D20/DELT

c      write(6,*) "Pw+Pc",(PPWT(NTG2)+PPCT(NTG2))*1.D6,"n*delta T*V",rntv
c      write(6,*) "c*epsilon*RabsE*S",RABSE**2*FACT1*8.8D0*VC/1.D12
c     & ,"Pw[W]",PPWT(NTG2)*1.D6
c      write(6,*) "Pc+PE+Pw",PPCT(NTG2)+PPET(NTG2)+PPWT(NTG2),
c     &     "(Pc+Pe+Pw)/Pc",(PPCT(NTG2)+PPET(NTG2)+PPWT(NTG2))/PPCT(NTG2)
c      write(6,*) "Pc+-Pc-/Pc",(PPCT(NTG2)-PPCT(NTG2-1))/PPCT(NTG2)
C     Spitzer

      Do NSA=1,NSAMAX
c      Do NSB=1,NSBMAX
         rnute=RNUD(1,1,NSA)*SQRT(2.D0)*RNFP(1,NSA)/RNFP0(NSA)
     &        *(PTFP0(NSA)/PTFD(1,NSA))**2
         resist2=RNUTE*AMFP(NSA)/RNFP(1,NSA)/AEFP(NSA)**2/1.D20
      resist = 1.D0/(RNFP0(NSA)/RNFP(1,NSA)*1.D20*AEFP(NSA)**2/AMFP(NSA)
     &    /RNUD(1,1,NSA)*RNFD(1,NSA)/RNFP0(NSA) ) *SQRT(2.D0)/rnfp0(NSA)
c      if(E0.ne.0.d0) 
c     &    write(6,*) "J/E*eta*1.D6", PIT(NSA,NTG2)/E0*1.D6/FACT1*resist2
c     &     ,resist,resist2
c     &     ,RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)
c     & ,"THETA0", (PTFP0(NSA)/(AMFP(NSA)*VC))**2
c     &     ,DCTT2(2,2,1,NSB,NSA),NSA,NSB
c      END DO
      END DO
      dtaue=1.09D16*(PTT(1,NTG2))**(1.5D0)/PNT(1,NTG2)/PZ(2)**2
     &     /15.D0*1.D-20
      if(nsamax.gt.1) then
      dtaui=6.6D17*(PTT(2,NTG2))**(1.5D0)/PNT(2,NTG2)/PZ(2)**4
     &     /15.D0*1.D-20*(AMFP(2)/AMFP(2))**0.5D0
      else
         dtaui=0.d0
      endif
c      write(6,*)"tau_e tau_i[ms]",dtaue*1.D3,dtaui*1.D3
      end if
c$$$c----end of conductivity check---------
c$$$c----check of beam power-------
c$$$      IF(NSBM.ne.0)THEN
c$$$c         su1=0.D0
c$$$c         if(NTG2.eq.1)then
c$$$c            sumt1=0.D0
c$$$c            sumt2=0.D0
c$$$c            sumt3=0.D0
c$$$c         end if
c$$$c         
c$$$c         Do NS=2,NSMAX
c$$$c            su1=su1+PN(NS)*PZ(NS)**2/PA(NS)
c$$$c         END DO
c$$$c         WBC=14.8D0*PTPP(1)*(1-R1)**2
c$$$c     &        *( PA(NSBM)**(3.D0/2.D0)*su1/PN(1) )**(2.D0/3.D0)
c$$$c         write(6,*)"Wbc[keV],Wb_init",WBC,PTT(1)
c$$$c         write(6,*)"Wb/Wbc, PC3=>1, PC3=>2"
c$$$c         sumt1=sumt1+PPCT2(NTG2,1)
c$$$c         sumt2=sumt2+PPCT2(NTG2,2)
c$$$c         sumt3=sumt3+PPCT(NTG2)
c$$$c         write(6,999)PTT(NTG2)/WBC
c$$$c     &        ,PPCT2(NTG2,1)/PPCT(NTG2),PPCT2(NTG2,2)/PPCT(NTG2)
c$$$c         write(6,*)" "
c$$$c         write(6,999)PTT(1)/WBC,sumt1/sumt3,sumt2/sumt3
c$$$      END IF
c$$$c-----end of beam check----
c$$$c      write(6,*)DCTT2(1,10,1,1),DCTT2(1,10,1,2)
c$$$c      write(6,*)FCPP2(1,10,1,1),FCPP2(1,10,1,2)
c$$$c      write(6,*)NTG2,PTPR(NSFP)*(1.D0-R1**2)
c$$$c      write(6,*)NTG2,DCPP2(1,10,1,1),DCPP2(1,10,1,2)
c$$$c-----slowing down-------------
c$$$      IF(NTG2.GE.NSMAX+1.and.NSFP.eq.NSFPMA+1)THEN
c$$$         rtaue=1.09D16*SQRT(PTT2(NTG3,1)**3)
c$$$     &        /PNT(NTG2)/1.d20/17.D0
c$$$         rtauei=rtaue*PA(2)/PA(1)*0.5D0
c$$$         rtaui =6.60D17*SQRT(PA(2)/PA(2)*PTT2(NTG3,2)**3)
c$$$     &        /PNT(NTG2)/17.D0/1.D20
c$$$      RTAUEI2=
c$$$     &     -(PTT2(NTG3-1,1)-PTT2(NTG3-1,2))/
c$$$     &     (PTT2(NTG3,1)-PTT2(NTG3-1,1))
c$$$     &*DELT
c$$$
c$$$      RTAUIE2=
c$$$     &     -(PTT2(NTG3-1,2)-PTT2(NTG3-1,1))/
c$$$     &     (PTT2(NTG3,2)-PTT2(NTG3-1,2))
c$$$     &*DELT
c$$$      write(6,1582)rtaue,rtaui,rtauei
c$$$      write(6,*)"slowing e-i from diff",rtauei2
c$$$      write(6,*)"slowing i-e from diff",rtauie2
c$$$      write(6,*)"ratio", rtauei2/rtauie2
c$$$c      write(7,*)NTG3,rtauei2,rtauie2
c$$$      END IF
c$$$ 1582 FORMAT("tau_e, tau_i, tau_ei [s]", 3E14.5)
c$$$ 999  FORMAT(f14.6,2E14.6)
c$$$ 1467 FORMAT(1H ,"T_",I1,"[keV]= ",E13.6)
c$$$c-----------------------------

c$$$      IF(TIMEFP.GT.0.D0) THEN
c$$$         IF(ABS(PPWT(NTG2)).GT.0.D0) THEN
c$$$            WRITE(6,104) PLHT(NTG2),PFWT(NTG2),PECT(NTG2),
c$$$     &                   PIT(NTG2)/PPWT(NTG2)
c$$$            DO NRAY=1,NRAYMX
c$$$               EFFIP =PIT(NTG2)/PPWT(NTG2)
c$$$               AVENE =PN(1)
c$$$               RGAMMEF=AVENE*RR*EFFIP
c$$$            WRITE(6,201) RAYIN(1,NRAY),RAYIN(6,NRAY),
c$$$     &                   RAYIN(7,NRAY),RGAMMEF
c$$$            ENDDO
c$$$C
c$$$            DO NS=1,NSMAX
c$$$            DO 1000 NR=1,NRMAX
c$$$               RPN=RPWT(NR,NTG1)*AMFP*1.D6
c$$$     &               /(RNFP0*RNFP(NR)*1.D20
c$$$     &                *PTFP(NR)*PTFP(NR)*RNUD(NR,NS))
c$$$               RJN=RJT(NR,NTG1)*AMFP*1.D6
c$$$     &               /(RNFP0*RNFP(NR)*1.D20*AEFP*PTFP(NR))
c$$$               IF(ABS(RPN).LT.1.D-70) THEN
c$$$                  RJP=0.D0
c$$$               ELSE
c$$$                  RJP=RJN/RPN
c$$$               ENDIF
c$$$               WRITE(6,106) RM(NR),RPN,RJN,RJP
c$$$ 1000       CONTINUE
c$$$            ENDDO
c$$$            DO 1100 NR=1,NRMAX
c$$$               WRITE(6,105) RM(NR),RPCT(NR,NTG1),RPWT(NR,NTG1),
c$$$     &                             RJT(NR,NTG1)
c$$$ 1100       CONTINUE
c$$$         ENDIF
c$$$C
c$$$         IF (ABS(PPET(NTG2)).GT.0.D0) THEN
c$$$            DO NS=1,NSMAX
c$$$            DO 2000 NR=1,NRMAX
c$$$               REN=RET(NR,NTG1)*AEFP/(RNUD(NR,NS)*PTFP(NR))
c$$$               RJN=RJT(NR,NTG1)*AMFP*1.D6
c$$$     &               /(RNFP0*RNFP(NR)*1.D20*AEFP*PTFP(NR))
c$$$               IF(ABS(REN).LT.1.D-70) THEN
c$$$                  RJE=0.D0
c$$$               ELSE
c$$$                  RJE=RJN/REN
c$$$               ENDIF
c$$$               WRITE(6,107) RM(NR),REN,RJN,RJE
c$$$ 2000       CONTINUE
c$$$            ENDDO
c$$$         ENDIF
c$$$      ENDIF
c$$$
c$$$      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
c$$$      WRITE(6,*)" "
c$$$      write(6,1333) THETA0, RTFP(1),PTT2(1,NSFP)/RTFP(1)
c$$$      write(6,1334) THETA0, RNFP(1),PNT2(NTG3,NSFP)/RNFP(1)
c$$$      write(6,*) RNFP(1), RN(1),RNT(1,1)
c$$$      WRITE(6,*)" "
c$$$ 1333 FORMAT(1H ,"THETA0=",E15.7,"   T_given=",E15.7,"   ratio=",E15.7)
c$$$ 1334 FORMAT(1H ,"THETA0=",E15.7,"   n_given=",E15.7,"   ratio=",E15.7)
c$$$
c$$$      IF(NSFP.eq.NSFPMA)THEN
c$$$         rtotI=0.D0
c$$$         rtotW=0.D0
c$$$         DO NS=1,NSMAX
c$$$            write(6,1467)  NS, PTT2(NTG3,NS)
c$$$            rtotI = rtotI + PIT2(NTG3,NS)
c$$$            rtotW = rtotW + PWT2(NTG3,NS)
c$$$         END DO
c$$$         PITT(NTG3) = rtotI
c$$$         PWTT(NTG3) = rtotW
c$$$         WRITE(6,1600)  PITT(NTG3)
c$$$         WRITE(6,1590)  PWTT(NTG3)
c$$$         WRITE(6,*) "NTG2=",NTG2,"     NTG3=",NTG3,"  NTG1=" ,NTG1
c$$$      END IF
c$$$
c$$$ 1600 FORMAT(1H ,'TOTAL CURRENT [MA]      = ',E15.7)
c$$$ 1590 FORMAT(1H ,'TOTAL STORED ENERGY [MJ]= ',E15.7)
c$$$
c$$$
c$$$      IF(mod(NTG2,NSFPMA-NSFPMI+1).eq.1)THEN
c$$$         write(6,*) "-------------------------------------------------"
c$$$      ELSE
c$$$         WRITE(6,*) "NTG2=",NTG2,"     NTG3=",NTG3,"  NTG1=" ,NTG1
c$$$         WRITE(6,*) "    "
c$$$      END IF
c$$$

      RETURN
C
  101 FORMAT(1H ,' TIME=',F12.3,' ms')
c$$$  102 FORMAT(1H ,' N [20]=',1PE11.4,' I [MA]=',1PE11.4,
c$$$     &           ' W [MJ]=',1PE11.4,' T[keV]=',1PE12.5)
c$$$  103 FORMAT(1H ,' PC[MW]=',1PE11.4,' PW[MW]=',1PE11.4,
c$$$     &           ' PE[MW]=',1PE11.4)
c$$$  104 FORMAT(1H ,' PLH   =',1PE11.4,' PFW   =',1PE11.4,
c$$$     &           ' PEC   =',1PE11.4,' I/P   =',1PE11.4)
c$$$  201 FORMAT(1H ,' F[MHz]=',1PE11.4,' THP[D]=',1PE11.4,
c$$$     &           ' THT[D]=',1PE11.4,' NOREFF=',1PE11.4)  
c$$$  105 FORMAT(1H ,' RHO   =',F6.3,5X,' PC    =',1PE11.4,
c$$$     &           ' PW    =',1PE11.4,' J     =',1PE11.4)
c$$$  106 FORMAT(1H ,' RHO   =',F6.3,5X,' PN    =',1PE11.4,
c$$$     &           ' JN    =',1PE11.4,' J/P   =',1PE11.4)
c$$$  107 FORMAT(1H ,' RHO   =',F6.3,5X,' EN    =',1PE11.4,
c$$$     &           ' JN    =',1PE11.4,' J/E   =',1PE11.4)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Subroutine FPNEWTON(NR,NSA,rtemp)

      INCLUDE 'fpcomm.inc'

C-----Average kinetic energy
      EAVE=RWS123(NR,NSA)*AMFP(NSA)*THETA0(NSA)
     &     /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)


C-----initial value of THETAL
      THETAL=2.d0*EAVE/3.d0
      xtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)

      CALL XNEWTON(EAVE,THETAL,ncount)

      rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
      xeave=AMFP(NSA)*VC**2*EAVE/(AEE*1.D3)
c      write(6,'(3I5,1P3E12.4)') NSA,NR,ncount,xeave,xtemp,rtemp

      RETURN

      CONTAINS

      SUBROUTINE xnewton(eave,thetal,ncount)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8),intent(in):: eave
      REAL(8),intent(inout):: thetal
      INTEGer,intent(out):: ncount
      REAL(8),parameter:: eps=1.d-10
      REAL(8):: delthetal,thetalnew,epsthetal

C--------iteration start
      ncount=0
      DO while(ncount.le.100)
         ncount=ncount+1
         delthetal=-(rfunc(thetal)-eave)/dfunc(thetal)
         thetalnew=thetal+delthetal
         epsthetal=ABS(delthetal/thetal)

         thetal=thetalnew
         IF(epsthetal.le.eps) EXIT
      END DO
      RETURN
      END SUBROUTINE xnewton

      FUNCTION rfunc(thetal)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8):: thetal,rfunc
      REAL(8):: z,dkbsl1,dkbsl2,BESEKN
      z=1.D0/thetal
      dkbsl1=BESEKN(1,Z)
      dkbsl2=BESEKN(2,Z)
      rfunc= (dkbsl1 /dkbsl2 -1.D0+3.D0/Z)
      RETURN
      END FUNCTION rfunc

      FUNCTION rfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfuncp
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      rfuncp= dkbsl1 /dkbsl2 -1.D0+3.D0/Z
      RETURN
      END FUNCTION rfuncp
      
      FUNCTION dfunc(thetal)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8):: thetal,dfunc
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3,BESEKN
      z=1.D0/thetal
      dkbsl0=BESEKN(0,z)
      dkbsl1=BESEKN(1,z)
      dkbsl2=BESEKN(2,z)
      dkbsl3=BESEKN(3,z)
      dfunc =( (dkbsl0 +dkbsl2 )/dkbsl2
     &          -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2
     &      +3.d0
      RETURN
      END FUNCTION dfunc

      FUNCTION dfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfuncp
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=1.D0 -  1.D0/8.D0/z +   9.D0/128.D0/z**2
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      dkbsl3=1.D0 + 35.D0/8.D0/z + 945.D0/128.D0/z**2
      dfuncp =( (dkbsl0 +dkbsl2 )/dkbsl2
     &         -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2
     &      +3.d0
      RETURN
      END FUNCTION dfuncp

      end






