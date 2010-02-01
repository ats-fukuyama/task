C
C ******************************************************
C                       DPDKDT
C      dielectric tensor with Drift Kinetic effects
C                 
C                      2010/02/01             
C                           programed by T.OKAMOTO
C ******************************************************
C
      SUBROUTINE DPDKDT(CW,CKPR,NS,NTH,NP,NR,NCH1,NCH2,MM,CLDISP)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION CLDISP(9),CLDISP1(9),CLDISP2(9)
C
      CALL DPDKDTR(CW,CKPR,NS,NTH,NP,NR,NCH1,NCH2,MM,CLDISP1)
      CALL DPDKDTI(CW,CKPR,NS,NTH,NP,NR,NCH1,NCH2,MM,CLDISP2)
      DO I=1,9
         CLDISP(I)=CLDISP1(I)+CLDISP2(I)
      ENDDO
      RETURN
      END
C
C ******************************************************
C                       DPDKDTR
C ******************************************************
C
      SUBROUTINE DPDKDTR(CW,CKPR,NS,NTH,NP,NR,NCH1,NCH2,MM,CLDISP1)
C
C------------------------------------------------------!
C     NS : PARTICLE SPECIES
C     NTH: PITCH ANGLE MESH
C     NP : MOMENTUM NUMBER
C     NR : NODE NUMBER (RADIAL POSITION)
C     NCH1,NCH2: POLOIDAL MESH
C     MM : POLOIDAL NUMBER
C------------------------------------------------------!

      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP(9)
      DIMENSION DRFP1(NTHMAX,NPMAX,NRMAX)
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION DRRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DRCWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION DPFP1(NTHMAX,NPMAX,NRMAX)
      DIMENSION DRPFP(NTHMAX,NPMAX)
      DIMENSION DGP(NTHMAX,NPMAX)
      DIMENSION DGT(NTHMAX,NPMAX)
      DIMENSION TCSM2(NTHMAX)
      DIMENSION PV(NPMAX)
C
      DELPL  = 0.5D0
C
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
C      
      RSL=DELR*(NR-0.5D0)

!--------derivation in radial direction & velocity---------------!
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP1(NTH,NP,NR) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         DRFP(NTH,NP) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         DRRFP(NTH,NP)=(DRFP1(NTH,NP,NR+1)-DRFP1(NTH,NP,NR-1))
     &        /(2*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH, NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP1(NTH,NP,NR) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=1,NPMAX-2
               DPFP1(NTH,NP,NR) = (FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/(2*DELP)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP1(NTH,NP,NR) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
         DRPFP(NTH,NP) = ((DPFP1(NTH,NP,NR+1)
     &        - DPFP1(NTH,NP,NR-1))/(2*DELP))/FP(NTH,NP,NR)
      ENDDO
C           
C  
!--------------diffrerential of g(v,theta),DGP,DGT----------!
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX
         PV(NP)=PG(NP)/AM
         DGP(NTH,NP)=DKPP*PV(NP)*TSNM(NTH)-(MM*AA*PV(NP)**2*TCSM(NTH)
     &              *TSNM(NTH)*CCHI2(NCH2))/RSL    
         DGT(NTH,NP)=-DKPP*PV(NP)*TCSM(NTH)+(MM*AA*PV(NP)
     &        *(1+TCSM2(NTH))*CCHI2(NCH2))/RSL
      ENDDO
      ENDDO
C------------------------------------------------------!

C*****************PRINCIPAL VALUE***********************
C
C*************SUM1********************
C
      CINTG111 = (0.D0,0.D0)
      CINTG112 = (0.D0,0.D0)
      CINTG113 = (0.D0,0.D0)
      CINTG121 = (0.D0,0.D0)
      CINTG122 = (0.D0,0.D0)
      CINTG123 = (0.D0,0.D0)
      CINTG131 = (0.D0,0.D0)
      CINTG132 = (0.D0,0.D0)
      CINTG133 = (0.D0,0.D0)
C
      DO NP=1,NPMAX-1
      DO NTH=1,NTHMAX
C        
         PV(NP)=PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
            CDENX = CW-DKPRX+MM*AA*PV(NP)**2*(1+TCSM2(NTH))*CCHI2(NCH2)
     &                                                     /(2*RSL)
            CDEN  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*DELP)**2
     &                             +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
            CPART11=DPFP(NTH,NP)-CWAST(NTH,NP)/CW
C
            CPART12=CI*(MM*AA**2*PV(NP)**4*(1+TCSM2(NTH))**2*CCHI2(NCH2)
     &           *SCHI2(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/(4*CDENX))
C
            PCART13=CI*((AA/2)*PV(NP)**2*(1+TCSM2(NTH))*SCHI2(NCH2))
     &           *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &           +CWAST(NTH,NP)/(CW*RSL))
C      
            PART14=2*PI*CI*PV(NP)**3*TSNM(NTH)*COEF
     &            *(CPART11+CPART12+CPART13)*DELTH*DELP
C

            PI1=(AA/2)**2*PV(NP)**4*(1+TCSM2(NTH))**2
            PI2=(AA/2)*PV(NP)**4*(1+TCSM2(NTH))
            PI3=PV(NP)**4
C
C
         CINTG111 = CINTG111 + CDEN*CPART14*PI1*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG112 = CINTG112 + CDEN*CPART14*PI1*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG113 = CINTG113 + CDEN*CPART14*PI2*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG121 = CINTG121 + CDEN*CPART14*PI1*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG122 = CINTG122 + CDEN*CPART14*PI1*SCHI1(NCH1)*SCHI2(NCH2)
         CINTG123 = CINTG123 + CDEN*CPART14*PI2*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG131 = CINTG131 + CDEN*CPART14*PI2*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG132 = CINTG132 + CDEN*CPART14*PI2*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG133 = CINTG133 + CDEN*CPART14*PI3*CCHI1(NCH1)*CCHI2(NCH2)
      ENDDO
      ENDDO
C
C*************SUM2********************
C
      CINTG211 = (0.D0,0.D0)
      CINTG212 = (0.D0,0.D0)
      CINTG213 = (0.D0,0.D0)
      CINTG221 = (0.D0,0.D0)
      CINTG222 = (0.D0,0.D0)
      CINTG223 = (0.D0,0.D0)
      CINTG231 = (0.D0,0.D0)
      CINTG232 = (0.D0,0.D0)
      CINTG233 = (0.D0,0.D0)
C
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX-1
C
         PV(NP)=PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
         CDENX = CW-DKPRX+MM*AA*PV(NP)**2*(1+TCSM2(NTH)*CCHI2(NCH2))
     &        /(2*RSL)
         CDEN  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*DELP)**2
     &        +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART2= 2*PI*CI*PV(NP)**3*COEF
     &        *TSNM(NTH)*(AA/2)*PV(NP)**2*(1+TCSM2(NTH))*SCHI2(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
     &        *DELTH*DELP
C 
            PI1=(AA/2)**2*PV(NP)**4*(1+TCSM2(NTH))**2
            PI2=(AA/2)*PV(NP)**4*(1+TCSM2(NTH))
            PI3=PV(NP)**4
C
         CINTG211 = CINTG211 + CDEN*CPART2*PI1*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG212 = CINTG212 + CDEN*CPART2*PI1*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG213 = CINTG213 + CDEN*CPART2*PI2*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG221 = CINTG221 + CDEN*CPART2*PI1*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG222 = CINTG222 + CDEN*CPART2*PI1*SCHI1(NCH1)*SCHI2(NCH2)
         CINTG223 = CINTG223 + CDEN*CPART2*PI2*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG231 = CINTG231 + CDEN*CPART2*PI2*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG232 = CINTG232 + CDEN*CPART2*PI2*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG233 = CINTG233 + CDEN*CPART2*PI3*CCHI1(NCH1)*CCHI2(NCH2)
C
      ENDDO
      ENDDO
C
C
         CLDISP(1) = CINTG111+CINTG211
         CLDISP(2) = CINTG112+CINTG212
         CLDISP(3) = CINTG113+CINTG213
         CLDISP(4) = CINTG121+CINTG221
         CLDISP(5) = CINTG122+CINTG222
         CLDISP(6) = CINTG123+CINTG223
         CLDISP(7) = CINTG131+CINTG231
         CLDISP(8) = CINTG132+CINTG232
         CLDISP(9) = CINTG133+CINTG233
C 
      RETURN
      END                     
C
C ******************************************************
C                       DPDKDTI
C ******************************************************
C
      SUBROUTINE DPDKDTI(CW,CKPR,NS,NP,NTH,NR,NCH1,NCH2,MM,CLDISP2)
C
C------------------------------------------------------! 
C     NS : PARTICLE SPECIES
C     NP : MOMENTUM NUMBER
C     NTH : PITCH ANGLE MESH
C     NR : NODE NUMBER (RADIAL POSITION)
C     NCH1,NCH2: POLOIDAL MESH
C     MM : POLOIDAL NUMBER
C------------------------------------------------------! 

      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'

      DIMENSION CLDISP(9)
      DIMENSION DRFP1(NTHMAX,NPMAX,NRMAX)
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION DRRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DRCWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION DPFP1(NTHMAX,NPMAX,NRMAX)
      DIMENSION DRPFP(NTHMAX,NPMAX)
      DIMENSION DGP(NTHMAX,NPMAX)
      DIMENSION DGT(NTHMAX,NPMAX)
      DIMENSION TCSM2(NTHMAX)
      DIMENSION PV(NPMAX)
C
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      AA=AM/(AE*BABS*RR)
      RSL=DELR*(NR-0.5D0)


!--------derivation in radial direction & velocity---------------!
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP1(NTH,NP,NR) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         DRFP(NTH,NP) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         DRRFP(NTH,NP)=(DRFP1(NTH,NP,NR+1)-DRFP1(NTH,NP,NR-1))
     &        /(2*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH, NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP1(NTH,NP,NR) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=1,NPMAX-2
               DPFP1(NTH,NP,NR) = (FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/(2*DELP)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP1(NTH,NP,NR) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
         DRPFP(NTH,NP) = ((DPFP1(NTH,NP,NR+1)
     &        - DPFP1(NTH,NP,NR-1))/(2*DELP))/FP(NTH,NP,NR)
      ENDDO
C      
C
!------------------------------------------------------!
C
C***************SINGULAR POINT***************************
C
C
C*****************SUM3************************
C
      CINTG311 = (0.D0,0.D0)
      CINTG312 = (0.D0,0.D0)
      CINTG313 = (0.D0,0.D0)
      CINTG321 = (0.D0,0.D0)
      CINTG322 = (0.D0,0.D0)
      CINTG323 = (0.D0,0.D0)
      CINTG331 = (0.D0,0.D0)
      CINTG332 = (0.D0,0.D0)
      CINTG333 = (0.D0,0.D0)
C
C
      DO NTH=1,NTHMAX
C
         PV(NP)=PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
         D = DKPP**2*TCSM2(NTH)-(2*MM*AA*CW
     &        *(1+TCSM2(NTH))*CCHI2(NCH2))/RSL
C
C PNEAR1
C 
         PNEAR1=(DKPP*TCSM(NTH)+SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI2(NCH2)/RSL)
C
         CPART31=DPFP(NTH,NP)-CWAST(NTH,NP)/CW
C        
         CPART32=(MM*AA**2*PNEAR1**4*(1+TCSM2(NTH))**2*CCHI2(NCH2)
     &        *SCHI2(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/(4*CDENX))
C
         CPART33=((AA/2)*PNEAR1**2*(1+TCSM2(NTH))*SCHI2(NCH2))
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))
C
C F_0
C
         CPART34=CPART31+CPART32+CPART33  
C
C F_1
C
         CPART35=(AA/2)*PENAR1**2*(1+TCSM2(NTH))*SCHI2(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C
         CPART36=2*PI*CI*PNEAR1**3*TSNM(NTH)*COEF
     &        *(CPART34+CPART35)*DELTH
     &        /(ABS(-DKPP*TCSM(NTH)
     &        +MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI2(NCH2)))
C
         PI11=(AA/2)**2*PNEAR1**4*(1+TCSM2(NTH))**2
         PI12=(AA/2)*PNEAR1**4*(1+TCSM2(NTH))
         PI13=PNEAR1**4
C
C
         CINTG311 = CINTG311 + CPART36*PI11*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG312 = CINTG312 + CPART36*PI11*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG313 = CINTG313 + CPART36*PI12*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG321 = CINTG321 + CPART36*PI11*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG322 = CINTG322 + CPART36*PI11*SCHI1(NCH1)*SCHI2(NCH2)
         CINTG323 = CINTG323 + CPART36*PI12*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG331 = CINTG331 + CPART36*PI12*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG332 = CINTG332 + CPART36*PI12*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG333 = CINTG333 + CPART36*PI13*CCHI1(NCH1)*CCHI2(NCH2)
       ENDDO
C
C
C*****************SUM4************************
C
      CINTG411 = (0.D0,0.D0)
      CINTG412 = (0.D0,0.D0)
      CINTG413 = (0.D0,0.D0)
      CINTG421 = (0.D0,0.D0)
      CINTG422 = (0.D0,0.D0)
      CINTG423 = (0.D0,0.D0)
      CINTG431 = (0.D0,0.D0)
      CINTG432 = (0.D0,0.D0)
      CINTG433 = (0.D0,0.D0)
C
      DO NTH=1,NTHMAX
C  
         PV(NP)=PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
CPNEAR2
C
         PNEAR2 = (DKPP*TCSM(NTH)-SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI2(NCH2)/RSL)
C
         CPART41=DPFP(NTH,NP)-CWAST(NTH,NP)/CW
C     
         CPART42=(MM*AA**2*PNEAR2**4*(1+TCSM2(NTH))**2*CCHI2(NCH2)
     &        *SCHI2(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/(4*CDENX))
C
         CPART43=((AA/2)*PNEAR2**2*(1+TCSM2(NTH))*SCHI2(NCH2))
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))
C
C    F_0
C
         CPART44=(CPART41+CPART42+CPART43) 
C
C    F_1
C
         CPART45=(AA/2)*PENAR2**2*(1+TCSM2(NTH))*SCHI2(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C
         CPART46 =2*PI*CI*PNEAR2**3*TSNM(NTH)*COEF
     &        *(CPART44+CPART45)*DELTH
     &        /ABS(-DKPP*TCSM(NTH)
     &        +MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI2(NCH2))
C
C
         PI21=(AA/2)**2*PNEAR2**4*(1+TCSM2(NTH))**2
         PI22=(AA/2)*PNEAR2**4*(1+TCSM2(NTH))
         PI23=PNEAR2**4
C
C
         CINTG411 = CINTG411 + CPART46*PI21*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG412 = CINTG412 + CPART46*PI21*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG413 = CINTG413 + CPART46*PI22*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG421 = CINTG421 + CPART46*PI21*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG422 = CINTG422 + CPART46*PI21*SCHI1(NCH1)*SCHI2(NCH2)
         CINTG423 = CINTG423 + CPART46*PI22*SCHI1(NCH1)*CCHI2(NCH2)
         CINTG431 = CINTG431 + CPART46*PI22*CCHI1(NCH1)*CCHI2(NCH2)
         CINTG432 = CINTG432 + CPART46*PI22*CCHI1(NCH1)*SCHI2(NCH2)
         CINTG433 = CINTG433 + CPART46*PI23*CCHI1(NCH1)*CCHI2(NCH2)
      ENDDO
C
C
      CLDISP(1) = -CI*PI*(CINTG311+CINTG411)
      CLDISP(2) = -CI*PI*(CINTG312+CINTG412)
      CLDISP(3) = -CI*PI*(CINTG313+CINTG413)
      CLDISP(4) = -CI*PI*(CINTG321+CINTG421)
      CLDISP(5) = -CI*PI*(CINTG322+CINTG422)
      CLDISP(6) = -CI*PI*(CINTG323+CINTG423)
      CLDISP(7) = -CI*PI*(CINTG331+CINTG431)
      CLDISP(8) = -CI*PI*(CINTG332+CINTG432)
      CLDISP(9) = -CI*PI*(CINTG333+CINTG433)
C
C 
      RETURN
      END                     
