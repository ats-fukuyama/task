C
C ******************************************************
C                       DPDKDT
C      dielectric tensor with Drift Kinetic effects
C                 
C                      2010/02/08             
C                           programed by T.OKAMOTO
C ******************************************************
C
      SUBROUTINE DPDKDT(CW,CKPR,NS,NTH,NP,NR,NCH1,NCH2,MM,CLDISP)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION CLDISP(9),CLDISP1(9),CLDISP2(9)
C
      CALL DPDKDTR(CW,CKPR,NP,NR,NCH1,NCH2,MM,CLDISP1)
      CALL DPDKDTI(CW,CKPR,NP,NR,NCH1,NCH2,MM,CLDISP2)
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
      SUBROUTINE DPDKDTR(CW,CKPR,NS,NR,NCH1,NCH2,MM,CLDISP)
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
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION DRRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DRCWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION DPFPRP(NTHMAX,NPMAX)
      DIMENSION DPFPRM(NTHMAX,NPMAX)
      DIMENSION DRPFP(NTHMAX,NPMAX)
      DIMENSION DGP(NTHMAX,NPMAX)
      DIMENSION DGT(NTHMAX,NPMAX)
      DIMENSION TCSM2(NTHMAX)
      DIMENSION PV(NPMAX)
      REAL*8 CHIL,CCHIL
C
      DELPL  = 0.5D0
C      
      RSL=RM(NR)
      CHIL= 0.5D0*( CHI(NCH1)+ CHI(NCH2))
      CCHIL=0.5D0*(CCHI(NCH1)+CCHI(NCH2))
      SCHIL=0.5D0*(SCHI(NCH1)+SCHI(NCH2))
      X=RR+RSL*CCHIL
      Y=0.D0
      Z=   RSL*SCHIL
      CALL PLMAG(X,Y,Z,RHON)
C
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
      DKPP=MM/RSL
C
C
      DO NP=1,NPMAX
         PV(NP)=PTH0*PG(NP)/AM
      ENDDO
      DO NTH=1,NTHMAX
          TCSM2(NTH)=TCSM(NTH)**2
      ENDDO

!--------derivation in radial direction & velocity---------------!
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP(NTH,NP)=      (FP(NTH,NP,NR+1)-FP(NTH,NP,NR-1))
     &                     /(2*DELR)
         DRRFP(NTH,NP)=     (     FP(NTH,NP,NR+1)
     &                      -2.D0*FP(NTH,NP,NR  )
     &                      +     FP(NTH,NP,NR-1))
     &                     /(DELR*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
C
      DO NTH=1,NTHMAX
         NP=1
            DPFP(NTH,NP)   = (FP(NTH,NP+1,NR  )-FP(NTH,NP,NR  ))
     &                      /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH,NP) = (FP(NTH,NP+1,NR+1)-FP(NTH,NP,NR+1))
     &                      /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH,NP) = (FP(NTH,NP+1,NR-1)-FP(NTH,NP,NR-1))
     &                      /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR-1))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FP(NTH,NP+1,NR  )-FP(NTH,NP-1,NR  ))
     &                      /(2*PV(NP)*PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH,NP) = (FP(NTH,NP+1,NR+1)-FP(NTH,NP-1,NR+1))
     &                      /(2*PV(NP)*PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH,NP) = (FP(NTH,NP+1,NR-1)-FP(NTH,NP-1,NR-1))
     &                      /(2*PV(NP)*PTH0*DELP*FP(NTH,NP,NR-1))
         ENDDO
         NP=NPMAX-1
            DPFP(NTH, NP)   = (FP(NTH,NP,NR  )-FP(NTH,NP-1,NR  ))
     &                       /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH, NP) = (FP(NTH,NP,NR+1)-FP(NTH,NP-1,NR+1))
     &                       /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH, NP) = (FP(NTH,NP,NR-1)-FP(NTH,NP-1,NR-1))
     &                       /(PV(NP)*PTH0*DELP*FP(NTH,NP,NR-1))
      ENDDO
C
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            DRPFP(NTH,NP) = (DPFPRP(NTH,NP)-DPFPRM(NTH,NP))
     &                      /(2*DELR)
         ENDDO
      ENDDO
C           
C  
!--------------diffrerential of g(v,theta),DGP,DGT----------!
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX
         DGP(NTH,NP)=DKPP*PV(NP)*TSNM(NTH)
     &            -(MM*AA*PV(NP)**2*TCSM(NTH)*TSNM(NTH)*CCHI(NCH2))/RSL
         DGT(NTH,NP)=-DKPP*TCSM(NTH)
     &            +(MM*AA*PV(NP)*(1+TCSM2(NTH))*CCHI(NCH2))/RSL
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
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
C
         VD=AA*PV(NP)**2*0.5D0*(1+TCSM2(NTH))
         VPARA=PV(NP)*TCSM(NTH)
         COEF=AE**2*FP(NTH,NP,NR)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)/RSL
         CDEN=1.D0/CDENX
C
C         CDENY  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*PTH0*DELP)**2
C     &                             +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART11=(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C         CPART12=(0.D0,0.D0)
C         CPART13=(0.D0,0.D0)
C
         CPART12= CI*MM*VD**2*CCHI(NCH2)*SCHI(NCH2)
     &          *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN**2
C
         CPART13= CI*VD*SCHI(NCH2)
     &           *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &           +CWAST(NTH,NP)/(CW*RSL))*CDEN
C      
         CPART14= CI*2*PI*PM(NP)**2*TSNM(NTH)*DELTH*DELP
     &            *PN0*1.D20
     &            *COEF*(CPART11+CPART12+CPART13)
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VPARA*VPARA
C
C
      CINTG111 = CINTG111+CDEN*CPART14*PI1*CCHI(NCH1)*CCHI(NCH2)
      CINTG112 = CINTG112+CDEN*CPART14*PI1*CCHI(NCH1)*SCHI(NCH2)
      CINTG113 = CINTG113+CDEN*CPART14*PI2*CCHI(NCH1)
      CINTG121 = CINTG121+CDEN*CPART14*PI1*SCHI(NCH1)*CCHI(NCH2)
      CINTG122 = CINTG122+CDEN*CPART14*PI1*SCHI(NCH1)*SCHI(NCH2)
      CINTG123 = CINTG123+CDEN*CPART14*PI2*SCHI(NCH1)
      CINTG131 = CINTG131+CDEN*CPART14*PI2           *CCHI(NCH2)
      CINTG132 = CINTG132+CDEN*CPART14*PI2           *SCHI(NCH2)
      CINTG133 = CINTG133+CDEN*CPART14*PI3
C 
      ENDDO
      ENDDO
C
C
C
C*************SUM2********************
C
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
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         COEF=AE**2*FP(NTH,NP,NR)
C
         VD=AA*PV(NP)**2*0.5D0*(1+TCSM2(NTH))
         VPARA=PV(NP)*TCSM(NTH)
C     
         CDENX = CW-DKPRX+MM*VD*CCHI(NCH2)/RSL
     &        /(2*RSL)
         CDEN=1.D0/CDENX
C         CDENY  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*PTH0*DELP)**2
C     &                          +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART2= -PI*PV(NP)**4*COEF
     &        *TSNM(NTH)*AA*(1+TCSM2(NTH))*SCHI(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN
     &        *DELTH*DELP
C
         PI1=VD**2
         PI2=VD*VPARA
         PI3=VPARA**2
C
         CINTG211 = CINTG211+CDEN*CPART2*PI1*CCHI(NCH1)*CCHI(NCH2)
         CINTG212 = CINTG212+CDEN*CPART2*PI1*CCHI(NCH1)*SCHI(NCH2)
         CINTG213 = CINTG213+CDEN*CPART2*PI2*CCHI(NCH1)
         CINTG221 = CINTG221+CDEN*CPART2*PI1*SCHI(NCH1)*CCHI(NCH2)
         CINTG222 = CINTG222+CDEN*CPART2*PI1*SCHI(NCH1)*SCHI(NCH2)
         CINTG223 = CINTG223+CDEN*CPART2*PI2*SCHI(NCH1)
         CINTG231 = CINTG231+CDEN*CPART2*PI2           *CCHI(NCH2)
         CINTG232 = CINTG232+CDEN*CPART2*PI2           *SCHI(NCH2)
         CINTG233 = CINTG233+CDEN*CPART2*PI3
      ENDDO
      ENDDO
C
C
C****************SUM3*********************************
C
      CINTG312a = (0.D0,0.D0)
      CINTG313a = (0.D0,0.D0)
      CINTG322a = (0.D0,0.D0)
      CINTG323a = (0.D0,0.D0)
      CINTG332a = (0.D0,0.D0)
      CINTG333a = (0.D0,0.D0)
C

      CINTG312 = (0.D0,0.D0)
      CINTG313 = (0.D0,0.D0)
      CINTG322 = (0.D0,0.D0)
      CINTG323 = (0.D0,0.D0)
      CINTG332 = (0.D0,0.D0)
      CINTG333 = (0.D0,0.D0)
C
      DO NP=1, NPMAX-1
      DO NTH=1, NTHMAX
C
         EPART1 = DRFP(NTH,NP)*PV(NP)**4*(1+TCSM2(NTH))*TSNM(NTH)
         EPART2 = EPART1
         EPART3 = DRFP(NTH,NP)*PV(NP)**3*TCSM(NTH)*TSNM(NTH)
C
      CINTG312a = CINTG312a+EPART1*DELTH*DELP
      CINTG313a = CINTG313a+EPART1*DELTH*DELP
      CINTG322a = CINTG322a+EPART2*DELTH*DELP
      CINTG323a = CINTG323a+EPART2*DELTH*DELP
      CINTG332a = CINTG332a+EPART3*DELTH*DELP
      CINTG333a = CINTG333a+EPART3*DELTH*DELP
C
      ENDDO
      ENDDO
C
      COEE1 = -CI*PI*AE*AA/(CW*BABS)
      COEE2 = COEE1
      COEE3 = -CI*2*PI*AE/(CW*BABS)
C
      CINTG312 = COEE1*CINTG312a*SCHI(NCH1)
      CINTG313 = COEE1*CINTG313a*SCHI(NCH1)
      CINTG322 = COEE2*CINTG322a*CCHI(NCH1)
      CINTG323 = COEE2*CINTG323a*CCHI(NCH1)
      CINTG332 = COEE3*CINTG332a
      CINTG333 = COEE3*CINTG333a
C
C****************SUM4*********************************
C
      CINTG412 = (0.D0,0.D0)
      CINTG421 = (0.D0,0.D0)
C
      DO NP=1, NPMAX-1
      DO NTH=1, NTHMAX
C
         FPART= 2*PI*PN0*1.D20*FP(NTH,NP,NR)*PM(NP)**2*TSNM(NTH)
C
      CINTG412 = CINTG412+FPART*DELTH*DELP
      CINTG421 = CINTG421+FPART*DELTH*DELP
C
      ENDDO
      ENDDO
C
      COEE4 = AE/BABS
C
      CINTG412 = COEE4*CINTG412
      CINTG421 = COEE4*CINTG421
C
C*****************************************************
C
      CLDISP(1) = CINTG111+CINTG211
      CLDISP(2) = CINTG112+CINTG212+CINTG312+CINTG412
      CLDISP(3) = CINTG113+CINTG213+CINTG313
      CLDISP(4) = CINTG121+CINTG221         +CINTG421
      CLDISP(5) = CINTG122+CINTG222+CINTG322
      CLDISP(6) = CINTG123+CINTG223+CINTG323
      CLDISP(7) = CINTG131+CINTG231
      CLDISP(8) = CINTG132+CINTG232+CINTG332
      CLDISP(9) = CINTG133+CINTG233+CINTG333
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
      RSL=RM(NR)
      CHIL=0.5D0*(CHI(NCH1)+CHI(NCH2))
      CCHIL=0.5D0*(CCHI(NCH1)+CCHI(NCH2))
      SCHIL=0.5D0*(SCHI(NCH1)+SCHI(NCH2))
      X=RR+RSL*CCHIL
      Y=0.D0
      Z=   RSL*SCHIL
      CALL PLMAG(X,Y,Z,RHON)
C
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
      DKPP=MM/RSL
C
      CWP=PN0**1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
C
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
            DPFP(NTH,NP) = (FP(NTH,NP+1,NR)-FP(NTH,NP,NR))
     &                    /(PTH0*DELP*FP(NTH,NP,NR))
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = (FP(NTH,NP+1,NR)- FP(NTH,NP-1,NR))
     &                       /(2*PHT0*DELP*FP(NTH,NP,NR))
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH, NP) = (FP(NTH,NP,NR)-FP(NTH,NP-1,NR))
     &                     /(PTH0*DELP*FP(NTH,NP,NR))
         END IF
      ENDDO
C
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP1(NTH,NP,NR) = (FP(NTH,NP+1,NR)-FP(NTH,NP,NR))
     &                        /(PTH0*DELP*FP(NTH,NP,NR))
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=1,NPMAX-2
               DPFP1(NTH,NP,NR) = (FP(NTH,NP+1,NR)-FP(NTH,NP-1,NR))
     &                           /(2*PTH0*DELP*FP(NTH,NP,NR))
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP1(NTH,NP,NR) = (FP(NTH,NP,NR)-FP(NTH,NP-1,NR))
     &                        /(PTH0*DELP*FP(NTH,NP,NR))
         END IF
         DRPFP(NTH,NP) = (DPFP1(NTH,NP,NR+1)-DPFP1(NTH,NP,NR-1))
     &                  /(2*PTH0*DELP*FP(NTH,NP,NR))
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
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
      CINTG311a = (0.D0,0.D0)
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
         PV(NP)=PTH0*PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
C
         D = DKPP**2*TCSM2(NTH)-(2*MM*AA*CW
     &        *(1+TCSM2(NTH))*CCHI(NCH2))/RSL
C
C PNEAR1
C 
         PNEAR1=(DKPP*TCSM(NTH)+SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI(NCH2)/RSL)
C
         CPART31=DPFP(NTH,NP)-CWAST(NTH,NP)/CW
C        
         CPART32=(MM*AA**2*PNEAR1**4*(1+TCSM2(NTH))**2*CCHI(NCH2)
     &        *SCHI(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/(4*CDENX))
C
         CPART33=((AA/2)*PNEAR1**2*(1+TCSM2(NTH))*SCHI(NCH2))
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))
C
C F_0
C
         CPART34=CPART31+CPART32+CPART33  
C
C COE1=MM*AA**2
C
C F_1
C
         CPART35=(AA/2)*PENAR1**2*(1+TCSM2(NTH))*SCHI(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C
         CPART36=2*PI*CI*PNEAR1**3*TSNM(NTH)*COEF
     &        *(CPART34+CPART35)*DELTH
     &        /(ABS(-DKPP*TCSM(NTH)
     &        +MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2)))
C
         VD=PNEAR1**2*(1+TCSM2(NTH))
         VPARA=PNEAR1*TCSM(NTH)
C
         PI11=VD*VD
         PI12=VD*VPARA
         PI13=VD*VPARA
C
C
         CINTG311a = CINTG311a + CPART36*PI11
         CINTG312a = CINTG312a + CPART36*PI11
         CINTG313a = CINTG313a + CPART36*PI12
         CINTG321a = CINTG321a + CPART36*PI11
         CINTG322a = CINTG322a + CPART36*PI11
         CINTG323a = CINTG323a + CPART36*PI12
         CINTG331a = CINTG331a + CPART36*PI12
         CINTG332a = CINTG332a + CPART36*PI12
         CINTG333a = CINTG333a + CPART36*PI13
       ENDDO
C
       CINTG311 = CINTG311a
       CINTG312 = CINTG312a
       CINTG313 = CINTG313a
       CINTG321 = CINTG321a
       CINTG322 = CINTG322a
       CINTG323 = CINTG323a
       CINTG331 = CINTG331a
       CINTG332 = CINTG332a
       CINTG333 = CINTG333a
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
         PV(NP)=PTH0*PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
C
CPNEAR2
C
         PNEAR2 = (DKPP*TCSM(NTH)-SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI(NCH2)/RSL)
C
         CPART41=DPFP(NTH,NP)-CWAST(NTH,NP)/CW
C     
         CPART42=(MM*AA**2*PNEAR2**4*(1+TCSM2(NTH))**2*CCHI(NCH2)
     &        *SCHI(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/(4*CDENX))
C
         CPART43=((AA/2)*PNEAR2**2*(1+TCSM2(NTH))*SCHI(NCH2))
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))
C
C    F_0
C
         CPART44=(CPART41+CPART42+CPART43) 
C
C    F_1
C
         CPART45=(AA/2)*PENAR2**2*(1+TCSM2(NTH))*SCHI(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C
         CPART46 =2*PI*CI*PNEAR2**3*TSNM(NTH)*COEF
     &        *(CPART44+CPART45)*DELTH
     &        /ABS(-DKPP*TCSM(NTH)
     &        +MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2))
C
C
         PI21=(AA/2)**2*PNEAR2**4*(1+TCSM2(NTH))**2
         PI22=(AA/2)*PNEAR2**4*(1+TCSM2(NTH))
         PI23=PNEAR2**4
C
C
         CINTG411 = CINTG411 + CPART46*PI21
         CINTG412 = CINTG412 + CPART46*PI21
         CINTG413 = CINTG413 + CPART46*PI22
         CINTG421 = CINTG421 + CPART46*PI21
         CINTG422 = CINTG422 + CPART46*PI21
         CINTG423 = CINTG423 + CPART46*PI22
         CINTG431 = CINTG431 + CPART46*PI22
         CINTG432 = CINTG432 + CPART46*PI22
         CINTG433 = CINTG433 + CPART46*PI23
      ENDDO
C
      FACT = 2.0D0*PI*DBLE(CWP)
C
      CLDISP(1) = -CI*PI*(CINTG311+CINTG411)*CCHI(NCH1)*CCHI(NCH2)*FACT
      CLDISP(2) = -CI*PI*(CINTG312+CINTG412)*CCHI(NCH1)*SCHI(NCH2)*FACT
      CLDISP(3) = -CI*PI*(CINTG313+CINTG413)*CCHI(NCH1)           *FACT 
      CLDISP(4) = -CI*PI*(CINTG321+CINTG421)*SCHI(NCH1)*CCHI(NCH2)*FACT
      CLDISP(5) = -CI*PI*(CINTG322+CINTG422)*SCHI(NCH1)*SCHI(NCH2)*FACT
      CLDISP(6) = -CI*PI*(CINTG323+CINTG423)*SCHI(NCH1)           *FACT
      CLDISP(7) = -CI*PI*(CINTG331+CINTG431)           *CCHI(NCH2)*FACT
      CLDISP(8) = -CI*PI*(CINTG332+CINTG432)           *SCHI(NCH2)*FACT
      CLDISP(9) = -CI*PI*(CINTG333+CINTG433)                      *FACT
C
C 
      RETURN
      END                     
