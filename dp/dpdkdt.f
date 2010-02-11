C
C ******************************************************
C                       DPDKDT
C      dielectric tensor with Drift Kinetic effects
C                 
C                      2010/02/08             
C                           programed by T.OKAMOTO
C ******************************************************
C
      SUBROUTINE DPDKDT(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNS)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION CDTNS(3,3),CDTNSR(3,3),CDTNSI(3,3)
C
      CALL DPDKDTR(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNSR)
C      CALL DPDKDTI(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNSI)
      DO J=1,3
         DO I=1,3
            CDTNS(I,J)=CDTNSR(I,J)+CDTNSI(I,J)
         ENDDO
      ENDDO
      RETURN
      END
C
C ******************************************************
C                       DPDKDTR
C ******************************************************
C
      SUBROUTINE DPDKDTR(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNSR)
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
      DIMENSION CDTNSR(3,3)
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

      DELPL  = 0.5D0
      
      PN0 = PN(NS)
      PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

      RSL=MAX(RM(NR),1.D-6)
      CHIL= 0.5D0*( CHI(NCH1)+ CHI(NCH2))
      CCHIL=0.5D0*(CCHI(NCH1)+CCHI(NCH2))
      SCHIL=0.5D0*(SCHI(NCH1)+SCHI(NCH2))
      X=RR+RSL*CCHIL
      Y=0.D0
      Z=   RSL*SCHIL
      CALL PLMAG(X,Y,Z,RHON)
C
      CALL GUFLSH
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
      DKPP=MM/RSL
C
      DO NP=1,NPMAX
         PV(NP)=PTH0*PG(NP)/AM
      ENDDO
      DO NTH=1,NTHMAX
          TCSM2(NTH)=TCSM(NTH)**2
      ENDDO

!--------derivation in radial direction & velocity---------------!
C
      IF(NR.EQ.1) THEN
         NRML=NR
         NRPL=NR+1
         DELRL=DELR
      ELSEIF(NR.EQ.NRMAX) THEN
         NRML=NR-1
         NRPL=NR
         DELRL=DELR
      ELSE
         NRML=NR-1
         NRPL=NR+1
         DELRL=2.D0*DELR
      ENDIF
      DO NP=1,NPMAX-1
      DO NTH=1,NTHMAX
         DRFP(NTH,NP)=      (FNS(NTH,NP,NRPL,NS)-FNS(NTH,NP,NRML,NS))
     &                     /(DELRL)
C         DRRFP(NTH,NP)=    ((FNS(NTH,NP,NRPL,NS)-FNS(NTH,NP,NR  ,NS))
C     &                     /DELR
C     &                     -(FNS(NTH,NP,NR  ,NS)-FNS(NTH,NP,NRML,NS))
C     &                     /DELR)/DELR
         CWAST(NTH,NP)  =MM*(DRFP(NTH,NP)
     &                  /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
C         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)
C     &                  /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
      ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         NP=1
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)-FNS(NTH,NP,NR,  NS))
     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,NS))
C            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NRPL,NS)-FNS(NTH,NP,NRPL,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRPL,NS))
C            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NRML,NS)-FNS(NTH,NP,NRML,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRML,NS))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)
     &                       -FNS(NTH,NP-1,NR,  NS))
     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,NS))
C            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NRPL,NS)
C     &                       -FNS(NTH,NP-1,NRPL,NS))
C     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NRPL,NS))
C            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NRML,NS)
C     &                       -FNS(NTH,NP-1,NRML,NS))
C     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NRML,NS))
         ENDDO
         NP=NPMAX-1
            DPFP(NTH, NP)   = (FNS(NTH,NP,  NR,  NS)
     &                        -FNS(NTH,NP-1,NR,  NS))
     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,  NS))
C            DPFPRP(NTH, NP) = (FNS(NTH,NP,  NRPL,NS)
C     &                        -FNS(NTH,NP-1,NRPL,NS))
C     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRPL,NS))
C            DPFPRM(NTH, NP) = (FNS(NTH,NP,  NRML,NS)
C     &                        -FNS(NTH,NP-1,NRML,NS))
C     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRML,NS))
      ENDDO

C      DO NP=1,NPMAX
C         DO NTH=1,NTHMAX
C            DRPFP(NTH,NP) = (DPFPRP(NTH,NP)-DPFPRM(NTH,NP))
C     &                      /(2*DELR)
C         ENDDO
C      ENDDO
           
!--------------diffrerential of g(v,theta),DGP,DGT----------!

C      DO NTH=1,NTHMAX
C      DO NP=1,NPMAX
C         DGP(NTH,NP)=DKPP*PV(NP)*TSNM(NTH)
C     &            -(MM*AA*PV(NP)**2*TCSM(NTH)*TSNM(NTH)*CCHI(NCH2))/RSL
C         DGT(NTH,NP)=-DKPP*TCSM(NTH)
C     &            +(MM*AA*PV(NP)*(1+TCSM2(NTH))*CCHI(NCH2))/RSL
C      ENDDO
C      ENDDO
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
         COEF=AE**2*FNS(NTH,NP,NR,NS)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)/RSL
         CDEN=1.D0/CDENX
C
C         CDENY  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*PTH0*DELP)**2
C     &                             +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART11=(-DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
         CPART12=(0.D0,0.D0)
         CPART13=(0.D0,0.D0)
         CPART14=(0.D0,0.D0)
C
C         CPART12= CI*MM*VD**2*CCHI(NCH2)*SCHI(NCH2)
C     &          *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN**2
C
C         CPART13= CI*VD*SCHI(NCH2)
C     &           *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
C     &           +CWAST(NTH,NP)/(CW*RSL))*CDEN
C      
C         CPART14= -PI*PV(NP)**4*COEF
C     &        *TSNM(NTH)*AA*(1+TCSM2(NTH))*SCHI(NCH2)
C     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN
C     &        *DELTH*DELP

         CPART1= CI*2*PI*PM(NP)**2*TSNM(NTH)*DELTH*DELP
     &            *PN0*1.D20
     &            *COEF*(CPART11+CPART12+CPART13+CPART14)
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VPARA*VPARA
C
C
      CINTG111 = CINTG111+CDEN*CPART1*PI1*CCHI(NCH1)*CCHI(NCH2)
      CINTG112 = CINTG112+CDEN*CPART1*PI1*CCHI(NCH1)*SCHI(NCH2)
      CINTG113 = CINTG113+CDEN*CPART1*PI2*CCHI(NCH1)
      CINTG121 = CINTG121+CDEN*CPART1*PI1*SCHI(NCH1)*CCHI(NCH2)
      CINTG122 = CINTG122+CDEN*CPART1*PI1*SCHI(NCH1)*SCHI(NCH2)
      CINTG123 = CINTG123+CDEN*CPART1*PI2*SCHI(NCH1)
      CINTG131 = CINTG131+CDEN*CPART1*PI2           *CCHI(NCH2)
      CINTG132 = CINTG132+CDEN*CPART1*PI2           *SCHI(NCH2)
      CINTG133 = CINTG133+CDEN*CPART1*PI3
C 
      ENDDO
      ENDDO
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
      COEE2 =  COEE1
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
         FPART= 2*PI*PN0*1.D20*FNS(NTH,NP,NR,NS)*PM(NP)**2*TSNM(NTH)
C
         CINTG412 = CINTG412+FPART*DELTH*DELP
         CINTG421 = CINTG421-FPART*DELTH*DELP
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
      CDTNSR(1,1)=CI*(CINTG111                  )/(CW*EPS0)
      CDTNSR(1,2)=CI*(CINTG112+CINTG312+CINTG412)/(CW*EPS0)
      CDTNSR(1,3)=CI*(CINTG113+CINTG313         )/(CW*EPS0)
      CDTNSR(2,1)=CI*(CINTG121         +CINTG421)/(CW*EPS0)
      CDTNSR(2,2)=CI*(CINTG122+CINTG322         )/(CW*EPS0)
      CDTNSR(2,3)=CI*(CINTG123+CINTG323         )/(CW*EPS0)
      CDTNSR(3,1)=CI*(CINTG131                  )/(CW*EPS0)
      CDTNSR(3,2)=CI*(CINTG132+CINTG332         )/(CW*EPS0)
      CDTNSR(3,3)=CI*(CINTG133+CINTG333         )/(CW*EPS0)
C
      RETURN
      END                     
C
C ******************************************************
C                       DPDKDTI
C ******************************************************
C
      SUBROUTINE DPDKDTI(CW,CKPR,NS,NP,NTH,NR,NCH1,NCH2,MM,CDTNSI)
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

      DIMENSION CDTNSI(3,3)
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION DRRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DRCWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION DPFPRP(NTHMAX,NPMAX)
      DIMENSION DPFPRM(NTHMAX,NPMAX)
      DIMENSION DRPFP(NTHMAX,NPMAX)
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
      DO NP=1,NPMAX
         PV(NP)=PTH0*PG(NP)/AM
      ENDDO
      DO NTH=1,NTHMAX
          TCSM2(NTH)=TCSM(NTH)**2
      ENDDO

      CWP=PN0**1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
!--------derivation in radial direction & velocity---------------!
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP(NTH,NP)=      (FNS(NTH,NP,NR+1,NS)-FNS(NTH,NP,NR-1,NS))
     &                     /(2*DELR)
         DRRFP(NTH,NP)=     (     FNS(NTH,NP,NR+1,NS)
     &                      -2.D0*FNS(NTH,NP,NR  ,NS )
     &                      +     FNS(NTH,NP,NR-1,nS))
     &                     /(DELR*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)
     &                /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)
     &                  /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         NP=1
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)-FNS(NTH,NP,NR,  NS))
     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,NS))
            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NR+1,NS)-FNS(NTH,NP,NR+1,NS))
     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NR-1,NS)-FNS(NTH,NP,NR-1,NS))
     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)
     &                       -FNS(NTH,NP-1,NR,  NS))
     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,NS))
            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NR+1,NS)
     &                       -FNS(NTH,NP-1,NR+1,NS))
     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NR-1,NS)
     &                       -FNS(NTH,NP-1,NR-1,NS))
     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
         ENDDO
         NP=NPMAX-1
            DPFP(NTH, NP)   = (FNS(NTH,NP,  NR,  NS)
     &                        -FNS(NTH,NP-1,NR,  NS))
     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR,  NS))
            DPFPRP(NTH, NP) = (FNS(NTH,NP,  NR+1,NS)
     &                        -FNS(NTH,NP-1,NR+1,NS))
     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
            DPFPRM(NTH, NP) = (FNS(NTH,NP,  NR-1,NS)
     &                        -FNS(NTH,NP-1,NR-1,NS))
     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
      ENDDO
C
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            DRPFP(NTH,NP) = (DPFPRP(NTH,NP)-DPFPRM(NTH,NP))
     &                      /(2*DELR)
         ENDDO
      ENDDO
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

         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         COEF=AE**2*FNS(NTH,NP,NR,NS)
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
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         COEF=AE**2*FNS(NTH,NP,NR,NS)
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
      FACT = 2.0D0*PI*DBLE(CWP)/(CW*EPS0)
C
      CDTNSI(1,1)=PI*(CINTG311+CINTG411)*CCHI(NCH1)*CCHI(NCH2)*FACT
      CDTNSI(1,2)=PI*(CINTG312+CINTG412)*CCHI(NCH1)*SCHI(NCH2)*FACT
      CDTNSI(1,3)=PI*(CINTG313+CINTG413)*CCHI(NCH1)           *FACT
      CDTNSI(2,1)=PI*(CINTG321+CINTG421)*SCHI(NCH1)*CCHI(NCH2)*FACT
      CDTNSI(2,2)=PI*(CINTG322+CINTG422)*SCHI(NCH1)*SCHI(NCH2)*FACT
      CDTNSI(2,3)=PI*(CINTG323+CINTG423)*SCHI(NCH1)           *FACT
      CDTNSI(3,1)=PI*(CINTG331+CINTG431)           *CCHI(NCH2)*FACT
      CDTNSI(3,2)=PI*(CINTG332+CINTG432)           *SCHI(NCH2)*FACT
      CDTNSI(3,3)=PI*(CINTG333+CINTG433)                      *FACT
C
C 
      RETURN
      END                     
