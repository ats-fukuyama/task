C
C ******************************************************
C                       DPDKDT
C      dielectric tensor with Drift Kinetic effects
C                 
C                      2010/02/16             
C                           programed by T.OKAMOTO
C ******************************************************
C
      SUBROUTINE DPDKDT(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNS)
C
      USE plcomm
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

      USE plcomm
      USE pllocal
      USE plprof,ONLY:pl_mag_old
      INCLUDE 'dpcomm.inc'
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
      CALL PL_MAG_OLD(X,Y,Z,RHON)
C
      CALL GUFLSH
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
      DKPP=MM/RSL
C
      DO NP=1,NPMAX
         PV(NP)=PTH0*PG(NP,NS)/AM
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
     &                      /(PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,NS))
C            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NRPL,NS)-FNS(NTH,NP,NRPL,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRPL,NS))
C            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NRML,NS)-FNS(NTH,NP,NRML,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NRML,NS))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)
     &                       -FNS(NTH,NP-1,NR,  NS))
     &                      /(2*PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,NS))
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
     &                       /(PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,  NS))
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
C
         CPART12= CI*MM*VD**2*CCHI(NCH2)*SCHI(NCH2)
     &          *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN**2
C
         CPART13= CI*VD*SCHI(NCH2)
     &           *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &           +CWAST(NTH,NP)/(CW*RSL))*CDEN
C
         CPART1= CI*2*PI*PM(NP,NS)**2*TSNM(NTH)*DELTH*DELP(NS)
     &            *PN0*1.D20
     &            *COEF*(CPART11+CPART12+CPART13)
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
      CINTG312a = CINTG312a+EPART1*DELTH*DELP(NS)
      CINTG313a = CINTG313a+EPART1*DELTH*DELP(NS)
      CINTG322a = CINTG322a+EPART2*DELTH*DELP(NS)
      CINTG323a = CINTG323a+EPART2*DELTH*DELP(NS)
      CINTG332a = CINTG332a+EPART3*DELTH*DELP(NS)
      CINTG333a = CINTG333a+EPART3*DELTH*DELP(NS)
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
         FPART= 2*PI*PN0*1.D20*FNS(NTH,NP,NR,NS)*PM(NP,NS)**2*TSNM(NTH)
C
         CINTG412 = CINTG412+FPART*DELTH*DELP(NS)
         CINTG421 = CINTG421-FPART*DELTH*DELP(NS)
C
      ENDDO
      ENDDO
C
      COEE4 = AE/BABS
C
      CINTG412 = COEE4*CINTG412
      CINTG421 = COEE4*CINTG421
C
C************SUM 5 : Polarization*****************************************
C
      CINTG511 = (0.D0,0.D0)
      CINTG522 = (0.D0,0.D0)
C
      DO NP=1, NPMAX-1
      DO NTH=1, NTHMAX
C
         FPART= 2*PI*PN0*1.D20*FNS(NTH,NP,NR,NS)*PM(NP,NS)**2*TSNM(NTH)
C
      CINTG511 = CINTG511+FPART*DELTH*DELP(NS)
      CINTG522 = CINTG522+FPART*DELTH*DELP(NS)
C
      ENDDO
      ENDDO
C
      COEE5 = AM*CW/BABS**2
C       
      CINTG511 = -CI*COEE5*CINTG511
      CINTG522 = -CI*COEE5*CINTG522
C
C*****************************************************
C
      CDTNSR(1,1)=CI*(CINTG111                  +CINTG511)/(CW*EPS0)
      CDTNSR(1,2)=CI*(CINTG112+CINTG312+CINTG412         )/(CW*EPS0)
      CDTNSR(1,3)=CI*(CINTG113+CINTG313                  )/(CW*EPS0)
      CDTNSR(2,1)=CI*(CINTG121         +CINTG421         )/(CW*EPS0)
      CDTNSR(2,2)=CI*(CINTG122+CINTG322         +CINTG522)/(CW*EPS0)
      CDTNSR(2,3)=CI*(CINTG123+CINTG323                  )/(CW*EPS0)
      CDTNSR(3,1)=CI*(CINTG131                           )/(CW*EPS0)
      CDTNSR(3,2)=CI*(CINTG132+CINTG332                  )/(CW*EPS0)
      CDTNSR(3,3)=CI*(CINTG133+CINTG333                  )/(CW*EPS0)
C
      RETURN
      END                     
C
C ******************************************************
C                       DPDKDTI
C ******************************************************
C
      SUBROUTINE DPDKDTI(CW,CKPR,NS,NR,NCH1,NCH2,MM,CDTNSI)
C
C------------------------------------------------------! 
C     NS : PARTICLE SPECIES
C     NP : MOMENTUM NUMBER
C     NTH : PITCH ANGLE MESH
C     NR : NODE NUMBER (RADIAL POSITION)
C     NCH1,NCH2: POLOIDAL MESH
C     MM : POLOIDAL NUMBER
C------------------------------------------------------! 

      USE plcomm
      USE pllocal
      USE plprof,ONLY: pl_mag_old
      INCLUDE 'dpcomm.inc'

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
      DELPL = 0.5D0
      PN0 = PN(NS)
      PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PN(NS))
      RSL=MAX(RM(NR),1.D0-6)
      CHIL=0.5D0*(CHI(NCH1)+CHI(NCH2))
      CCHIL=0.5D0*(CCHI(NCH1)+CCHI(NCH2))
      SCHIL=0.5D0*(SCHI(NCH1)+SCHI(NCH2))
      X=RR+RSL*CCHIL
      Y=0.D0
      Z=   RSL*SCHIL
      CALL PL_MAG_OLD(X,Y,Z,RHON)
C
      CALL GUFLSH
      AM=PA(NS)*AMP             
      AE=PZ(NS)*AEE             
      AA=AM/(AE*BABS*RR)
      DKPP=MM/RSL
C
      DO NP=1,NPMAX
         PV(NP)=PTH0*PG(NP,NS)/AM
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
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP(NTH,NP)=      (FNS(NTH,NP,NRPL,NS)-FNS(NTH,NP,NRML,NS))
     &                     /(DELR)
C         DRRFP(NTH,NP)=     (     FNS(NTH,NP,NR+1,NS)
C     &                      -2.D0*FNS(NTH,NP,NR  ,NS )
C     &                      +     FNS(NTH,NP,NR-1,nS))
C     &                     /(DELR*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)
     &                /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
C         DRCWAST(NTH,NP)=MM*(DRRFP(NTH,NP)
C     &                  /FNS(NTH,NP,NR,NS))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         NP=1
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)-FNS(NTH,NP,NR,  NS))
     &                      /(PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,NS))
C            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NR+1,NS)-FNS(NTH,NP,NR+1,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
C            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NR-1,NS)-FNS(NTH,NP,NR-1,NS))
C     &                      /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FNS(NTH,NP+1,NR,  NS)
     &                       -FNS(NTH,NP-1,NR,  NS))
     &                      /(2*PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,NS))
C            DPFPRP(NTH,NP) = (FNS(NTH,NP+1,NR+1,NS)
C     &                       -FNS(NTH,NP-1,NR+1,NS))
C     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
C            DPFPRM(NTH,NP) = (FNS(NTH,NP+1,NR-1,NS)
C     &                       -FNS(NTH,NP-1,NR-1,NS))
C     &                      /(2*PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
         ENDDO
C
        NP=NPMAX-1
            DPFP(NTH, NP)   = (FNS(NTH,NP,  NR,  NS)
     &                        -FNS(NTH,NP-1,NR,  NS))
     &                       /(PV(NP)*PTH0*DELP(NS)*FNS(NTH,NP,NR,  NS))
C            DPFPRP(NTH, NP) = (FNS(NTH,NP,  NR+1,NS)
C     &                        -FNS(NTH,NP-1,NR+1,NS))
C     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR+1,NS))
C            DPFPRM(NTH, NP) = (FNS(NTH,NP,  NR-1,NS)
C     &                        -FNS(NTH,NP-1,NR-1,NS))
C     &                       /(PV(NP)*PTH0*DELP*FNS(NTH,NP,NR-1,NS))
      ENDDO
C
C      DO NP=1,NPMAX
C         DO NTH=1,NTHMAX
C            DRPFP(NTH,NP) = (DPFPRP(NTH,NP)-DPFPRM(NTH,NP))
C     &                      /(2*DELR)
C         ENDDO
C      ENDDO
CC
!------------------------------------------------------!
C
C***************SINGULAR POINT***************************
C
C
C*****************SUM1************************
C
      CINTG611 = (0.D0,0.D0)
      CINTG612 = (0.D0,0.D0)
      CINTG613 = (0.D0,0.D0)
      CINTG621 = (0.D0,0.D0)
      CINTG622 = (0.D0,0.D0)
      CINTG623 = (0.D0,0.D0)
      CINTG631 = (0.D0,0.D0)
      CINTG632 = (0.D0,0.D0)
      CINTG633 = (0.D0,0.D0)
C
C
      DO NTH=1,NTHMAX
C
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*TCSM(NTH) 
         DKPRY=DBLE(CKPRY)
         COEF=AE**2*FNS(NTH,NP,NR,NS)
C
         D = DKPRY**2*TCSM2(NTH)-(2*MM*AA*CW
     &        *(1+TCSM2(NTH))*CCHI(NCH2))/RSL
C
C PNEAR1
C 
         PNEAR1=(DKPRY+SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI(NCH2)/RSL)

         VD=AA*PNEAR1**2*0.5D0*(1+TCSM2(NTH))
         VPARA=PNEAR1*TCSM(NTH)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)/RSL
         CDEN=1.D0/CDENX
C
C F_0
C
         CPART31=(-DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C        
         CPART32=CI*MM*VD**2*CCHI(NCH2)*SCHI(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN**2
C
         CPART33=CI*VD*SCHI(NCH2)
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))*CDEN
C
         CPART34=CI*2*PI*PNEAR1**2*TSNM(NTH)*COEF
     &        *(CPART31+CPART32+CPART33)*DELTH*PN0*1.D20
     &        /(ABS(-DKPRY+MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2)))
C
C F_1
         CPART35=VD*SCHI(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN
C
         CPART36=-PI*COEF*CPART35*PN0*1.D20*DELTH
     &        /(ABS(-DKPRY+MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2)))
C
C F_0+F_1
         CPART37=CPART34+CPART36
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VD*VPARA
C
         CINTG611 = CINTG611 + CPART37*PI1*CCHI(NCH1)*CCHI(NCH2)
         CINTG612 = CINTG612 + CPART37*PI1*CCHI(NCH1)*SCHI(NCH2)
         CINTG613 = CINTG613 + CPART37*PI2*CCHI(NCH1)
         CINTG621 = CINTG621 + CPART37*PI1*SCHI(NCH1)*CCHI(NCH2)
         CINTG622 = CINTG622 + CPART37*PI1*SCHI(NCH1)*SCHI(NCH2)
         CINTG623 = CINTG623 + CPART37*PI2*SCHI(NCH1)
         CINTG631 = CINTG631 + CPART37*PI2           *CCHI(NCH2)
         CINTG632 = CINTG632 + CPART37*PI2           *SCHI(NCH2)
         CINTG633 = CINTG633 + CPART37*PI3
       ENDDO
C
C
C
C*****************SUM2************************
C
      CINTG711 = (0.D0,0.D0)
      CINTG712 = (0.D0,0.D0)
      CINTG713 = (0.D0,0.D0)
      CINTG721 = (0.D0,0.D0)
      CINTG722 = (0.D0,0.D0)
      CINTG723 = (0.D0,0.D0)
      CINTG731 = (0.D0,0.D0)
      CINTG732 = (0.D0,0.D0)
      CINTG733 = (0.D0,0.D0)
C
      DO NTH=1,NTHMAX
C
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*TCSM(NTH)
         DKPRY=DBLE(CKPRY)
         COEF=AE**2*FNS(NTH,NP,NR,NS)
C
C PNEAR2
C
C
         PNEAR2=(DKPRY-SQRT(D))/(MM*AA*(1+TCSM2(NTH))
     &        *CCHI(NCH2)/RSL)

         VD=AA*PNEAR1**2*0.5D0*(1+TCSM2(NTH))
         VPARA=PNEAR1*TCSM(NTH)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)/RSL
         CDEN=1.D0/CDENX
C
C F_0
         CPART41=(-DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C
         CPART42=CI*MM*VD**2*CCHI(NCH2)*SCHI(NCH2)
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN**2
C
         CPART43=CI*VD*SCHI(NCH2)
     &        *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &        +CWAST(NTH,NP)/(CW*RSL))*CDEN
C
         CPART44=CI*2*PI*PNEAR1**2*TSNM(NTH)*COEF
     &        *(CPART41+CPART42+CPART43)*DELTH*PN0*1.D20
     &        /(ABS(-DKPRY+MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2)))
C
C F_1
         CPART45=VD*SCHI(NCH2)*(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)*CDEN
C
         CPART46=-PI*COEF*CPART45*PN0*1.D20*DELTH
     &        /(ABS(-DKPRY+MM*AA*PNEAR1*(1+TCSM2(NTH))*CCHI(NCH2)))
C
C F_0+F_1
         CPART47=CPART44+CPART46
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VD*VPARA
C
C
         CINTG711 = CINTG711 + CPART47*PI1*CCHI(NCH1)*CCHI(NCH2)
         CINTG712 = CINTG712 + CPART47*PI1*CCHI(NCH1)*SCHI(NCH2)
         CINTG713 = CINTG713 + CPART47*PI2*CCHI(NCH1)
         CINTG721 = CINTG721 + CPART47*PI1*SCHI(NCH1)*CCHI(NCH2)
         CINTG722 = CINTG722 + CPART47*PI1*SCHI(NCH1)*SCHI(NCH2)
         CINTG723 = CINTG723 + CPART47*PI2*SCHI(NCH1)
         CINTG731 = CINTG731 + CPART47*PI2           *CCHI(NCH2)           
         CINTG732 = CINTG732 + CPART47*PI2           *SCHI(NCH2)
         CINTG733 = CINTG733 + CPART47*PI3
      ENDDO
C
C
      CDTNSI(1,1)=CI*(-CI*PI*(CINTG611+CINTG711))/(CW*EPS0)
      CDTNSI(1,2)=CI*(-CI*PI*(CINTG612+CINTG712))/(CW*EPS0)
      CDTNSI(1,3)=CI*(-CI*PI*(CINTG613+CINTG713))/(CW*EPS0)
      CDTNSI(2,1)=CI*(-CI*PI*(CINTG621+CINTG721))/(CW*EPS0)
      CDTNSI(2,2)=CI*(-CI*PI*(CINTG622+CINTG722))/(CW*EPS0)
      CDTNSI(2,3)=CI*(-CI*PI*(CINTG623+CINTG723))/(CW*EPS0)
      CDTNSI(3,1)=CI*(-CI*PI*(CINTG631+CINTG731))/(CW*EPS0)
      CDTNSI(3,2)=CI*(-CI*PI*(CINTG632+CINTG732))/(CW*EPS0)
      CDTNSI(3,3)=CI*(-CI*PI*(CINTG633+CINTG733))/(CW*EPS0)
C   
C 
      RETURN
      END                     
