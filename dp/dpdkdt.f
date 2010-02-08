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
      CWP=PN0**1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
C
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
     &                      /(PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH,NP) = (FP(NTH,NP+1,NR+1)-FP(NTH,NP,NR+1))
     &                      /(PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH,NP) = (FP(NTH,NP+1,NR-1)-FP(NTH,NP,NR-1))
     &                      /(PTH0*DELP*FP(NTH,NP,NR-1))
         DO NP=2,NPMAX-2
            DPFP(NTH,NP)   = (FP(NTH,NP+1,NR  )-FP(NTH,NP-1,NR  ))
     &                      /(2*PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH,NP) = (FP(NTH,NP+1,NR+1)-FP(NTH,NP-1,NR+1))
     &                      /(2*PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH,NP) = (FP(NTH,NP+1,NR-1)-FP(NTH,NP-1,NR-1))
     &                      /(2*PTH0*DELP*FP(NTH,NP,NR-1))
         ENDDO
         NP=NPMAX-1
            DPFP(NTH, NP)   = (FP(NTH,NP,NR  )-FP(NTH,NP-1,NR  ))
     &                       /(PTH0*DELP*FP(NTH,NP,NR))
            DPFPRP(NTH, NP) = (FP(NTH,NP,NR+1)-FP(NTH,NP-1,NR+1))
     &                       /(PTH0*DELP*FP(NTH,NP,NR+1))
            DPFPRM(NTH, NP) = (FP(NTH,NP,NR-1)-FP(NTH,NP-1,NR-1))
     &                       /(PTH0*DELP*FP(NTH,NP,NR-1))
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
         PV(NP)=PTH0*PG(NP)/AM
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
      CINTG111a = (0.D0,0.D0)
      CINTG112a = (0.D0,0.D0)
      CINTG113a = (0.D0,0.D0)
      CINTG121a = (0.D0,0.D0)
      CINTG122a = (0.D0,0.D0)
      CINTG123a = (0.D0,0.D0)
      CINTG131a = (0.D0,0.D0)
      CINTG132a = (0.D0,0.D0)
      CINTG133a = (0.D0,0.D0)
C
      CINTG111b = (0.D0,0.D0)
      CINTG112b = (0.D0,0.D0)
      CINTG113b = (0.D0,0.D0)
      CINTG121b = (0.D0,0.D0)
      CINTG122b = (0.D0,0.D0)
      CINTG123b = (0.D0,0.D0)
      CINTG131b = (0.D0,0.D0)
      CINTG132b = (0.D0,0.D0)
      CINTG133b = (0.D0,0.D0)
C
      CINTG111c = (0.D0,0.D0)
      CINTG112c = (0.D0,0.D0)
      CINTG113c = (0.D0,0.D0)
      CINTG121c = (0.D0,0.D0)
      CINTG122c = (0.D0,0.D0)
      CINTG123c = (0.D0,0.D0)
      CINTG131c = (0.D0,0.D0)
      CINTG132c = (0.D0,0.D0)
      CINTG133c = (0.D0,0.D0)
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
         PV(NP)=PTH0*PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
C
         VD=PV(NP)**2*(1+TCSM2(NTH))
         VPARA=PV(NP)*TCSM(NTH)
C         COEF=AE**2*FP(NTH,NP,NR)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)*0.5D0*AA/RSL
C         CDEN=1.D0/CDENX
C
         CDEN  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*PTH0*DELP)**2
     &                             +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART11=PV(NP)**2*TSNM(NTH)*FP(NTH,NP,NR)
     &           *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)
C         CPART12=(0.D0,0.D0)
C         CPART13=(0.D0,0.D0)
C
         CPART12=FP(NTH,NP,NR)*TSNM(NTH)*PV(NP)**6*(1+TCSM2(NTH))**2
     &          *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/CDENX**2
C
         CPART13=FP(NTH,NP,NR)*TSNM(NTH)*PV(NP)**4*(1+TCSM2(NTH))
     &           *(DRPFP(NTH,NP)-DRCWAST(NTH,NP)/CW
     &           +CWAST(NTH,NP)/(CW*RSL))/CDENX
C      
C         CPART14= CI*2*PI*PM(NP)**2*TSNM(NTH)*DELTH*DELP
C     &            *COEF*(CPART11+CPART12+CPART13)*PN0*1.D20
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VPARA*VPARA
C
         CINTG111a = CINTG111a + CDEN*CPART11*PI1*DELTH*PTH0*DELP
         CINTG112a = CINTG111a 
         CINTG113a = CINTG113a + CDEN*CPART11*PI2*DELTH*PTH0*DELP
         CINTG121a = CINTG111a 
         CINTG122a = CINTG111a 
         CINTG123a = CINTG113a 
         CINTG131a = CINTG113a 
         CINTG132a = CINTG113a 
         CINTG133a = CINTG133a + CDEN*CPART11*PI3*DELTH*PTH0*DELP
C
         CINTG111b = CINTG111b + CDEN*CPART12*PI1*DELTH*PTH0*DELP
         CINTG112b = CINTG111b 
         CINTG113b = CINTG113b + CDEN*CPART12*PI2*DELTH*PTH0*DELP
         CINTG121b = CINTG111b 
         CINTG122b = CINTG111b 
         CINTG123b = CINTG113b 
         CINTG131b = CINTG113b 
         CINTG132b = CINTG113b 
         CINTG133b = CINTG133b + CDEN*CPART12*PI3*DELTH*PTH0*DELP
C
         CINTG111c = CINTG111c + CDEN*CPART13*PI1*DELTH*PTH0*DELP
         CINTG112c = CINTG111c 
         CINTG113c = CINTG113c + CDEN*CPART13*PI2*DELTH*PTH0*DELP
         CINTG121c = CINTG111c 
         CINTG122c = CINTG111c 
         CINTG123c = CINTG113c 
         CINTG131c = CINTG113c 
         CINTG132c = CINTG113c 
         CINTG133c = CINTG133c + CDEN*CPART13*PI3*DELTH*PTH0*DELP
C
      ENDDO
      ENDDO
C
      COEA = CI*2*PI*AE**2
      COEB = -AE**2*PI*MM*AA**2*0.5D0*CCHI(NCH2)*SCHI(NCH2)
      COEC = -AE**2*PI*AA*SCHI(NCH2)
C
      FACT=2.D0*PI*DBLE(CWP)
C
      CINTG111 = (COEA*CINTG111a+COEB*CINTG111b+COEC*CINTG111c)
     &     *CCHI(NCH1)*CCHI(NCH2)*FACT
      CINTG112 = (COEA*CINTG112a+COEB*CINTG112b+COEC*CINTG112c)
     &     *CCHI(NCH1)*SCHI(NCH2)*FACT
      CINTG113 = (COEA*CINTG113a+COEB*CINTG113b+COEC*CINTG113c)
     &     *CCHI(NCH1)*FACT
      CINTG121 = (COEA*CINTG121a+COEB*CINTG121b+COEC*CINTG121c)
     &     *SCHI(NCH1)*CCHI(NCH2)*FACT
      CINTG122 = (COEA*CINTG122a+COEB*CINTG122b+COEC*CINTG122c)
     &     *SCHI(NCH1)*SCHI(NCH2)*FACT
      CINTG123 = (COEA*CINTG123a+COEB*CINTG123b+COEC*CINTG123c)
     &     *SCHI(NCH1)*FACT
      CINTG131 = (COEA*CINTG131a+COEB*CINTG131b+COEC*CINTG131c)
     &     *CCHI(NCH2)*FACT
      CINTG132 = (COEA*CINTG132a+COEB*CINTG132b+COEC*CINTG132c)           
     &     *SCHI(NCH2)*FACT
      CINTG133 = (COEA*CINTG133a+COEB*CINTG133b+COEC*CINTG133c)*FACT
C
C
C*************SUM2********************
C
      CINTG211a = (0.D0,0.D0)
      CINTG212a = (0.D0,0.D0)
      CINTG213a = (0.D0,0.D0)
      CINTG221a = (0.D0,0.D0)
      CINTG222a = (0.D0,0.D0)
      CINTG223a = (0.D0,0.D0)
      CINTG231a = (0.D0,0.D0)
      CINTG232a = (0.D0,0.D0)
      CINTG233a = (0.D0,0.D0)
C

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
         PV(NP)=PTH0*PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         VD=PV(NP)**2*(1+TCSM2(NTH))
         VPARA=PV(NP)*TCSM(NTH)
C         COEF=AE**2*FP(NTH,NP,NR)
C
         CDENX=CW-DKPRX+MM*VD*CCHI(NCH2)*0.5D0*AA/RSL
C         CDENX = CW-DKPRX+MM*AA*PV(NP)**2*(1+TCSM2(NTH)*CCHI(NCH2))
C     &        /(2*RSL)
         CDEN  = CDENX/(CDENX**2+DELPL*(DGP(NP,NTH)*PTH0*DELP)**2
     &                          +DELPL*(DGT(NP,NTH)*DELTH)**2)
C
         CPART2= PV(NP)**4*TSNM(NTH)*(1+TCSM2(NTH))
     &        *(DPFP(NTH,NP)-CWAST(NTH,NP)/CW)/CDENX
C
         PI1=VD*VD
         PI2=VD*VPARA
         PI3=VPARA*VPARA
C
         CINTG211a = CINTG211a + CDEN*CPART2*PI1*DELTH*PTH0*DELP
         CINTG212a = CINTG211a
         CINTG213a = CINTG213a + CDEN*CPART2*PI2*DELTH*PTH0*DELP
         CINTG221a = CINTG211a 
         CINTG222a = CINTG211a 
         CINTG223a = CINTG213a 
         CINTG231a = CINTG213a 
         CINTG232a = CINTG213a 
         CINTG233a = CINTG233a + CDEN*CPART2*PI3*DELTH*PTH0*DELP
      ENDDO
      ENDDO
C
      COED = -PI*AE**2*AA*SCHI(NCH2)
C
      FACT=2.D0*PI*DBLE(CWP)
C
      CINTG211 = COED*CINTG211a*CCHI(NCH1)*CCHI(NCH2)*FACT
      CINTG212 = COED*CINTG212a*CCHI(NCH1)*SCHI(NCH2)*FACT
      CINTG213 = COED*CINTG213a*CCHI(NCH1)           *FACT
      CINTG221 = COED*CINTG221a*SCHI(NCH1)*CCHI(NCH2)*FACT
      CINTG222 = COED*CINTG222a*SCHI(NCH1)*SCHI(NCH2)*FACT
      CINTG223 = COED*CINTG223a*SCHI(NCH1)           *FACT
      CINTG231 = COED*CINTG231a           *CCHI(NCH2)*FACT
      CINTG232 = COED*CINTG232a           *SCHI(NCH2)*FACT
      CINTG233 = COED*CINTG233a*FACT
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
         PV(NP)=PTH0*PG(NP)/AM
         TCSM2(NTH)=TCSM(NTH)*2
         EPART1= DRFP(NTH,NP)*PV(NP)**4*(1+TCSM2(NTH))*TSNM(NTH)
         EPART2= EPART1
         EPART3= DRFP(NTH,NP)*PV(NP)**3*TCSM(NTH)*TSNM(NTH)
C
      CINTG312a = CINTG312a+EPART2*DELTH*PTH0*DELP
      CINTG313a = CINTG313a+EPART2*DELTH*PTH0*DELP
      CINTG322a = CINTG322a+EPART2*DELTH*PTH0*DELP
      CINTG323a = CINTG323a+EPART2*DELTH*PTH0*DELP
      CINTG332a = CINTG332a+EPART2*DELTH*PTH0*DELP
      CINTG333a = CINTG333a+EPART3*DELTH*PTH0*DELP
C
      ENDDO
      ENDDO
C
      COEE1 = -CI*PI*AM/(CW*RR)*SCHI(NCH1)
      COEE2 = -CI*PI*AM/(CW*RR)*CCHI(NCH1)
      COEE3 = -CI*2*PI*AE*BABS/CW
C
      CINTG312 = COEE1*CINTG312a
      CINTG313 = COEE2*CINTG313a
      CINTG322 = COEE2*CINTG322a
      CINTG323 = COEE2*CINTG323a
      CINTG332 = COEE2*CINTG332a
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
         PV(NP)=PTH0*PG(NP)/AM
         FPART= FP(NTH,NP,NR)*PV(NP)**2*TSNM(NTH)
C
      CINTG412 = CINTG412+FPART*DELTH*PTH0*DELP
      CINTG421 = CINTG421+FPART*DELTH*PTH0*DELP
C
      ENDDO
      ENDDO
C
      COEE4 = 2*PI*AE/BABS
C
      CINTG412 = COEE4*CINTG412
      CINTG421 = COEE4*CINTG421
C
C*****************************************************
C
      CLDISP(1) = CINTG111+CINTG211
      CLDISP(2) = CINTG112+CINTG212+CINTG312+CINTG412
      CLDISP(3) = CINTG113+CINTG213+CINTG313
      CLDISP(4) = CINTG121+CINTG221+CINTG421
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
         CKPRY=CKPR*PV(NP)*TSNM(NTH)
         DKPRY=DBLE(CKPRY)
         TCSM2(NTH)=TCSM(NTH)**2
         COEF=AE**2*FP(NTH,NP,NR)
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
         PI11=(AA/2)**2*PNEAR1**4*(1+TCSM2(NTH))**2
         PI12=(AA/2)*PNEAR1**4*(1+TCSM2(NTH))
         PI13=PNEAR1**4
C
C
         CINTG311 = CINTG311 + CPART36*PI11*CCHI(NCH1)*CCHI(NCH2)
         CINTG312 = CINTG312 + CPART36*PI11*CCHI(NCH1)*SCHI(NCH2)
         CINTG313 = CINTG313 + CPART36*PI12*CCHI(NCH1)*CCHI(NCH2)
         CINTG321 = CINTG321 + CPART36*PI11*SCHI(NCH1)*CCHI(NCH2)
         CINTG322 = CINTG322 + CPART36*PI11*SCHI(NCH1)*SCHI(NCH2)
         CINTG323 = CINTG323 + CPART36*PI12*SCHI(NCH1)*CCHI(NCH2)
         CINTG331 = CINTG331 + CPART36*PI12*CCHI(NCH1)*CCHI(NCH2)
         CINTG332 = CINTG332 + CPART36*PI12*CCHI(NCH1)*SCHI(NCH2)
         CINTG333 = CINTG333 + CPART36*PI13*CCHI(NCH1)*CCHI(NCH2)
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
         PV(NP)=PTH0*PG(NP)/AM
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
         CINTG411 = CINTG411 + CPART46*PI21*CCHI(NCH1)*CCHI(NCH2)
         CINTG412 = CINTG412 + CPART46*PI21*CCHI(NCH1)*SCHI(NCH2)
         CINTG413 = CINTG413 + CPART46*PI22*CCHI(NCH1)*CCHI(NCH2)
         CINTG421 = CINTG421 + CPART46*PI21*SCHI(NCH1)*CCHI(NCH2)
         CINTG422 = CINTG422 + CPART46*PI21*SCHI(NCH1)*SCHI(NCH2)
         CINTG423 = CINTG423 + CPART46*PI22*SCHI(NCH1)*CCHI(NCH2)
         CINTG431 = CINTG431 + CPART46*PI22*CCHI(NCH1)*CCHI(NCH2)
         CINTG432 = CINTG432 + CPART46*PI22*CCHI(NCH1)*SCHI(NCH2)
         CINTG433 = CINTG433 + CPART46*PI23*CCHI(NCH1)*CCHI(NCH2)
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
