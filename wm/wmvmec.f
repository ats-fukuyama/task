C     $Id$
C
C    ****** RADIAL MESH AND METRIC TENSOR ******
C
C
      SUBROUTINE WMXRZV(IERR)
C
      INCLUDE 'wmcomm.inc'
C     
      IERR=0
C
C      IF(ID.EQ.4) THEN
         CALL WMHFRD
C      ELSE
C         CALL WMHGRD(ID,IERR)
C         IF(IERR.NE.0) RETURN
C      ENDIF
      CALL WMHCRZ
      CALL WMHMTR
      CALL WMHCBB
C
      RETURN
      END
C
C    ****** READ WOUT FILE ******
C
      SUBROUTINE WMHFRD
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      CHARACTER FNAME*80
C
C 1002    WRITE(6,*) 'READ FILE NAME ?'
C         READ(5,'(A80)',ERR=1002,END=9000) FNAME
          FNAME='wout'
C
      OPEN(8,FILE=FNAME,FORM='UNFORMATTED',STATUS='OLD')
C
      READ(8) (R1DATA(I),I=1,3),MNMAX,NTOR0,MN0,N1,N2
      VOLI=DBLE(R1DATA(1))
      RGAM=DBLE(R1DATA(2))
      RC  =DBLE(R1DATA(3))
C
      NSRMAX=31
      WRITE(6,*) 'MNMAX,NTOR0,MN0=',MNMAX,NTOR0,MN0
      WRITE(6,*) 'VOLI,RGAM,RC   =',VOLI,RGAM,RC
C
      DO MN=1,NSRMAX*MNMAX
         READ(8) (R1DATA(I),I=1,12)
         XM(MN)   =R1DATA( 1)
         XN(MN)   =R1DATA( 2)
         RMNC(MN) =R1DATA( 3)
         ZMNS(MN) =R1DATA( 4)
C         WRITE(6,*) MN,XM(MN),XN(MN),RMNC(MN),ZMNS(MN)
         RLMNS(MN)=R1DATA( 5)
         BMOD(MN) =R1DATA( 6)
         RGMOD(MN)=R1DATA( 7)
C         BSUBU    =R1DATA( 8)
C         BSUBV    =R1DATA( 9)
C         BSUBS    =R1DATA(10)
         BSU(MN)  =R1DATA(11)
         BSV(MN)  =R1DATA(12)
      ENDDO
C
      READ(8) ((R2DATA(I,NSR),I=1,11),NSR=2,NSRMAX)
      DO NSR=2,NSRMAX
         RIOTAS(NSR) =R2DATA( 1,NSR)
         RMASS(NSR)  =R2DATA( 2,NSR)
         PRES(NSR)   =R2DATA( 3,NSR)
         PHIPS(NSR)  =R2DATA( 4,NSR)
         BPCO(NSR)   =R2DATA( 5,NSR)
         BACO(NSR)   =R2DATA( 6,NSR)
         PHI(NSR)    =R2DATA( 7,NSR)
         VP(NSR)     =R2DATA( 8,NSR)
         RJTHETA(NSR)=R2DATA( 9,NSR)
         RJZETA(NSR) =R2DATA(10,NSR)
         SPECW(NSR)  =R2DATA(11,NSR)
      ENDDO
C
      READ(8) ((R3DATA(I,J),I=1,2),J=1,100)
      DO J=1,100
         FSQT(I)=R3DATA(1,J)
         WDOT(I)=R3DATA(2,J)
      ENDDO
C
      CLOSE(8)
C
C     ----- Setup normalized radius for file data -----
C
      PSIA =-PHI(NSRMAX)
      DO NSR=1,NSRMAX
         S=-PHI(NSR)/PSIA
         IF(S.LT.0.D0) THEN
            XS(NSR)=0.D0
         ELSE
            XS(NSR)=SQRT(S)
         ENDIF
         IF(NSR.EQ.1) THEN
            XSH(NSR)=0.D0
         ELSE
            XSH(NSR)=0.5D0*(XS(NSR)+XS(NSR-1))
         ENDIF
C         WRITE(6,'(I5,2X,1P3E12.5)') NSR,PHI(NSR),XS(NSR),PRES(NSR)
      ENDDO
C
      RHOB=RB/RA
      DRHO=RHOB/NRMAX
      DO NR=1,NRMAX+1
         XRHO(NR)=DRHO*(NR-1)
         XR(NR)  =RB*XRHO(NR)
      ENDDO
C
      RETURN
      END
C
C    ****** READ WOUT FILE ******
C
      SUBROUTINE WMHGRD(ID,IERR)
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      CHARACTER FNAME*80
C
      IERR=0
C
C 1002 WRITE(6,*) 'READ FILE NAME ?'
C      READ(5,'(A80)',ERR=1002,END=9000) FNAME
      IF(ID.EQ.5) THEN
         FNAME='woutS'
      ELSEIF(ID.EQ.6) THEN
         FNAME='woutT'
      ELSE
         WRITE(6,*) 'XX WMHGRD: UNKNOWN ID :',ID
         IERR=1
         RETURN
      ENDIF
C
      OPEN(8,FILE=FNAME,FORM='UNFORMATTED',STATUS='OLD')
C
      READ(8) VOLI,RGAM,RC,D2,D3,N1,
     &        MNMAX,NSRMAX,NTHET1,NZETA1,
     &        N2,N3,N4,N5,NTOR0,MN0
C
C      WRITE(6,*) 'NSRMAX=',NSRMAX
C      WRITE(6,*) 'MNMAX =',MNMAX
C      WRITE(6,*) 'NTHET1=',NTHET1
C      WRITE(6,*) 'NZETA1=',NZETA1
C
      NTHET2=1+NTHET1/2
      NZNT=NZETA1*NTHET2
C
C
      MN=0
      DO NSR=1,NSRMAX
         DO I=1,MNMAX
            MN=MN+1
            READ(8) XM(MN),XN(MN),
     &              RMNC(MN),ZMNS(MN),RLMNS(MN),
     &              BMOD(MN),RGMOD(MN),
     &              D1,D2,D3,
     &              BSU(MN),BSV(MN)
C            IF(MN.LE.100) WRITE(6,'(A,I5,1P5E12.4)')
C     &              'MN,M,N,R,Z,B=',
C     &               MN,XM(MN),XN(MN),RMNC(MN),ZMNS(MN),BMOD(MN)
         ENDDO
C
         READ(8) (D1,D2,D3,D4,D5,D6,D7,K=1,NZNT)
      ENDDO
C
      READ(8) (RIOTAS(NSR),PRES(NSR),VP(NSR),PHIPS(NSR),
     &         BPCO(NSR),BACO(NSR),
     &         PHI(NSR),RCHI(NSR),RJTHETA(NSR),RJZETA(NSR),
     &         SPECW(NSR),NSR=2,NSRMAX)
C
      IF(ID.EQ.6) THEN
         READ(8) (AMVM(K),K=0,5)
      ENDIF
      READ(8) (FSQT(I),WDOT(I),I=1,100)

C
      CLOSE(8)
C
C     ----- Setup normalized radius for file data -----
C
      PSIA =-PHI(NSRMAX)
      DO NSR=1,NSRMAX
         S=-PHI(NSR)/PSIA
         IF(S.LT.0.D0) THEN
            XS(NSR)=0.D0
         ELSE
            XS(NSR)=SQRT(S)
         ENDIF
         IF(NSR.EQ.1) THEN
            XSH(NSR)=0.D0
         ELSE
            XSH(NSR)=0.5D0*(XS(NSR)+XS(NSR-1))
         ENDIF
C         WRITE(6,'(I5,2X,1P3E12.5)') NSR,PHI(NSR),XS(NSR),PRES(NSR)
      ENDDO
C
      RHOB=RB/RA
      DRHO=RHOB/NRMAX
      DO NR=1,NRMAX+1
         XRHO(NR)=DRHO*(NR-1)
         XR(NR)  =RB*XRHO(NR)
      ENDDO
C
      RETURN
      END
C
C    ****** CALCULATE R AND Z ******
C
      SUBROUTINE WMHCRZ
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C      ***** SPLINE PSIPS *****
C
      CALL SPL1D(XS,PHI,FX1,U1,NSRMAX,0,IERR)
C
      DO NR=1,NRMAX+1
         CALL SPL1DF(XRHO(NR),PSIPS(NR),XS,U1,NSRMAX,IERR)
      ENDDO
C     
C      ***** SPLINE RMNC(S),ZMNS(S) DRMNC(S) DZMNS(S) *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            YRBS(NSR)=RMNC(MN+(NSR-1)*MNMAX)
            YZBS(NSR)=ZMNS(MN+(NSR-1)*MNMAX)
         ENDDO
C         IF(MYRANK.EQ.0) THEN
C            WRITE(6,'( I5,1P6E12.4)') MN,XM(MN),XN(MN),
C     &                                YRBS(1),YRBS(2),YRBS(3),YRBS(4)
C            WRITE(6,'(29X,1P4E12.4)') YZBS(1),YZBS(2),YZBS(3),YZBS(4)
C         ENDIF
         CALL SPL1D(XS,YRBS,FX1,U1,NSRMAX,0,IERR)
         CALL SPL1D(XS,YZBS,FX2,U2,NSRMAX,0,IERR)
C
         DO NR=1,NRMAX+1
            CALL SPL1DD(XRHO(NR),SRMNC(MN,NR),DRMNC(MN,NR),
     &                  XS,U1,NSRMAX,IERR)
C            WRITE(6,'(3I5,1P2E12.4)') MN,NR,IERR,XRHO(NR),XS(NSRMAX)
C            IF(IERR.EQ.2) THEN
            IF(XRHO(NR)-XS(NSRMAX).GT.0.D0) THEN
               DRMNC(MN,NR)=(YRBS(NSRMAX)-YRBS(NSRMAX-1))
     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
               SRMNC(MN,NR)=YRBS(NSRMAX)
     &                     +DRMNC(MN,NR)*(XRHO(NR)-XS(NSRMAX))
            ENDIF
            CALL SPL1DD(XRHO(NR),SZMNS(MN,NR),DZMNS(MN,NR),
     &                  XS,U2,NSRMAX,IERR)
C            IF(IERR.EQ.2) THEN
            IF(XRHO(NR)-XS(NSRMAX).GT.0.D0) THEN
               DZMNS(MN,NR)=(YZBS(NSRMAX)-YZBS(NSRMAX-1))
     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
               SZMNS(MN,NR)=YZBS(NSRMAX)
     &                     +DZMNS(MN,NR)*(XRHO(NR)-XS(NSRMAX))
            ENDIF
         ENDDO
      ENDDO
C
C     ***** CALCULATE R,Z *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            RMNCC(MN,NSR)=RMNC(MN+(NSR-1)*MNMAX)
            ZMNSS(MN,NSR)=ZMNS(MN+(NSR-1)*MNMAX)
         ENDDO
      ENDDO
      RETURN
      END
C
C    ****** CALCULTE METRIC ******
C
      SUBROUTINE WMHMTR
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C      ***** CULCULATE METRIC TENSOR AND JACOBIAN*****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NPHMAX)
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
         RPSS=0.D0
         ZPSS=0.D0
         DRS=0.D0
         DZS=0.D0
         DRTH=0.D0
         DZTH=0.D0
         DRPH=0.D0
         DZPH=0.D0
         DPHPH=0.D0
C
         TH=DTH*(NTH-1)
         PH=DPH*(NPH-1)
         DO MN=1,MNMAX
            RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
            RCOS=COS(XM(MN)*TH-XN(MN)*PH)
            RPSS =RPSS +       SRMNC(MN,NR)*RCOS
            ZPSS =ZPSS +       SZMNS(MN,NR)*RSIN
            DRS  =DRS  +       DRMNC(MN,NR)*RCOS
            DZS  =DZS  +       DZMNS(MN,NR)*RSIN
            DRTH =DRTH -XM(MN)*SRMNC(MN,NR)*RSIN
            DZTH =DZTH +XM(MN)*SZMNS(MN,NR)*RCOS
            DRPH =DRPH +XN(MN)*SRMNC(MN,NR)*RSIN
            DZPH =DZPH -XN(MN)*SZMNS(MN,NR)*RCOS
            DPHPH=DPHPH+       SRMNC(MN,NR)*RCOS
         ENDDO
C
C  *****  DRS,DZS,DPHS,DRTH,DZTH,DPHTH,DRPH,DZPH,DPHPH GRAPH *****
C            
         RPST(  NTH,NPH,NR)=RPSS
         ZPST(  NTH,NPH,NR)=ZPSS
C
         IF(NR.NE.1) THEN
            DRS=DRS/(2.D0*PSIA*XRHO(NR))
            DZS=DZS/(2.D0*PSIA*XRHO(NR))
            RG11(NTH,NPH,NR)=(DRS *DRS +DZS *DZS )*XRHO(NR)**2
            RG12(NTH,NPH,NR)=(DRS *DRTH+DZS *DZTH)
            RG13(NTH,NPH,NR)=(DRS *DRPH+DZS *DZPH)*XRHO(NR)
            RG22(NTH,NPH,NR)=(DRTH*DRTH+DZTH*DZTH)/XRHO(NR)**2
            RG23(NTH,NPH,NR)=(DRTH*DRPH+DZTH*DZPH)/XRHO(NR)
            RG33(NTH,NPH,NR)= DRPH*DRPH+DZPH*DZPH+DPHPH*DPHPH
            RJ  (NTH,NPH,NR)=DPHPH*(DRS*DZTH-DZS*DRTH)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         RG11(NTH,NPH,1)=RG11(NTH,NPH,2)
         RG12(NTH,NPH,1)=RG12(NTH,NPH,2)
         RG13(NTH,NPH,1)=RG13(NTH,NPH,2)
         RG22(NTH,NPH,1)=RG22(NTH,NPH,2)
         RG23(NTH,NPH,1)=RG23(NTH,NPH,2)
         RG33(NTH,NPH,1)=RG33(NTH,NPH,2)
         RJ(NTH,NPH,1)  =RJ(NTH,NPH,2)
      ENDDO
      ENDDO
C
      RETURN
      END
C
C    ***** SPLINE POLOIDAL AND TOROIDAL MAGNETIC FIELD *****
C
      SUBROUTINE WMHCBB
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            BSUS(NSR)=BSU(MN+(NSR-1)*MNMAX)
            BSVS(NSR)=BSV(MN+(NSR-1)*MNMAX)
         ENDDO
C         IF(MYRANK.EQ.0) THEN
C            WRITE(6,'( I5,1PE10.2,1P5E12.4)') MN,XM(MN),
C     &                      BSUS(1),BSUS(2),BSUS(3),BSUS(4),BSUS(5)
C            WRITE(6,'( 5X,1PE10.2,1P5E12.4)') XN(MN),
C     &                      BSVS(1),BSVS(2),BSVS(3),BSVS(4),BSVS(5)
C         ENDIF
         IF(XM(MN).EQ.0) THEN
            BSUS(1)=(4*BSUS(2)-BSUS(3))/3.D0
            BSVS(1)=(4*BSVS(2)-BSVS(3))/3.D0
C            BSUS(1)=(2*BSUS(2)-BSUS(3))/1.D0
C            BSVS(1)=(2*BSVS(2)-BSVS(3))/1.D0
            FX1(1)=0.D0
            FX2(1)=0.D0
            CALL SPL1D(XS,BSUS,FX1,U1,NSRMAX,1,IERR)
            CALL SPL1D(XS,BSVS,FX2,U2,NSRMAX,1,IERR)
         ELSEIF(XM(MN).EQ.1) THEN
C            BSUS(1)=(SQRT(2.D0)*BSUS(2)-BSUS(3))/(SQRT(2.D0)-1.D0)
C            BSVS(1)=(SQRT(2.D0)*BSVS(2)-BSVS(3))/(SQRT(2.D0)-1.D0)
            BSUS(1)=(2*BSUS(2)-BSUS(3))/1.D0
            BSVS(1)=(2*BSVS(2)-BSVS(3))/1.D0
C            BSUS(1)= 0.D0
C            BSVS(1)= 0.D0
            CALL SPL1D(XS,BSUS,FX1,U1,NSRMAX,0,IERR)
            CALL SPL1D(XS,BSVS,FX2,U2,NSRMAX,0,IERR)
         ELSE
            BSUS(1)= 0.D0
            BSVS(1)= 0.D0
C            FX1(1)=0.D0
C            FX2(1)=0.D0
            CALL SPL1D(XS,BSUS,FX1,U1,NSRMAX,0,IERR)
            CALL SPL1D(XS,BSVS,FX2,U2,NSRMAX,0,IERR)
         ENDIF
C
         DO NR=1,NRMAX+1
            CALL SPL1DF(XRHO(NR),BSTHL,XS,U1,NSRMAX,IERR)
C            IF(IERR.EQ.2) THEN
            IF(XRHO(NR)-XS(NSRMAX).GT.0.D0) THEN
               BSTHL=BSUS(NSRMAX)+(BSUS(NSRMAX)-BSUS(NSRMAX-1))
     &                           /(XS(NSRMAX)-XS(NSRMAX-1))
     &                           *(XRHO(NR)-XS(NSRMAX))
            ENDIF
            BSTH(MN,NR)=BSTHL
            IF(XRHO(NR)-XS(NSRMAX).GT.0.D0) THEN
               BSPHL=BSVS(NSRMAX)+(BSVS(NSRMAX)-BSVS(NSRMAX-1))
     &                           /(XS(NSRMAX)-XS(NSRMAX-1))
     &                           *(XRHO(NR)-XS(NSRMAX))
            ELSE
               CALL SPL1DF(XRHO(NR),BSPHL,XS,U2,NSRMAX,IERR)
            ENDIF
            BSPH(MN,NR)=BSPHL
         ENDDO
      ENDDO
C
C      ***** CULCULATE MAGNETIC FIELD *****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NPHMAX)
C      DO NR=2,NRMAX+1
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         TH=DTH*(NTH-1)
         PH=DPH*(NPH-1)
C
         SBTH=0.D0
         SBPH=0.D0
         DO MN=1,MNMAX
            RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
            RCOS=COS(XM(MN)*TH-XN(MN)*PH)
            SBTH=SBTH+BSTH(MN,NR)*RCOS
            SBPH=SBPH+BSPH(MN,NR)*RCOS
         ENDDO
         BFLD(2,NTH,NPH,NR)=SBTH
         BFLD(3,NTH,NPH,NR)=SBPH
      ENDDO
      ENDDO
      ENDDO
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C         BFLD(2,NTH,NPH,1)=(4*BFLD(2,NTH,NPH,2)-BFLD(2,NTH,NPH,3))/3
C         BFLD(3,NTH,NPH,1)=(4*BFLD(3,NTH,NPH,2)-BFLD(3,NTH,NPH,3))/3
      ENDDO
      ENDDO
C     
C    ***************************************
C
      P0=0.D0
      DO NS=1,NSMAX
         P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
      ENDDO
      P0=P0*1.D20*AEE*1.D3/1.D6
C
      DO NR=1,NRMAX+1
         RHOL=XRHO(NR)
         IF(RHOL.LE.1.D0) THEN
            IF(PN(1).LE.0.D0) THEN
               FACTN=0.D0
            ELSE
               FEDGE=PNS(1)/PN(1)
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1)**PROFN2+FEDGE
            ENDIF
            PT=(PTPR(1)+2*PTPP(1))/3.D0
            IF(PT.LE.0.D0) THEN
               FACTT=0.D0
            ELSE
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1)**PROFT2+FEDGE
            ENDIF
            PPS(NR)=P0*FACTN*FACTT
         ELSE
            PPS(NR)=0.D0
         ENDIF
      ENDDO
C
      NSUMAX=31
      DTHU=2.D0*PI/(NSUMAX-1)
      DO NSU=1,NSUMAX
         DO NPH=1,NPHMAX
            RSU(NSU,NPH)=0.D0
            ZSU(NSU,NPH)=0.D0
            TH=DTHU*(NSU-1)
            PH=DPH*(NPH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSU(NSU,NPH)=RSU(NSU,NPH)+RMNCC(MN,NSRMAX)*RCOS
               ZSU(NSU,NPH)=ZSU(NSU,NPH)+ZMNSS(MN,NSRMAX)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
      NSWMAX=31
      DTHW=2.D0*PI/(NSWMAX-1)
      DO NSW=1,NSWMAX
         DO NPH=1,NPHMAX
            RSW(NSW,NPH)=0.D0
            ZSW(NSW,NPH)=0.D0
            TH=DTHW*(NSW-1)
            PH=DPH*(NPH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSW(NSW,NPH)=RSW(NSW,NPH)+SRMNC(MN,NRMAX+1)*RCOS
               ZSW(NSW,NPH)=ZSW(NSW,NPH)+SZMNS(MN,NRMAX+1)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
C     ***** SPLINE IOTAS AND CULCULATE QPS *****
C
      RIOTAS(1)=2.D0*RIOTAS(2)-RIOTAS(3)
C
      CALL SPL1D(XS,RIOTAS,FX3,U3,NSRMAX,0,IERR)
C
      DO NR=1,NRMAX+1
         CALL SPL1DF(XRHO(NR),RIOTASL,XS,U3,NSRMAX,IERR)
         IF(IERR.EQ.2) THEN
            RIOTASL=RIOTAS(NSRMAX)+(RIOTAS(NSRMAX)-RIOTAS(NSRMAX-1))
     &                      /(XS(NSRMAX)-XS(NSRMAX-1))
     &                      *(XRHO(NR)-XS(NSRMAX))
         ENDIF
         QPS(NR)=2.D0*PI/RIOTASL
      ENDDO                     
C
C     ***** COMPUTE R,Z MAGNETIC AXES *****
C
      RGMIN=RSW(1,1)
      RGMAX=RSW(1,1)
      ZGMIN=ZSW(1,1)
      ZGMAX=ZSW(1,1)
      DO NPH=1,NPHMAX
         DO NSW=1,NSWMAX
            RGMIN=MIN(RGMIN,RSW(NSW,NPH))
            RGMAX=MAX(RGMAX,RSW(NSW,NPH))
            ZGMIN=MIN(ZGMIN,ZSW(NSW,NPH))
            ZGMAX=MAX(ZGMAX,ZSW(NSW,NPH))
         ENDDO
      ENDDO
      RR=0.5D0*(RGMIN+RGMAX)
C
      RETURN
      END
C
C   ***** DRAW R,Z GRAPH *****
C
      SUBROUTINE WMGS1D(XG,YG,NXGMAX,KTITL,NGD)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION XG(NXGMAX),YG(NXGMAX)
      DIMENSION GX(NRM),GY(NRM)
      CHARACTER KTITL*6
C
      DIMENSION GP(4,4)
      DATA GP/ 3.0, 10.8,  9.5, 16.5,
     &         3.0, 10.8,  1.0,  8.0,
     &        13.8, 21.6,  9.5, 16.5,
     &        13.8, 21.6,  1.0,  8.0/
C
      DO NXG=1,NXGMAX
         GX(NXG)=GUCLIP(XG(NXG))
         GY(NXG)=GUCLIP(YG(NXG))
      ENDDO
      GXMIN=0.0
      GXMAX=GUCLIP(RB/RA)
C
      CALL WMGSUB(NXGMAX,GX,GXMIN,GXMAX,GY,1,GP(1,NGD),KTITL)
C
      RETURN
      END
