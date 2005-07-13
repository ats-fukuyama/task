C     $Id$
C
C   ************************************** 
C   **       CALCULATE EQUILIBRIUM      **
C   ************************************** 
C
      SUBROUTINE EQCALC(IERR)
C      
      INCLUDE '../eq/eqcomc.inc'
C
      IERR=0
      CALL EQMESH
      CALL EQPSIN
      CALL EQDEFB
      CALL EQLOOP(IERR)
C         IF(IERR.NE.0) RETURN
      CALL EQTORZ
      CALL EQCALP
      RETURN
      END
C
C   ************************************************ 
C   **       Mesh Definition (sigma, theta)       **
C   ************************************************ 
C
      SUBROUTINE EQMESH
C      
      INCLUDE '../eq/eqcomc.inc'
C
      DSG=1.D0/NSGMAX
      DTG=2.D0*PI/NTGMAX 
C
      DO NSG=1,NSGMAX
         SIGM(NSG)=DSG*(NSG-0.5D0)
      ENDDO
      DO NSG=1,NSGMAX+1
         SIGG(NSG)=DSG*(NSG-1.D0)
      ENDDO
      DO NTG=1,NTGMAX
         THGM(NTG)=DTG*(NTG-0.5D0)
      ENDDO
      DO NTG=1,NTGMAX+1
         THGG(NTG)=DTG*(NTG-1.D0)
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **                Initial Psi                 **
C   ************************************************
C
      SUBROUTINE EQPSIN
C
      INCLUDE '../eq/eqcomc.inc'
C
      RAXIS=RR
      ZAXIS=0.0D0
C
C     --- assuming elliptic crosssection, flat current profile ---
C
      PSI0=-0.5D0*RMU0*RIP*1.D6*RR
      PSITS=PI*RKAP*RA**2*BB
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
          PSI(NTG,NSG)=PSI0*(1-SIGM(NSG)*SIGM(NSG))
          DELPSI(NTG,NSG)=0.D0
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************  
C   **          Boundary Definition               **
C   ************************************************
C
      SUBROUTINE EQDEFB
C      
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION DRHOM(NTGM),DRHOG(NTGMP)
      EXTERNAL EQFBND
C
C     ------ Define criteria for the Brent method ------
C
      EPSZ=1.D-8
C
C     ------ Define boundary radius RHOM/G ------
C
      DO NTG=1,NTGMAX
         ZBRF=TAN(THGM(NTG))
         THDASH=ZBRENT(EQFBND,THGM(NTG)-1.0D0,THGM(NTG)+1.0D0,EPSZ)
         RHOM(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                    +RKAP**2*SIN(THDASH)**2)
C
         ZBRF=TAN(THGG(NTG))
         THDASH=ZBRENT(EQFBND,THGG(NTG)-1.0D0,THGG(NTG)+1.0D0,EPSZ)
         RHOG(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                    +RKAP**2*SIN(THDASH)**2)
      ENDDO
      RHOG(NTGMAX+1)=RHOG(1)
C
C     ------ Define Delta theta on the boundary ------
C
      DRHOM(1)=(RHOM(2)-RHOM(NTGMAX))/(2*DTG)
      DO NTG=2,NTGMAX-1
         DRHOM(NTG)=(RHOM(NTG+1)-RHOM(NTG-1))/(2*DTG)
      ENDDO
      DRHOM(NTGMAX)=(RHOM(1)-RHOM(NTGMAX-1))/(2*DTG)
C
      DRHOG(1)=(RHOG(2)-RHOG(NTGMAX))/(2*DTG)
      DO NTG=2,NTGMAX-1
         DRHOG(NTG)=(RHOG(NTG+1)-RHOG(NTG-1))/(2*DTG)
      ENDDO
      DRHOG(NTGMAX)=(RHOG(1)-RHOG(NTGMAX-1))/(2*DTG)
      DRHOG(NTGMAX+1)=DRHOG(1)
C
C     ------ Calculate major radius on grid points ------
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         RMG(NSG,NTG)=RR+SIGM(NSG)*RHOG(NTG)*COS(THGG(NTG))
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         RGM(NSG,NTG)=RR+SIGG(NSG)*RHOM(NTG)*COS(THGM(NTG))       
      ENDDO
      ENDDO
C
C     ------ Calculate factors for coefficient matrix ------
C
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         AA(NSG,NTG)=SIGG(NSG)/RGM(NSG,NTG)
     &              +(DRHOM(NTG)*DRHOM(NTG)*SIGG(NSG))
     &              /(RGM(NSG,NTG)*RHOM(NTG)*RHOM(NTG))  
         AB(NSG,NTG)=-DRHOM(NTG)/(RGM(NSG,NTG)*RHOM(NTG))
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         AC(NSG,NTG)=-DRHOG(NTG)/(RMG(NSG,NTG)*RHOG(NTG))
         AD(NSG,NTG)=1/(RMG(NSG,NTG)*SIGM(NSG))
      ENDDO
      ENDDO
C
C      DO NSG=1,2
C         DO NTG=1,NTGMAX
C            WRITE(6,'(2I5,1P4E12.4)') NSG,NTG,AA(NSG,NTG),AB(NSG,NTG),
C     &                                       AC(NSG,NTG),AD(NSG,NTG)
C         ENDDO
C      ENDDO
C      PAUSE

      RETURN
      END
C
C   ************************************************
C   **               Iteration Loop               **
C   ************************************************
C
      SUBROUTINE EQLOOP(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      IERR=0
C
      DO NLOOP=1,NLOOP_MAX_EQ
         CALL EQBAND
         CALL EQRHSV(IERR)
         IF(IERR.NE.0) RETURN
         CALL EQSOLV
C
         SUM0=0.D0
         SUM1=0.D0
         DO NSG=1,NSGMAX
            DO NTG=1,NTGMAX
               SUM0=SUM0+PSI(NTG,NSG)**2
               SUM1=SUM1+DELPSI(NTG,NSG)**2
            ENDDO
         ENDDO
         SUM=SQRT(SUM1/SUM0)
C
         IF(NPRINT.GE.1) THEN
            WRITE(6,'(A,1P4E14.6)')
     &        'SUM,RAXIS,ZAXIS,PSI0=',SUM,RAXIS,ZAXIS,PSI0
         ENDIF
         IF(SUM.LT.EPSEQ) RETURN
      ENDDO
C
      WRITE(6,*) 'XX EQLOOP: NLOOP exceeds NLOOP_MAX_EQ'
      IERR=100
      RETURN
      END
C
C   ***********************************************
C   **            Matrix calculation             **
C   ***********************************************
C
      SUBROUTINE EQBAND
C
      INCLUDE '../eq/eqcomc.inc'
C
C     ------ Define matrix length MMAX and matrix half width NBND ------
C
      MMAX=NSGMAX*NTGMAX
      NBND=2*NTGMAX
C
C     ------ Initialize band matrix coefficients Q ------
C
      DO N=1,MMAX
      DO M=1,2*NBND-1
          Q(M,N)=0.D0
      ENDDO
      ENDDO
C
C     ------ Calculate band matrix coefficients Q ------
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         I=(NSG-1)*NTGMAX+NTG
         IF(NTG.EQ.1)THEN
            Q(NBND         -1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
            Q(NBND  +NTGMAX-1,I)= (AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND+2*NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
         ELSE
            Q(NBND  -NTGMAX-1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
            Q(NBND         -1,I)=(AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND  +NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
         ENDIF
         Q(NBND-NTGMAX,I)=  AA(NSG,NTG)/(DSG*DSG)
     &                   -(AC(NSG,NTG+1)-AC(NSG,NTG))/(4*DSG*DTG)
         Q(NBND       ,I)=-(AA(NSG+1,NTG)+AA(NSG,NTG))/(DSG*DSG)
     &                   -(AD(NSG,NTG+1)+AD(NSG,NTG))/(DTG*DTG)
         Q(NBND+NTGMAX,I)=  AA(NSG+1,NTG)/(DSG*DSG)
     &                   +(AC(NSG,NTG+1)-AC(NSG,NTG))/(4*DSG*DTG)
C
C     ------ Set periodic condition ------
C
         IF(NTG.EQ.NTGMAX)THEN	
            Q(NBND-2*NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))
     &                            /(4*DSG*DTG)
            Q(NBND-  NTGMAX+1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                            /(4*DSG*DTG)
     &                          +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND         +1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                            /(4*DSG*DTG)
         ELSE
            Q(NBND-NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))/(4*DSG*DTG)
            Q(NBND       +1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                          /(4*DSG*DTG)
     &                        +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND+NTGMAX+1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                          /(4*DSG*DTG)
         ENDIF
         IF(NSG.EQ.1) THEN
            WRITE(6,'(2I5,1P3E12.4)') NSG,NTG,
     &      Q(NBND         -1,I),Q(NBND+NTGMAX,I),Q(NBND+2*NTGMAX,I)
         ENDIF
         IF(NSG.EQ.2) THEN
            WRITE(6,'(2I5,1P3E12.4)') NSG,NTG,
     &      Q(NBND-NTGMAX,I),Q(NBND         -1,I),Q(NBND+NTGMAX,I)
         ENDIF
      ENDDO
      ENDDO
C
C     ------ Set radial boundary condition ------
C
      DO I=1,NTGMAX
      DO J=1,NTGMAX/2
         Q(NBND+J-I,I)=Q(NBND+J-I,I)
     &                +Q(NBND+J-I-NTGMAX/2,I)
         Q(NBND+J-I-NTGMAX/2,I)=0.D0
         Q(NBND+J-I+NTGMAX/2,I)=Q(NBND+J-I+NTGMAX/2,I)
     &                         +Q(NBND+J-I-NTGMAX,I)
         Q(NBND+J-I-NTGMAX,I)=0.D0
      ENDDO
      ENDDO
      DO I=1,NTGMAX
      DO J=1,NTGMAX
         Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   =Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   -Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)
         Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)=0.D0
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **                RHS vector                  **
C   ************************************************
C
      SUBROUTINE EQRHSV(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      EXTERNAL EQPSID
C
      IERR=0
C
C     ----- calculate PSIRZ(R,Z) from PSI(sigma, theta) -----
C
      CALL EQTORZ
C
C     ----- calculate spline coef for psi(R,Z) -----
C
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,URZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IER)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX EQRHSV: SPL2D for PSIRZ: IER=',IER
         IERR=101
         RETURN
      ENDIF
C
C     ----- calculate PSIITB -----
C
      CALL EQCNVA(RHOITB**2,PSIITB)
C
C     ----- calculate position of magnetic axis -----
C
      DELT=1.D-8
      EPS=1.D-4
      ILMAX=40
      LIST=0
      RINIT=RAXIS
      ZINIT=ZAXIS
      RSAVE=RAXIS
      ZSAVE=ZAXIS
      CALL NEWTN(EQPSID,RINIT,ZINIT,RAXIS,ZAXIS,
     &           DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5,1P2E12.4)')
     &        'XX EQRHSV: NEWTN ERROR: IER=',IER,RSAVE,ZSAVE
         RAXIS=RSAVE
         ZAXIS=ZSAVE
      ENDIF
      IF(RAXIS.LE.RR+RB.AND.
     &   RAXIS.GE.RR-RB.AND.
     &   ZAXIS.LE.RKAP*RB.AND.
     &   ZAXIS.GE.-RKAP*RB) THEN
         PSI0=PSIG(RAXIS,ZAXIS)
         SAXIS=PSIG(RAXIS,ZAXIS)
      ELSE
         WRITE(6,'(A)') 'XX EQRHSV: AXIS OUT OF PLASMA:'
         IERR=102
         RETURN
      ENDIF
C
C     ----- calculate R(sigma, theta) and Z(sigma,theta) -----
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         RMM(NTG,NSG)=RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG))
         ZMM(NTG,NSG)=   SIGM(NSG)*RHOM(NTG)*SIN(THGM(NTG))
      ENDDO
      ENDDO
C
C     ----- calculate right hand side vector -----
C
      IMDLEQF=MOD(MDLEQF,5)
      JMDLEQF=MDLEQF/5
C
C     ----- Given pressure and toroidal current profiles -----
C
      IF(IMDLEQF.EQ.0) THEN
         RRC=RAXIS
         FJP=0.D0
         FJT=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQJPSI(PSIN,HJPSID,HJPSI)
            CALL EQTPSI(PSIN,TPSI,DTPSI)
            CALL EQOPSI(PSIN,OMGPSI,DOMGPSI)
            PP(NTG,NSG)=PPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-HJPSI
            HJP2A=EXP(RMM(NTG,NSG)**2*OMGPSI**2*AMP/(2.D0*TPSI))
            HJP2B=EXP(RRC**2         *OMGPSI**2*AMP/(2.D0*TPSI))
            HJP2C=HJP2A-(RRC**2/RMM(NTG,NSG)**2)*HJP2B
            HJP2D=HJP2A-(RRC**4/RMM(NTG,NSG)**4)*HJP2B
            HJP2E=-0.5D0*2.D0*PI*PPSI*RMM(NTG,NSG)**3
            HJP2F=AMP*(2.D0*OMGPSI*DOMGPSI/TPSI
     &           -DTPSI*OMGPSI**2/TPSI**2)
            HJP2(NTG,NSG)=HJP2C*HJP1(NTG,NSG)+HJP2D*HJP2E*HJP2F
            HJT2(NTG,NSG)=(RRC/RMM(NTG,NSG))*HJT1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
C
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT=FJT+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         TJ=(RIP*1.D6-FJP)/FJT
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT2(NTG,NSG)
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQJPSI(PSIN,HJPSID,HJPSI)
            CALL EQTPSI(PSIN,TPSI,DTPSI)
            CALL EQOPSI(PSIN,OMGPSI,DOMGPSI)
            TT(NTG,NSG)=SQRT((2.D0*PI*BB*RR)**2
     &                 +2.D0*RMU0*RRC
     &                 *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                 *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
            RHO(NTG,NSG)=(PPSI*AMP/TPSI)
     &                  *EXP(RMM(NTG,NSG)**2*OMGPSI**2*AMP/(2.D0*TPSI))
         ENDDO
         ENDDO
C
C     ----- Given pressure and poloidal current profiles -----
C
      ELSEIF(IMDLEQF.EQ.1) THEN
         FJP=0.D0
         FJT1=0.D0
         FJT2=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQFPSI(PSIN,FPSI,DFPSI)
            PP(NTG,NSG)=PPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-2.D0*PI*BB*RR       *DFPSI
     &                                   /(RMU0*RMM(NTG,NSG))
            HJT2(NTG,NSG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                                   /(RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
     &                   *EXP(OTC*RMM(NTG,NSG)**2*AMP/2.D0)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP =FJP +HJP2(NTG,NSG)*DVOL
            FJT1=FJT1+HJT1(NTG,NSG)*DVOL
            FJT2=FJT2+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         TJ=(-FJT1-SQRT(FJT1**2+4.D0*FJT2*(RIP*1.D6-FJP)))
     &      /(2.D0*FJT2)
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT1(NTG,NSG)
     &                                +TJ*TJ*HJT2(NTG,NSG)
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQFPSI(PSIN,FPSI,DFPSI)
            CALL EQTPSI(PSIN,TPSI,DTPSI)
            TT(NTG,NSG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            RHO(NTG,NSG)=(PPSI*AMP/TPSI)
     &                  *EXP(OTC*RMM(NTG,NSG)**2*AMP/2.D0)
         ENDDO
         ENDDO 
C
C     ----- Given pressure and parallel current profiles -----
C
      ELSEIF(IMDLEQF.EQ.2) THEN
         FJP=0.D0
         FJT=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQJPSI(PSIN,HJPSID,HJPSI)
            CALL EQGPSI(PSIN,RUGQL,RUGJL,RUG1L,RUG2L,RUG3L,RUGBL)
            PP(NTG,NSG)=PPSI
            HJP1(NTG,NSG)=2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=HJPSI
            HJP2C=1.D0-RR**2*RUGJL**2/(RMM(NTG,NSG)**2*RUGBL)
            HJP2(NTG,NSG)=HJP2C*HJP1(NTG,NSG)
            HJT2(NTG,NSG)=RR*RUGJL/(RMM(NTG,NSG)*RUGBL)*HJT1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
C
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT=FJT+HJT2(NTG,NSG)*DVOL
C
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT2(NTG,NSG)
            TT(NTG,NSG)=BB*RR
            RHO(NTG,NSG)=0.D0
         ENDDO
         ENDDO
         RIP=(FJP+FJT)*1.D-6
C
C     ----- Given pressure and safety factor profiles ------
C
      ELSEIF(IMDLEQF.EQ.3) THEN
         FJP=0.D0
         FJT1=0.D0
         FJT2=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQQPSI(PSIN,QPSI,DQPSI)
            PP(NTG,NSG)=PPSI
            HJP1(NTG,NSG)=RMM(NTG,NSG)*DPPSI
            CALL FNFQT(PSIN,FQTL,DFQTL)
            TTL=QPSI/FQTL
            DTT=(DQPSI*FQTL-QPSI*DFQTL)/FQTL**2
            HJT1(NTG,NSG)=BB*RR*DTT/(RMU0*RMM(NTG,NSG))
            HJT2(NTG,NSG)=(TTL-BB*RR)*DTT/(RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
     &                   *EXP(OTC*RMM(NTG,NSG)**2*AMP/2.D0)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP =FJP +HJP2(NTG,NSG)*DVOL
            FJT1=FJT1+HJT1(NTG,NSG)*DVOL
            FJT2=FJT2+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         TJ=(-FJT1-SQRT(FJT1**2-4.D0*FJT2*(RIP*1.D6+FJP)))
     &      /(2.D0*FJT2)
C         TJ=(-RIP*1.D6-FJP)/FJT
         WRITE(6,'(A,1P4E12.4)') 'TJ,FJP,FJT1,FJT2=',
     &                            TJ,FJP,FJT1,FJT2
         TJ1=TJ
         IF(TJ.GT.0.01D0) TJ1=0.01D0
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT1(NTG,NSG)
     &                                +TJ*TJ*HJT2(NTG,NSG)
            PSIN=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQQPSI(PSIN,QPSI,DQPSI)
            CALL FNFQT(PSIN,FQTL,DFQTL)
            TT(NTG,NSG)=BB*RR+TJ1*(QPSI/FQTL-BB*RR)
            RHO(NTG,NSG)=(PPSI*AMP/TPSI)
     &                  *EXP(OTC*RMM(NTG,NSG)**2*AMP/2.D0)
         ENDDO
         ENDDO 
         STOP
      ENDIF
C
      RETURN
      END
C
C   ************************************************
C   **              Matrix Solver                 **
C   ************************************************
C
      SUBROUTINE EQSOLV
C
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION FJT(MLM),PSIOLD(NTGM,NSGM)
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIOLD(NTG,NSG)=PSI(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         I=(NSG-1)*NTGMAX
         DO NTG=1,NTGMAX
            FJT(I+NTG)=2.D0*PI*RMU0*HJT(NTG,NSG)
     &                 *SIGM(NSG)*RHOM(NTG)*RHOM(NTG)
         ENDDO
      ENDDO
C
C      DO NSG=1,2
C         DO NTG=1,NTGMAX
C            WRITE(6,'(2I5,1PE12.4)') NSG,NTG,HJT(NSG,NTG)
C         ENDDO
C      ENDDO
C      PAUSE
C
      DO I=1,3
         WRITE(6,'(1p5E12.4)') FJT(I),(Q(J,I),J=1,4*NTGMAX-1)
      ENDDO
      pause
C
      CALL BANDRD(Q,FJT,NTGMAX*NSGMAX,4*NTGMAX-1,MWM,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX EQSOLV: BANDRD ERROR: IERR = ',IERR
         ENDIF
C
       DO NSG=1,NSGMAX
          I=(NSG-1)*NTGMAX
          DO NTG=1,NTGMAX
             PSI(NTG,NSG)=FJT(I+NTG)
          ENDDO
       ENDDO
C     
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         DELPSI(NTG,NSG)=PSI(NTG,NSG)-PSIOLD(NTG,NSG)
      ENDDO
      ENDDO
      RETURN
      END
C
C   ******************************************
C   ***** Calculate Spline Coeff for Psi *****
C   ******************************************
C
      SUBROUTINE EQSETF
C
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION PSISX(NTGPM,NSGPM),PSITX(NTGPM,NSGPM)
      DIMENSION PSISTX(NTGPM,NSGPM)
C
C     ----- mesh extended in sigma (radius) and theta (periodic) -----
C
      SIGMX(1)=0.D0
      DO NSG=1,NSGMAX
         SIGMX(NSG+1)=SIGM(NSG)
      ENDDO
      SIGMX(NSGMAX+2)=1.D0
C
      THGMX(1)=0.D0
      DO NTG=1,NTGMAX
         THGMX(NTG+1)=THGM(NTG)
      ENDDO
      THGMX(NTGMAX+2)=2.D0*PI
C
      SUM=0.D0
      DO NTG=1,NTGMAX
         SUM=SUM+(9*PSI(NTG,1)-PSI(NTG,2))/8.D0
      ENDDO
      PSI1=SUM/NTGMAX
      DO NTG=1,NTGMAX+2
         PSIX(NTG,1)=PSI1
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIX(NTG+1,NSG+1)=PSI(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         PSIX(       1,NSG+1)=(9*PSI(     1,NSG)-PSI(       2,NSG))
     &                        /16.D0
     &                       +(9*PSI(NTGMAX,NSG)-PSI(NTGMAX-1,NSG))
     &                        /16.D0
         PSIX(NTGMAX+2,NSG+1)=PSIX(1,NSG+1)
      ENDDO
C
      DO NTG=1,NTGMAX+2
         PSIX(NTG,NSGMAX+2)=0.D0
      ENDDO
C
      CALL SPL2D(THGMX,SIGMX,PSIX,PSITX,PSISX,PSISTX,U,
     &           NTGPM,NTGMAX+2,NSGMAX+2,4,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for PSIX: IERR=',IERR
      RETURN
      END
C
C   ************************************************
C   **          sigma,theta to R,Z                **
C   ************************************************
C
      SUBROUTINE EQTORZ
C
      INCLUDE '../eq/eqcomc.inc'
      EXTERNAL EQFBND
C
      RMIN= RR-RB
      RMAX= RR+RB
      ZMIN=-RKAP*RB
      ZMAX= RKAP*RB
C
      SIG1=1.D0
      SIG2=0.95D0
C
      DTRG=(RMAX-RMIN)/(NRGMAX-1)
      DTZG=(ZMAX-ZMIN)/(NZGMAX-1)
      EPSZ=1.D-8
C
      DO NRG=1,NRGMAX
         RG(NRG)=RMIN+DTRG*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=ZMIN+DTZG*(NZG-1)
      ENDDO
C
      CALL EQSETF
C
      DO NRG=1,NRGMAX
      DO NZG=1,NZGMAX
         THL=ATAN2(ZG(NZG),RG(NRG)-RR)
         IF(THL.LT.0.D0) THL=THL+2.D0*PI
         ZBRF=TAN(THL)
         THDASH=ZBRENT(EQFBND,THL-1.0D0,THL+1.0D0,EPSZ)
         RHOL=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &               +RKAP**2*SIN(THDASH)**2)
         SIGL=SQRT((RG(NRG)-RR)**2+ZG(NZG)**2)/RHOL
         IF(SIGL.LT.1.D0) THEN
            PSIRZ(NRG,NZG)=PSIF(SIGL,THL)
         ELSE
            PSI1=PSIF(SIG1,THL)
            PSI2=PSIF(SIG2,THL)
            PSIRZ(NRG,NZG)=PSI2+(PSI1-PSI2)*(SIGL-SIG2)/(SIG1-SIG2)
         ENDIF
      ENDDO
      ENDDO
C
      RETURN
      END
C
C   ************************************************
C   **         Boundary shape function            **
C   ************************************************
C
      FUNCTION EQFBND(X)
C      
      INCLUDE 'eqcomc.inc'
C
      EQFBND=ZBRF*COS(X+RDLT*SIN(X))-RKAP*SIN(X)
      RETURN
      END
C
C   ************************************************
C   **      CALCULATE pp,tt,temp,omega,rho        **
C   ************************************************
C
      SUBROUTINE EQCALP
C
      INCLUDE '../eq/eqcomc.inc'
C
      IMDLEQF=MOD(MDLEQF,5)
      DPS=PSI0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
         PSIN=1.D0-PSIPS(NPS)/PSI0
         CALL EQPPSI(PSIN,PPSI,DPPSI)
         CALL EQJPSI(PSIN,HJPSID,HJPSI)
         CALL EQTPSI(PSIN,TPSI,DTPSI)
         CALL EQOPSI(PSIN,OMGPSI,DOMGPSI)
C
         IF (IMDLEQF.EQ.0) THEN
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQJPSI(PSIN,HJPSID,HJPSI)
            CALL EQTPSI(PSIN,TPSI,DTPSI)
            CALL EQOPSI(PSIN,OMGPSI,DOMGPSI)
            PPPS(NPS)=PPSI
            OMPS(NPS)=OMGPSI
            TTPS(NPS)=SQRT((2.D0*PI*BB*RR)**2
     &                  +2.D0*RMU0*RRC
     &                  *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                  *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
            TEPS(NPS)=TPSI/(AEE*1.D3)
         ELSEIF (IMDLEQF.EQ.1) THEN
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQFPSI(PSIN,FPSI,DFPSI)
            CALL EQOPSI(PSIN,OMGPSI,DOMGPSI)
            CALL EQTPSI(PSIN,TPSI,DTPSI)
            PPPS(NPS)=PPSI
            TTPS(NPS)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            OMPS(NPS)=OMGPSI
            TEPS(NPS)=TPSI/(AEE*1.D3)
         ELSEIF (IMDLEQF.EQ.2) THEN
            CALL EQPPSI(PSIN,PPSI,DPPSI)
            CALL EQJPSI(PSIN,HJPSID,HJPSI)
            TTPS(NPS)=SQRT(BB**2*RR**2
     &                    +2.D0*RMU0*RRC*(TJ*HJPSID-RRC*PPSI))
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ELSEIF (IMDLEQF.EQ.3) THEN
            CALL EQQPSI(PSIN,QPSI,DQPSI)
            CALL FNFQT(PSIN,FQTL,DFQTL)
            TTPS(NPS)=QPSI/FQTL
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ENDIF
      ENDDO
      RETURN
      END
C
C     ***** CALCULATE TTHETA *****
C
      SUBROUTINE EQCALT(TTH)
C
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION TTH(NUGM)
C
      DRHO=1.D0/DBLE(NUGMAX-1)
      TTH(NR)=2.D0*PI*BB*RR
      DO NU=NUGMAX,2,-1
         RHOP= NU   *DRHO
         RHON=(NU-1)*DRHO
         PSINP=RHOP**2
         PSINN=RHON**2
         CALL EQGPSI(PSINP,RUGQPL,RUGJP,RUG1P,RUG2P,RUG3P,RUGBP)
         CALL EQGPSI(PSINN,RUGQNL,RUGJN,RUG1N,RUG2N,RUG3N,RUGBN)
         RG5=0.5D0*(RG5P+RG5N)
         PQL=0.5D0*(PQLP+PQLN)
         HJL=0.5D0*(HJLP+HJLN)
         VA= 0.5D0*RMU*(PLP-PLN)/(BB**2*RG5)
         VB=-RMU*HJL*PSITS*(PSINP-PSINN)/(PQL*BB*RG5)
         TTH(NU-1)=((1.D0+VA)*TTH(NR)-VB)/(1.D0-VA)
      ENDDO
      RETURN
      END
C
C     ***** CALCULATE METRIC *****
C
      SUBROUTINE EQLPSI(PSIN,PL,HJL,PQL,RG5)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL EQPPSI(PSIN,PL,DPPSI)
      CALL EQJPSI(PSIN,HJL,DHJPSI)
      PQL=FNQPS(PSIN)
      RG5=FNG5S(PSIN)
      RETURN
      END
C
C     ***** CALCULATE Q *****
C
      SUBROUTINE EQPSIQ
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION XA(NNM),YA(2,NNM)
      DIMENSION DERIV(NPSM)
      EXTERNAL PSIAX,EQDERV
C
      REDGE=ZBRENT(PSIAX,RR,RR+RB,1.D-8)
      DR=(REDGE-RAXIS)/(NRMAX-1)
      NMAX=200
      IF(NMAX.GT.NNM) NMAX=NNM
C
      PSQ(1)=PSIG(RAXIS,ZAXIS)
C
      DO NR=2,NRMAX
         RINIT=RAXIS+DR*(NR-1)
         ZINIT=ZAXIS
         PSQ(NR)=PSIG(RINIT,ZINIT)
C
         CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
         SUML=0.D0
         SUMV=0.D0
         SUMQ=0.D0
         SUMR1=0.D0
         SUMR2=0.D0
         SUMR3=0.D0
         DO N=2,NA
            CALL EQPSID(YA(1,N),YA(2,N),DPSIDR,DPSIDZ)
            R=YA(1,N)
            H=XA(N)-XA(N-1)
            BPRL=SQRT(DPSIDR**2+DPSIDZ**2)
C
            SUMV=SUMV+H*R/BPRL
            SUMQ=SUMQ+H/(BPRL*R)
            SUMR1=SUMR1+H*R
            SUMR2=SUMR2+H*R*BPRL
            SUMR3=SUMR3+H*BPRL/R
         ENDDO
         FVT(NR)=SUMV/(2.D0*PI)
         FQT(NR)=SUMQ/(2.D0*PI)
         AVR1(NR)=SUMR1/SUMV
         AVR2(NR)=SUMR2/SUMV
C         WRITE(6,*) NR,PSQ(NR),FQT(NR)
      ENDDO
C
      FQT(1)=FQT(2)-(FQT(3)-FQT(2))*(PSQ(2)-PSQ(1))/(PSQ(3)-PSQ(2))
      WRITE(6,'(1P3E12.4)') FQT(1),FQT(2),FQT(NRMAX)
C
      CALL SPL1D(PSQ,FQT,DERIV,UFQT,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for FQT: IERR=',IERR
      RETURN
      END
C
C   *******************************************
C   ***** Calculate Psi at (sigma, theta) *****
C   *******************************************
C
      FUNCTION PSIF(RSIG,RTHG)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DF(RTHG,RSIG,PSIL,THGMX,SIGMX,U,
     &            NTGPM,NTGMAX+2,NSGMAX+2,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PSIF: SPL2DF ERROR : IERR=',IERR
      PSIF=PSIL
      RETURN
      END
