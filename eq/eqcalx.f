C     $Id$
C
C     ****** FREE BOUNDARY EQUILIBRIUM SOLVER ******
C
      SUBROUTINE EQCALX(ID,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
C
      DRG=(RGMAX-RGMIN)/(NRGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=RGMIN+DRG*(NRG-1)
      ENDDO
C
      DZG=(ZGMAX-ZGMIN)/(NZGMAX-1)
      DO NZG=1,NZGMAX
         ZG(NZG)=ZGMIN+DZG*(NZG-1)
      ENDDO
C
      IF(ID.EQ.0) CALL EQPSIX
C
      MWMAX=8*(NRGMAX+2)-1
      MLMAX=4*NRGMAX*NZGMAX
      NBND=4*(NRGMAX+2)
C
      DO NLOOP=1,NLPMAX
         CALL EQCALFMA
         CALL EQCALFVB(ID)
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIXOLD(NRG,NZG)=PSIX(NRG,NZG)
         ENDDO
         ENDDO
C
         CALL BANDRD(FMA,FVB,MLMAX,MWMAX,MWM,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'BANDRD ERROR: IERR =',IERR
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            N=4*((NZG-1)*NRGMAX+NRG-1)+1
            PSIX(NRG,NZG)=FVB(N)
            UPSIX(1,NRG,NZG)=FVB(N)
            UPSIX(2,NRG,NZG)=FVB(N+1)
            UPSIX(3,NRG,NZG)=FVB(N+2)
            UPSIX(4,NRG,NZG)=FVB(N+3)
         ENDDO
         ENDDO
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            DELPSIX(NRG,NZG)=PSIX(NRG,NZG)-PSIXOLD(NRG,NZG)
            HJTRZ(NRG,NZG)=FJRZ(NRG,NZG)
         ENDDO
         ENDDO
C
         DO I=1,21
            R=2.9D0+0.01D0*(I-1)
            Z=0.D0
            PSIL=PSIXH(R,Z)
            CALL PSIXHD(R,Z,DPSIDR,DPSIDZ)
            WRITE(6,'(A,1P4E12.4)') 
     &           'R,PSI,DR,DZ=',R,PSIL,DPSIDR,DPSIDZ
         ENDDO
C
         DO I=1,11
            X=0.1D0*(I-1)
            WRITE(6,'(A,1P5E12.4)') 
     &           'R,H0,H1,DH0,DH1=',X,HRMT0(X),HRMT1(X),
     &                                DHRMT0(X),DHRMT1(X)
         ENDDO
C
         DO I=1,21
            R=3.D0
            Z=-0.1D0+0.01D0*(I-1)
            PSIL=PSIXH(R,Z)
            CALL PSIXHD(R,Z,DPSIDR,DPSIDZ)
            WRITE(6,'(A,1P4E12.4)') 
     &           'Z,PSI,DR,DZ=',Z,PSIL,DPSIDR,DPSIDZ
         ENDDO
C
C         DO I=1,10
C            READ(5,*) R,Z
C            CALL PSIXHD(R,Z,D1,D2)
C            WRITE(6,'(A,1P5E12.4)') 'R,Z,PSI,DR,DZ=',
C     &           R,Z,PSIXH(R,Z),D1,D2
C         ENDDO
C
         CALL EQAXISX(IERR)
C
         PSI0=PSIX(1,1)
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSI0=MIN(PSI0,PSIX(NRG,NZG))
         ENDDO
         ENDDO
         PSIPA=-PSI0
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=PSIX(NRG,NZG)-PSI0
            HJTRZ(NRG,NZG)=FJRZ(NRG,NZG)
         ENDDO
         ENDDO
C
         SUM=0.D0
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            SUM=SUM+DELPSIX(NRG,NZG)**2
         ENDDO
         ENDDO
         WRITE(6,*) 'SUM=',SUM
         IF(SUM.LT.EPSEQ) GOTO 1000
      ENDDO
      WRITE(6,*) 'XX NO CONVERGENCE IN EQCALX'
C
 1000 CONTINUE
      DPS=1.D0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPNL=DPS*(NPS-1)
         PSIPS(NPS)=PSIPA*PSIPNL
         CALL EQPPSI(PSIPNL,PPSI,DPPSI)
         CALL EQFPSI(PSIPNL,FPSI,DFPSI)
         PPPS(NPS)=PPSI
         TTPS(NPS)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         WRITE(6,'(I5,1P5E12.4)')
     &        NPS,PSIPNL,PPSI,FPSI,PPPS(NPS),TTPS(NPS)
      ENDDO
C
      CALL EQGOUT(1)
C
      RETURN
      END
C
C     ******************************************
C     *     Set FEM element integral table     *
C     ******************************************
C
      SUBROUTINE EQTABL
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION K1AL(4,4),M1AL(4,4),L1AL(4,4),N1AL(4,4)
C
      DATA K1AL/1,1,2,2,
     &          1,1,2,2,
     &          3,3,4,4,
     &          3,3,4,4/
      DATA M1AL/1,2,1,2,
     &          3,4,3,4,
     &          1,2,1,2,
     &          3,4,3,4/
      DATA L1AL/1,1,2,2,
     &          1,1,2,2,
     &          3,3,4,4,
     &          3,3,4,4/
      DATA N1AL/1,2,1,2,
     &          3,4,3,4,
     &          1,2,1,2,
     &          3,4,3,4/
C
      DO J=1,4
      DO I=1,4
         K1A(I,J)=K1AL(I,J)
         L1A(I,J)=L1AL(I,J)
         M1A(I,J)=M1AL(I,J)
         N1A(I,J)=N1AL(I,J)
      ENDDO
      ENDDO
C
      RK(1,1,1,1)= 3.D0/5.D0
      RK(1,1,1,2)= 3.D0/5.D0
      RK(1,1,2,1)=-3.D0/5.D0
      RK(1,1,2,2)=-3.D0/5.D0
      RK(1,2,1,1)=-3.D0/5.D0
      RK(1,2,1,2)=-3.D0/5.D0
      RK(1,2,2,1)= 3.D0/5.D0
      RK(1,2,2,2)= 3.D0/5.D0
c
      RK(2,1,1,1)= 1.D0/70.D0
      RK(2,1,1,2)=-4.D0/35.D0
      RK(2,1,2,1)=-4.D0/35.D0
      RK(2,1,2,2)= 1.D0/70.D0
      RK(2,2,1,1)=-1.D0/70.D0
      RK(2,2,1,2)= 4.D0/35.D0
      RK(2,2,2,1)= 4.D0/35.D0
      RK(2,2,2,2)=-1.D0/70.D0
c
      RK(3,1,1,1)= 1.D0/70.D0
      RK(3,1,1,2)=-4.D0/35.D0
      RK(3,1,2,1)=-1.D0/70.D0
      RK(3,1,2,2)= 4.D0/35.D0
      RK(3,2,1,1)=-4.D0/35.D0
      RK(3,2,1,2)= 1.D0/70.D0
      RK(3,2,2,1)= 4.D0/35.D0
      RK(3,2,2,2)=-1.D0/70.D0
c
      RK(4,1,1,1)= 43.D0/420.D0
      RK(4,1,1,2)= 13.D0/420.D0
      RK(4,1,2,1)=-1.D0/60.D0
      RK(4,1,2,2)=-1.D0/60.D0
      RK(4,2,1,1)=-1.D0/60.D0
      RK(4,2,1,2)=-1.D0/60.D0
      RK(4,2,2,1)=13.D0/420.D0
      RK(4,2,2,2)=43.D0/420.D0
C
      RKK(1,1,1,1)=-9.D0/70.D0
      RKK(1,1,1,2)= 9.D0/70.D0
      RKK(1,1,2,1)= 9.D0/70.D0
      RKK(1,1,2,2)=-9.D0/70.D0
      RKK(1,2,1,1)= 9.D0/70.D0
      RKK(1,2,1,2)=-9.D0/70.D0
      RKK(1,2,2,1)=-9.D0/70.D0
      RKK(1,2,2,2)= 9.D0/70.D0
c
      RKK(2,1,1,1)= 1.D0/140.D0
      RKK(2,1,1,2)=-3.D0/140.D0
      RKK(2,1,2,1)= 3.D0/140.D0
      RKK(2,1,2,2)=-1.D0/140.D0
      RKK(2,2,1,1)=-1.D0/140.D0
      RKK(2,2,1,2)= 3.D0/140.D0
      RKK(2,2,2,1)=-3.D0/140.D0
      RKK(2,2,2,2)= 1.D0/140.D0
c
      RKK(3,1,1,1)= 1.D0/140.D0
      RKK(3,1,1,2)=-3.D0/140.D0
      RKK(3,1,2,1)=-1.D0/140.D0
      RKK(3,1,2,2)= 3.D0/140.D0
      RKK(3,2,1,1)= 3.D0/140.D0
      RKK(3,2,1,2)=-1.D0/140.D0
      RKK(3,2,2,1)=-3.D0/140.D0
      RKK(3,2,2,2)= 1.D0/140.D0
c
      RKK(4,1,1,1)=-1.D0/120.D0
      RKK(4,1,1,2)= 1.D0/1680.D0
      RKK(4,1,2,1)= 1.D0/840.D0
      RKK(4,1,2,2)=-1.D0/840.D0
      RKK(4,2,1,1)= 1.D0/840.D0
      RKK(4,2,1,2)=-1.D0/840.D0
      RKK(4,2,2,1)=-1.D0/168.D0
      RKK(4,2,2,2)= 1.D0/120.D0
C
      RL(1,1,1,1)= 43.D0/140.D0
      RL(1,1,1,2)= 9.D0/140.D0
      RL(1,1,2,1)= 9.D0/140.D0
      RL(1,1,2,2)= 9.D0/140.D0
      RL(1,2,1,1)= 9.D0/140.D0
      RL(1,2,1,2)= 9.D0/140.D0
      RL(1,2,2,1)= 9.D0/140.D0
      RL(1,2,2,2)= 43.D0/140.D0
c
      RL(2,1,1,1)=-97.D0/2520.D0
      RL(2,1,1,2)=-1.D0/72.D0
      RL(2,1,2,1)=43.D0/2520.D0
      RL(2,1,2,2)= 1.D0/72.D0
      RL(2,2,1,1)=-1.D0/72.D0
      RL(2,2,1,2)=-43.D0/2520.D0
      RL(2,2,2,1)= 1.D0/72.D0
      RL(2,2,2,2)=97.D0/2520.D0
c
      RL(3,1,1,1)= -97.D0/2520.D0
      RL(3,1,1,2)=-1.D0/72.D0
      RL(3,1,2,1)=-1.D0/72.D0
      RL(3,1,2,2)=-43.D0/2520.D0
      RL(3,2,1,1)=43.D0/2520.D0
      RL(3,2,1,2)= 1.D0/72.D0
      RL(3,2,2,1)= 1.D0/72.D0
      RL(3,2,2,2)=97.D0/2520.D0
C
      RL(4,1,1,1)= 2.D0/315.D0
      RL(4,1,1,2)= 1.D0/315.D0
      RL(4,1,2,1)=-1.D0/280.D0
      RL(4,1,2,2)=-1.D0/280.D0
      RL(4,2,1,1)=-1.D0/280.D0
      RL(4,2,1,2)=-1.D0/280.D0
      RL(4,2,2,1)= 1.D0/315.D0
      RL(4,2,2,2)= 2.D0/315.D0
C
      RLL(1,1,1,1)=-97.D0/2520.D0
      RLL(1,1,1,2)= 43.D0/2520.D0
      RLL(1,1,2,1)=-1.D0/72.D0
      RLL(1,1,2,2)= 1.D0/72.D0
      RLL(1,2,1,1)=-1.D0/72.D0
      RLL(1,2,1,2)= 1.D0/72.D0
      RLL(1,2,2,1)=-43.D0/2520.D0
      RLL(1,2,2,2)= 97.D0/2520.D0
C
      RLL(2,1,1,1)= 2.D0/315.D0
      RLL(2,1,1,2)=-1.D0/280.D0
      RLL(2,1,2,1)= -1.D0/280.D0
      RLL(2,1,2,2)=1.D0/350.D0
      RLL(2,2,1,1)= 1.D0/350.D0
      RLL(2,2,1,2)=-1.D0/280.D0
      RLL(2,2,2,1)=-1.D0/280.D0
      RLL(2,2,2,2)= 2.D0/315.D0
C
      RLL(3,1,1,1)= 2.D0/315.D0
      RLL(3,1,1,2)=-1.D0/280.D0
      RLL(3,1,2,1)= 1.D0/315.D0
      RLL(3,1,2,2)=-1.D0/280.D0
      RLL(3,2,1,1)=-1.D0/280.D0
      RLL(3,2,1,2)= 1.D0/3150.D0
      RLL(3,2,2,1)=-1.D0/280.D0
      RLL(3,2,2,2)= 2.D0/315.D0
C
      RLL(4,1,1,1)=-1.D0/840.D0
      RLL(4,1,1,2)= 1.D0/1260.D0
      RLL(4,1,2,1)= 1.D0/1260.D0
      RLL(4,1,2,2)=-1.D0/1260.D0
      RLL(4,2,1,1)= 1.D0/1260.D0
      RLL(4,2,1,2)=-1.D0/1260.D0
      RLL(4,2,2,1)=-1.D0/1260.D0
      RLL(4,2,2,2)= 1.D0/840.D0
C
      RM(1,1,1)= 13.D0/ 35.D0
      RM(1,1,2)=  9.D0/ 70.D0
      RM(1,2,1)=  9.D0/ 70.D0
      RM(1,2,2)= 13.D0/ 35.D0
      RM(2,1,1)=-11.D0/210.D0
      RM(2,1,2)= 13.D0/420.D0
      RM(2,2,1)=-13.D0/420.D0
      RM(2,2,2)= 11.D0/210.D0
      RM(3,1,1)=-11.D0/210.D0
      RM(3,1,2)=-13.D0/420.D0
      RM(3,2,1)= 13.D0/420.D0
      RM(3,2,2)= 11.D0/210.D0
      RM(4,1,1)=  1.D0/105.D0
      RM(4,1,2)=- 1.D0/140.D0
      RM(4,2,1)=- 1.D0/140.D0
      RM(4,2,2)=  1.D0/105.D0
C
      RN(1,1,1)= 6.D0/ 5.D0
      RN(1,1,2)=-6.D0/ 5.D0
      RN(1,2,1)=-6.D0/ 5.D0
      RN(1,2,2)= 6.D0/ 5.D0
      RN(2,1,1)=-1.D0/10.D0
      RN(2,1,2)=-1.D0/10.D0
      RN(2,2,1)= 1.D0/10.D0
      RN(2,2,2)= 1.D0/10.D0
      RN(3,1,1)=-1.D0/10.D0
      RN(3,1,2)= 1.D0/10.D0
      RN(3,2,1)=-1.D0/10.D0
      RN(3,2,2)= 1.D0/10.D0
      RN(4,1,1)= 2.D0/15.D0
      RN(4,1,2)=-1.D0/30.D0
      RN(4,2,1)=-1.D0/30.D0
      RN(4,2,2)= 2.D0/15.D0
      RETURN
      END
C
C   ************************************************ 
C   **              Initial Psi                 **
C   ************************************************
C
      SUBROUTINE EQPSIX
C
      INCLUDE '../eq/eqcomx.inc'
      DIMENSION DERIV(NRVM)
C
      RAXIS=RR
      ZAXIS=0.0D0
C
C     --- assuming elliptic crosssection, flat current profile ---
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MOD(MDLEQFL,5).EQ.4) THEN
         PSITA=PI*RKAP*RA**2*BB
         PSIPA=PSITA/QA
         PSI0=-PSIPA
      ELSE
         PSI0=-0.5D0*RMU0*RIP*1.D6*RR
         PSIPA=-PSI0
         PSITA=PI*RKAP*RA**2*BB
      ENDIF
C
      DRHO=1.D0/(NRVMAX-1)
      DO NRV=1,NRVMAX
         RHOL=(NRV-1)*DRHO
         PSIPNV(NRV)=RHOL*RHOL
         PSIPV(NRV)=PSIPA*RHOL*RHOL
         PSITV(NRV)=PSITA*RHOL*RHOL
         QPV(NRV)=PSITA/PSIPA
         TTV(NRV)=2.D0*PI*RR*BB
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         RHON2=((RG(NRG)-RR)**2+(ZG(NZG)/RKAP)**2)/RA**2
         PSIX(NRG,NZG)=PSI0*(1-RHON2)
         DELPSIX(NRG,NZG)=0.D0
      ENDDO
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C      WRITE(6,'(A,I5,1P4E12.4)')
C     &     ('NRV:',NRV,PSIPNV(NRV),PSITV(NRV),QPV(NRV),TTV(NRV),
C     &      NRV=1,NRVMAX)
C
      RETURN
      END
C
C     ********************************************
C     *     Calculate FEM coefficient matrix     *
C     ********************************************
C
      SUBROUTINE EQCALFMA
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION FACTK(4),FACTL(4),FACTM(4),FACTN(4)
C
      DRG=(RGMAX-RGMIN)/(NRGMAX-1)
      DZG=(ZGMAX-ZGMIN)/(NZGMAX-1)
C
      DO N=1,MLMAX
      DO M=1,MWMAX
         FMA(M,N)=0
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         DRG=RG(NRG)-RG(NRG+1)
         FACTK(1)=1.D0/DRG**2
         FACTK(2)=1.D0/DRG
         FACTK(3)=1.D0/DRG
         FACTK(4)=1.D0
         FACTL(1)=1.D0
         FACTL(2)=1.D0*DRG
         FACTL(3)=1.D0*DRG
         FACTL(4)=1.D0*DRG**2
C
         DZG=ZG(NZG)-ZG(NZG+1)
         FACTM(1)=1.D0
         FACTM(2)=1.D0*DZG
         FACTM(3)=1.D0*DZG
         FACTM(4)=1.D0*DZG**2
         FACTN(1)=1.D0/DZG**2
         FACTN(2)=1.D0/DZG
         FACTN(3)=1.D0/DZG
         FACTN(4)=1.D0
C
         DO I1=1,2
         DO I2=1,2
         DO I3=1,2
         DO J1=1,2
         DO J2=1,2
         DO K=1,4
         DO L=1,4
            N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
            M=4*((NZG-1)*NRGMAX+NRG-1)+L+4*(I2-1)+4*NRGMAX*(J2-1)
C
            K1=K1A(K,L)
            M1=M1A(K,L)
            L1=L1A(K,L)
            N1=N1A(K,L)
            FMA(NBND+M-N,N)=FMA(NBND+M-N,N)
     &                     +(RK(K1,I1,I2,I3)/RG(NRG+I3-1)
     &                      -RKK(K1,I1,I2,I3)/RG(NRG+I3-1)**2)
     &                                   *FACTK(K1)
     &                      *RM(M1,J1,J2)*FACTM(M1)
     &                     +(RL(L1,I1,I2,I3)/RG(NRG+I3-1)
     &                      -RLL(L1,I1,I2,I3)/RG(NRG+I3-1)**2)
     &                                   *FACTL(L1)
     &                      *RN(N1,J1,J2)*FACTN(N1) 
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX,NZGMAX-1
      DO NRG=1,NRGMAX
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         DO MW=1,MWMAX
            FMA(MW,N)=0.D0
         ENDDO
         FMA(NBND,N)=1.D0
      ENDDO
      ENDDO
C
      DO NRG=1,NRGMAX,NRGMAX-1
      DO NZG=2,NZGMAX-1
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         DO MW=1,MWMAX
            FMA(MW,N)=0.D0
         ENDDO
         FMA(NBND,N)=1.D0
      ENDDO
      ENDDO
C
C      OPEN(16,FILE='eqxdata',FORM='FORMATTED')
C      DO 9999 ML=1,MLMAX
C         WRITE(16,*) 'ML=',ML
C         WRITE(16,'(10F7.2)') (FMA(MW,ML),MW=1,MWMAX)
C 9999 CONTINUE
C      CLOSE(16)
C
      RETURN
      END
C
C
C     **********************************
C     *     Calulate FEM RHS vector    *
C     **********************************
C
      SUBROUTINE EQCALFVB(ID)
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION FACTM(4),FACTH(4)
      DIMENSION RJ(4),PSIBRZ(4)
      DIMENSION HJ0(NRGM,NZGM),HJ1(NRGM,NZGM),HJ2(NRGM,NZGM)
C
      DRG=(RGMAX-RGMIN)/(NRGMAX-1)
      DZG=(ZGMAX-ZGMIN)/(NZGMAX-1)
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
        FJRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         IF(PSIX(NRG,NZG)*PSI0.GT.0.D0) THEN
            PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            PPRZ(NRG,NZG)=PPSI
         ELSE
            PPRZ(NRG,NZG)=0.D0
         ENDIF
      ENDDO
      ENDDO
C
      FJ0=0.D0
      FJ1=0.D0
      FJ2=0.D0
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &        ZG(NZG).LE.ZLIMP.AND.
     &        ZG(NZG).GE.ZLIMM) THEN
            PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
            HJ1(NRG,NZG)=-2.D0*PI*BB*RR*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
            HJ2(NRG,NZG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
         ELSE
            HJ0(NRG,NZG)=0.D0
            HJ1(NRG,NZG)=0.D0
            HJ2(NRG,NZG)=0.D0
         ENDIF
C
         DVOL=DRG*DZG
         FJ0=FJ0+HJ0(NRG,NZG)*DVOL
         FJ1=FJ1+HJ1(NRG,NZG)*DVOL
         FJ2=FJ2+HJ2(NRG,NZG)*DVOL
      ENDDO
      ENDDO
C
      IF(FJ1.GT.0.D0) THEN
         TJ=(-FJ1+SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
      ELSE
         TJ=(-FJ1-SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
      ENDIF
      WRITE(6,'(A,1P4E12.4)') 'FJ0,FJ1,FJ2,TJ=',FJ0,FJ1,FJ2,TJ
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         FJRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)+TJ*HJ1(NRG,NZG)
     &                                        +TJ*TJ*HJ2(NRG,NZG))
         IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &        ZG(NZG).LE.ZLIMP.AND.
     &        ZG(NZG).GE.ZLIMM) THEN
            PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            TTRZ(NRG,NZG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSE
            TTRZ(NRG,NZG)=2.D0*PI*BB*RR
         ENDIF
      ENDDO
      ENDDO
C
      DO NPFC=1,NPFCMAX
         NZC=INT((ZPFC(NPFC)-ZGMIN)/DZG)+1
         IF(NZC.GE.1.AND.NZC.LT.NZGMAX) THEN
C            DO NRG=1,NRGMAX
C            DO NZG=1,NZC+2
C               FJRZ(NRG,NZG)=0.D0
C            ENDDO
C            ENDDO
            DO I=1,3
               IF(I.EQ.1) THEN
                  RPFCL=RPFC(NPFC)
                  RIPFCL=RIPFC(NPFC)
               ELSEIF(I.EQ.2) THEN
                  RPFCL=RPFC(NPFC)-WPFC(NPFC)
                  RIPFCL=-0.5D0*RIPFC(NPFC)
               ELSEIF(I.EQ.3) THEN
                  RPFCL=RPFC(NPFC)+WPFC(NPFC)
                  RIPFCL=-0.5D0*RIPFC(NPFC)
               ENDIF
               NRC=INT((RPFCL-RGMIN)/DRG)+1
               IF(NRC.GE.1.AND.NRC.LT.NRGMAX) THEN
                  FACTR=(RPFCL-RG(NRC))/DRG
                  FACTZ=(ZPFC(NPFC)-ZG(NZC))/DZG
                  FACTRC=1.D0-FACTR
                  FACTZC=1.D0-FACTZ
C
                  FJRZ(NRC  ,NZC  )=FJRZ(NRC  ,NZC  )
     &                             +FACTRC*FACTZC*RIPFCL/(DRG*DZG)
                  FJRZ(NRC+1,NZC  )=FJRZ(NRC+1,NZC  )
     &                             +FACTR *FACTZC*RIPFCL/(DRG*DZG)
                  FJRZ(NRC  ,NZC+1)=FJRZ(NRC  ,NZC+1)
     &                             +FACTRC*FACTZ *RIPFCL/(DRG*DZG)
                  FJRZ(NRC+1,NZC+1)=FJRZ(NRC+1,NZC+1)
     &                             +FACTR *FACTZ *RIPFCL/(DRG*DZG)
               ELSE
                  WRITE(6,*) 'XX (RDC) OUT OF REGION'
               ENDIF
            ENDDO
         ELSE
            WRITE(6,*) 'XX (ZPFC) OUT OF REGION: ZPFC=',ZPFC(NPFC)
         ENDIF
      ENDDO
C
      DO N=1,MLMAX
         FVB(N)=0
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         DRG=RG(NRG)-RG(NRG+1)
         FACTH(1)=1.D0
         FACTH(2)=1.D0*DRG
         FACTH(3)=1.D0*DRG
         FACTH(4)=1.D0*DRG**2
C
         DZG=ZG(NZG)-ZG(NZG+1)
         FACTM(1)=1.D0
         FACTM(2)=1.D0*DZG
         FACTM(3)=1.D0*DZG
         FACTM(4)=1.D0*DZG**2
C
         DO I1=1,2
         DO I2=1,2
         DO J1=1,2
         DO J2=1,2
            RJ(1)=FJRZ(NRG+I2-1,NZG+J2-1)
            IF(NRG+I2-1.GT.1.AND.
     &         NRG+I2-1.LT.NRGMAX) THEN
               RJ(2)=(FJRZ(NRG+I2,NZG+J2-1)-FJRZ(NRG+I2-2,NZG+J2-1))
     &              /(2.D0*DRG)
            ELSE  
               RJ(2)=0
            ENDIF
            IF(NZG+J2-1.GT.1.AND.
     &         NZG+J2-1.LT.NZGMAX) THEN
               RJ(3)=(FJRZ(NRG+I2-1,NZG+J2)-FJRZ(NRG+I2-1,NZG+J2-2))
     &              /(2.D0*DZG)
            ELSE  
               RJ(3)=0
            ENDIF
            RJ(4)=0
            DO K=1,4
            DO L=1,4
               N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
C
               K1=K1A(K,L)
               M1=M1A(K,L)
               FVB(N)=FVB(N)+RM(K1,I1,I2)*FACTH(K1)
     &                      *RM(M1,J1,J2)*FACTM(M1)*RMU0*1.D6*RJ(L)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C      WRITE(6,*) (FVB(ML),ML=1,MLMAX)
C
      DO NZG=1,NZGMAX,NZGMAX-1
      DO NRG=1,NRGMAX
         R1=RG(NRG)
         Z1=ZG(NZG)
         CALL EQPSIB(R1,Z1,PSIBRZ)
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         FVB(N)=PSIBRZ(K)
      ENDDO
      ENDDO
C
      DO  NRG=1,NRGMAX,NRGMAX-1
      DO NZG=2,NZGMAX-1
         R1=RG(NRG)
         Z1=ZG(NZG)
         CALL EQPSIB(R1,Z1,PSIBRZ)
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         FVB(N)=PSIBRZ(K)
      ENDDO
      ENDDO
      RETURN
      END
C
C     ****** BOUNDARY PSI ON WALL ******
C
      SUBROUTINE EQPSIB(R,Z,PSIBRZ)
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION PSIBRZ(4),PSIV(4,0:5)
C
      R01=RR
      R02=RR**2
      R03=RR**3
      R04=RR**4
      R05=RR**5
C
      PSIV(1,0)= 1.D0
      PSIV(2,0)= 0.D0
      PSIV(3,0)= 0.D0
      PSIV(4,0)= 0.D0
C
      PSIV(1,1)=(R**2-R02)/(2.D0*R01)
      PSIV(2,1)= R/R01
      PSIV(3,1)= 0.D0
      PSIV(4,1)= 0.D0
C
      PSIV(1,2)=-(R**2*Z**2
     &           -(R**2-R02)**2/4.D0)/R02
      PSIV(2,2)=-(2.D0*R*Z**2
     &           -R*(R**2-R02))/R02
      PSIV(3,2)=-2.D0*R**2*Z/R02
      PSIV(4,2)=-4.D0*R*Z/R02
C
      PSIV(1,3)=(R**2*Z**4
     &          -1.5D0*(R**2-R02)*R**2*Z**2
     &          +0.125D0*(R**2-R02)**3)/R03
      PSIV(2,3)=(2.D0*R*Z**4
     &          -3.D0*R**3.D0*Z**2
     &          -3.D0*R*(R**2-R02)*Z**2
     &          +0.75D0*R*(R**2-R02)**2)/R03
      PSIV(3,3)=(4.D0*R**2*Z**3
     &          -3.D0*(R**2-R03)*R**2*Z)/R03
      PSIV(4,3)=(8.D0*R*Z**3
     &          -6.D0*R**3*Z
     &          -6.D0*(R**2-R02)*R*Z)/R03
C
      PSIV(1,4)=-(0.2D0*R**2*Z**6
     &           -0.25D0*(3.D0*R**2-2.D0*R02)*R**2*Z**4
     &           +0.375D0*(R**2-R02)**2*R**2*Z**2
     &           -0.015625D0*(R**2-R02)**4)/R04
      PSIV(2,4)=-(0.4D0*R*Z**6-1.5D0*R**3*Z**4
     &           -0.5D0*R*(3.D0*R**2-2.D0*R02)*Z**4
     &           +1.5D0*R**3*(R**2-R02)*Z**2
     &           +0.75D0*R*(R**2-R02)**2*Z**2
     &           -0.125D0*R*(R**2-R02)**3)/R04
      PSIV(3,4)=-(1.2D0*R**2*Z**5
     &           -(3.D0*R**2-2.D0*R02)*R**2*Z**3
     &           +0.75D0*(R**2-R02)**2*R**2*Z)/R04
      PSIV(4,4)= (2.4D0*R*Z**5
     &           -6.D0*R**3*Z**3
     &           -2.D0*(3.D0*R**2-2.D0*R02)*Z**3
     &           +3.D0*R**3*(R**2-R02)*Z
     &           +1.5D0*R*(R**2-R02)*Z)/R04
C
      PSIV(1,5)= (4.D0/7.D0*R**2*Z**8
     &           -2.D0*(2.D0*R**2-R02)*R**2*Z**6
     &           +2.5D0*(2.D0*R**2-R02)*(R**2-R02)*R**2*Z**4
     &           -1.25D0*(R**2-R02)**3*R**2*Z**2
     &           +0.03125D0*(R**2-R02)**5)/R05
      PSIV(2,5)= (8.D0/7.D0*R*Z**8
     &           -8*R**3*Z**6
     &           -4.D0*R*(2.D0*R**2-R02)*Z**6
     &           +5.D0*R**3*(4.D0*R**2-3.D0*R02)*Z**4
     &           +5.D0*R*(2.D0*R**2-R02)*(R**2-R02)*Z**4
     &           -7.5D0*R**3*(R**2-R02)**2*Z**2
     &           -2.5D0*R*(R**2-R02)**3*Z**2
     &           +0.3125D0*R*(R**2-R02)**4)/R05
      PSIV(3,5)= (32.D0/7.D0*R**2*Z**7
     &           -12.D0*(2.D0*R**2-R02)*R**2*Z**5
     &           +10.D0*(2.D0*R**2-R02)*(R**2-R02)*R**2*Z**3
     &           -2.5D0*(R**2-R02)**3*R**2*Z)/R05
      PSIV(1,5)= (64.D0/7.D0*R*Z**7
     &           -48.D0*R**3*Z**5
     &           -24.D0*R*(2.D0*R**2-R02)*Z**5
     &           +20.D0*R**3*(4.D0*R**2-3.D0*R02)*Z**3
     &           +20.D0*R*(2.D0*R**2-R02)*(R**2-R02)*Z**3
     &           -15.D0*R**3*(R**2-R02)**2*Z
     &           -5.D0*R*(R**2-R02)**3*Z)/R05
C
      PSIBRZ(1)=0.D0
      PSIBRZ(2)=0.D0
      PSIBRZ(3)=0.D0
      PSIBRZ(4)=0.D0
      DO I=0,5
         PSIBRZ(1)=PSIBRZ(1)+2.D0*PI*PSIV(1,I)*PSIB(I)
         PSIBRZ(2)=PSIBRZ(2)+2.D0*PI*PSIV(2,I)*PSIB(I)
         PSIBRZ(3)=PSIBRZ(3)+2.D0*PI*PSIV(3,I)*PSIB(I)
         PSIBRZ(4)=PSIBRZ(4)+2.D0*PI*PSIV(4,I)*PSIB(I)
      ENDDO
C
      RETURN
      END
C
C     ***** CALCULATE MAGNETIC AXIS AND EDGE *****
C
      SUBROUTINE EQAXISX(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      EXTERNAL PSIXHD,PSIXZ0
C
      IERR=0
C
C     ----- calculate position of magnetic axis -----
C
      DELT=1.D-8
      EPS=1.D-4
      ILMAX=40
      LIST=1
      RINIT=RAXIS
      ZINIT=ZAXIS
      RSAVE=RAXIS
      ZSAVE=ZAXIS
      CALL NEWTN(PSIXHD,RINIT,ZINIT,RAXIS,ZAXIS,
     &           DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5,1P2E12.4)')
     &        'XX EQAXISX: NEWTN ERROR: IER=',IER,RSAVE,ZSAVE
         WRITE(6,'(A)') 'XX EQAXISX: AXIS NOT FOUND:'
         IERR=102
         RETURN
      ENDIF
C
      WRITE(6,*) RAXIS,ZAXIS,PSIXH(RAXIS,ZAXIS)
C
      IF(RAXIS.LE.RR+RB.AND.
     &   RAXIS.GE.RR-RB.AND.
     &   ZAXIS.LE.RKAP*RB.AND.
     &   ZAXIS.GE.-RKAP*RB) THEN
         PSI0=PSIXH(RAXIS,ZAXIS)
         PSIPA=-PSI0
      ELSE
         WRITE(6,'(A)') 'XX EQAXIS: AXIS OUT OF PLASMA:'
         IERR=103
         RETURN
      ENDIF
C
C     ----- calculate outer plasma surface -----
C
      RMAX=MAX(RR+RB,RGMAX)
      REDGE=FBRENT(PSIXZ0,RR-0.5D0*RA,RMAX,1.D-8)
      WRITE(6,*) REDGE,PSIG(REDGE,ZAXIS)
C
      RETURN
      END
C
      FUNCTION PSIXZ0(R)
      INCLUDE '../eq/eqcomx.inc'
      PSIXZ0=PSIXH(R,0.D0)
      RETURN
      END
C
      FUNCTION PSIXH(R,Z)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
      PSIXH=0.D0
C
      IF(RG(NRGMAX).EQ.RG(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(ZG(NZGMAX).EQ.ZG(1)) THEN
         IERR=8
         RETURN
      ENDIF
C
      FRG=1.D0/(RG(NRGMAX)-RG(1))
      NRG=NINT((R-RG(1))*FRG*(NRGMAX-1))+1
      IF(NRG.LT.1) THEN
         IERR=1
         NRG=2
      ENDIF
      IF(NRG.GT.NRGMAX) THEN
         IERR=2
         NRG=NRGMAX
      ENDIF
      FZG=1.D0/(ZG(NZGMAX)-ZG(1))
      NZG=NINT((Z-ZG(1))*FZG*(NZGMAX-1))+1
      IF(NZG.LT.1) THEN
         IERR=3
         NZG=2
      ENDIF
      IF(NZG.GT.NZGMAX) THEN
         IERR=4
         NZG=NZGMAX
      ENDIF
C
 5001 IF(NRG.GE.NRGMAX) GOTO 5002
      IF((R-RG(NRG  ))*FRG.LE.0.D0) GOTO 5002
         NRG=NRG+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NRG.LE.2) GOTO 5004
      IF((R-RG(NRG-1))*FRG.GE.0.D0) GOTO 5004
         NRG=NRG-1
         GOTO 5003
 5004 CONTINUE
      IF(NRG.LT.2)     NRG=2
C
 5005 IF(NZG.GE.NZGMAX) GOTO 5006
      IF((Z-ZG(NZG  ))*FZG.LE.0.D0) GOTO 5006
         NZG=NZG+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NZG.LE.2) GOTO 5008
      IF((Z-ZG(NZG-1))*FZG.GE.0.D0) GOTO 5008
         NZG=NZG-1
         GOTO 5007
 5008 CONTINUE
      IF(NZG.LT.2)     NZG=2
C
      DRG=(RG(NRG)-R)/(RG(NRG)-RG(NRG-1))
      DZG=(ZG(NZG)-Z)/(ZG(NZG)-ZG(NZG-1))
C
      HR0=HRMT0(DRG)
      HZ0=HRMT0(DZG)
      HR1=HRMT1(DRG)
      HZ1=HRMT1(DZG)
C
      HR0C=HRMT0(1.D0-DRG)
      HZ0C=HRMT0(1.D0-DZG)
      HR1C=HRMT1(1.D0-DRG)
      HZ1C=HRMT1(1.D0-DZG)
C
      PSIXH=UPSIX(1,NRG-1,NZG-1)*HR0 *HZ0
     &     +UPSIX(2,NRG-1,NZG-1)*HR0 *HZ1
     &     +UPSIX(3,NRG-1,NZG-1)*HR1 *HZ0
     &     +UPSIX(4,NRG-1,NZG-1)*HR1 *HZ1
     &     +UPSIX(1,NRG  ,NZG-1)*HR0C*HZ0
     &     +UPSIX(2,NRG  ,NZG-1)*HR0C*HZ1
     &     +UPSIX(3,NRG  ,NZG-1)*HR1C*HZ0
     &     +UPSIX(4,NRG  ,NZG-1)*HR1C*HZ1
     &     +UPSIX(1,NRG-1,NZG  )*HR0 *HZ0C
     &     +UPSIX(2,NRG-1,NZG  )*HR0 *HZ1C
     &     +UPSIX(3,NRG-1,NZG  )*HR1 *HZ0C
     &     +UPSIX(4,NRG-1,NZG  )*HR1 *HZ1C
     &     +UPSIX(1,NRG  ,NZG  )*HR0C*HZ0C
     &     +UPSIX(2,NRG  ,NZG  )*HR0C*HZ1C
     &     +UPSIX(3,NRG  ,NZG  )*HR1C*HZ0C
     &     +UPSIX(4,NRG  ,NZG  )*HR1C*HZ1C
C
      RETURN
      END
C
      SUBROUTINE PSIXHD(R,Z,DPSIDR,DPSIDZ)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
      IF(RG(NRGMAX).EQ.RG(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(ZG(NZGMAX).EQ.ZG(1)) THEN
         IERR=8
         RETURN
      ENDIF
C
      FRG=1.D0/(RG(NRGMAX)-RG(1))
      NRG=NINT((R-RG(1))*FRG*(NRGMAX-1))+1
      IF(NRG.LT.1) THEN
         IERR=1
         NRG=2
      ENDIF
      IF(NRG.GT.NRGMAX) THEN
         IERR=2
         NRG=NRGMAX
      ENDIF
      FZG=1.D0/(ZG(NZGMAX)-ZG(1))
      NZG=NINT((Z-ZG(1))*FZG*(NZGMAX-1))+1
      IF(NZG.LT.1) THEN
         IERR=3
         NZG=2
      ENDIF
      IF(NZG.GT.NZGMAX) THEN
         IERR=4
         NZG=NZGMAX
      ENDIF
C
 5001 IF(NRG.GE.NRGMAX) GOTO 5002
      IF((R-RG(NRG  ))*FRG.LE.0.D0) GOTO 5002
         NRG=NRG+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NRG.LE.2) GOTO 5004
      IF((R-RG(NRG-1))*FRG.GE.0.D0) GOTO 5004
         NRG=NRG-1
         GOTO 5003
 5004 CONTINUE
      IF(NRG.LT.2)     NRG=2
C
 5005 IF(NZG.GE.NZGMAX) GOTO 5006
      IF((Z-ZG(NZG  ))*FZG.LE.0.D0) GOTO 5006
         NZG=NZG+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NZG.LE.2) GOTO 5008
      IF((Z-ZG(NZG-1))*FZG.GE.0.D0) GOTO 5008
         NZG=NZG-1
         GOTO 5007
 5008 CONTINUE
      IF(NZG.LT.2)     NZG=2
C
C      WRITE(6,'(A,2I5,1P4E12.4)') 'NRG,NZG,RG,ZG=',NRG,NZG,
C     &     RG(NRG-1),RG(NRG),ZG(NZG-1),ZG(NZG)
C
      DRG=(RG(NRG)-R)/(RG(NRG)-RG(NRG-1))
      DZG=(ZG(NZG)-Z)/(ZG(NZG)-ZG(NZG-1))
C
C      WRITE(6,'(A,1P2E12.4)') 'DRG,DZG=',DRG,DZG
C
      HR0=HRMT0(DRG)
      HZ0=HRMT0(DZG)
      HR1=HRMT1(DRG)
      HZ1=HRMT1(DZG)
C      WRITE(6,'(A,1P4E12.4)') 'HR0...=',HR0,HZ0,HR1,HZ1
C
      HR0C=HRMT0(1.D0-DRG)
      HZ0C=HRMT0(1.D0-DZG)
      HR1C=HRMT1(1.D0-DRG)
      HZ1C=HRMT1(1.D0-DZG)
C      WRITE(6,'(A,1P4E12.4)') 'HR0C..=',HR0C,HZ0C,HR1C,HZ1C
C
      DHR0=-DHRMT0(DRG)/(RG(NRG)-RG(NRG-1))
      DHZ0=-DHRMT0(DZG)/(ZG(NZG)-ZG(NZG-1))
      DHR1=-DHRMT1(DRG)/(RG(NRG)-RG(NRG-1))
      DHZ1=-DHRMT1(DZG)/(ZG(NZG)-ZG(NZG-1))
C      WRITE(6,'(A,1P4E12.4)') 'DHR0..=',DHR0,DHZ0,DHR1,DHZ1
C
      DHR0C=DHRMT0(1.D0-DRG)/(RG(NRG)-RG(NRG-1))
      DHZ0C=DHRMT0(1.D0-DZG)/(ZG(NZG)-ZG(NZG-1))
      DHR1C=DHRMT1(1.D0-DRG)/(RG(NRG)-RG(NRG-1))
      DHZ1C=DHRMT1(1.D0-DZG)/(ZG(NZG)-ZG(NZG-1))
C      WRITE(6,'(A,1P4E12.4)') 'DHR0C.=',DHR0C,DHZ0C,DHR1C,DHZ1C
C
      DPSIDR=UPSIX(1,NRG-1,NZG-1)*DHR0 *HZ0
     &      +UPSIX(2,NRG-1,NZG-1)*DHR0 *HZ1
     &      +UPSIX(3,NRG-1,NZG-1)*DHR1 *HZ0
     &      +UPSIX(4,NRG-1,NZG-1)*DHR1 *HZ1
     &      +UPSIX(1,NRG  ,NZG-1)*DHR0C*HZ0
     &      +UPSIX(2,NRG  ,NZG-1)*DHR0C*HZ1
     &      +UPSIX(3,NRG  ,NZG-1)*DHR1C*HZ0
     &      +UPSIX(4,NRG  ,NZG-1)*DHR1C*HZ1
     &      +UPSIX(1,NRG-1,NZG  )*DHR0 *HZ0C
     &      +UPSIX(2,NRG-1,NZG  )*DHR0 *HZ1C
     &      +UPSIX(3,NRG-1,NZG  )*DHR1 *HZ0C
     &      +UPSIX(4,NRG-1,NZG  )*DHR1 *HZ1C
     &      +UPSIX(1,NRG  ,NZG  )*DHR0C*HZ0C
     &      +UPSIX(2,NRG  ,NZG  )*DHR0C*HZ1C
     &      +UPSIX(3,NRG  ,NZG  )*DHR1C*HZ0C
     &      +UPSIX(4,NRG  ,NZG  )*DHR1C*HZ1C
C
      DPSIDZ=UPSIX(1,NRG-1,NZG-1)*HR0 *DHZ0
     &      +UPSIX(2,NRG-1,NZG-1)*HR0 *DHZ1
     &      +UPSIX(3,NRG-1,NZG-1)*HR1 *DHZ0
     &      +UPSIX(4,NRG-1,NZG-1)*HR1 *DHZ1
     &      +UPSIX(1,NRG  ,NZG-1)*HR0C*DHZ0
     &      +UPSIX(2,NRG  ,NZG-1)*HR0C*DHZ1
     &      +UPSIX(3,NRG  ,NZG-1)*HR1C*DHZ0
     &      +UPSIX(4,NRG  ,NZG-1)*HR1C*DHZ1
     &      +UPSIX(1,NRG-1,NZG  )*HR0 *DHZ0C
     &      +UPSIX(2,NRG-1,NZG  )*HR0 *DHZ1C
     &      +UPSIX(3,NRG-1,NZG  )*HR1 *DHZ0C
     &      +UPSIX(4,NRG-1,NZG  )*HR1 *DHZ1C
     &      +UPSIX(1,NRG  ,NZG  )*HR0C*DHZ0C
     &      +UPSIX(2,NRG  ,NZG  )*HR0C*DHZ1C
     &      +UPSIX(3,NRG  ,NZG  )*HR1C*DHZ0C
     &      +UPSIX(4,NRG  ,NZG  )*HR1C*DHZ1C
C
      RETURN
      END
C
      FUNCTION HRMT0(X)
      IMPLICIT NONE
      REAL*8 HRMT0,X
      HRMT0=X*X*(3.D0-2.D0*X)
      RETURN
      END
C
      FUNCTION HRMT1(X)
      IMPLICIT NONE
      REAL*8 HRMT1,X
      HRMT1=X*X*(X-1.D0)
      RETURN
      END
C
      FUNCTION DHRMT0(X)
      IMPLICIT NONE
      REAL*8 DHRMT0,X
      DHRMT0=6.D0*X*(1.D0-X)
      RETURN
      END
C
      FUNCTION DHRMT1(X)
      IMPLICIT NONE
      REAL*8 DHRMT1,X
      DHRMT1=X*(3.D0*X-2.D0)
      RETURN
      END
