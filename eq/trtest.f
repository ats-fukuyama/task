C     $Id$
C
      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      PARAMETER (NTRM=51)
      DIMENSION RHOM(NTRM),RHOG(NTRM)
      DIMENSION PSIP(NTRM),PSIT(NTRM)
      DIMENSION PPSI(NTRM),HJPSI(NTRM)
      DIMENSION VTPSI(NTRM),TPSI(NTRM)
      DIMENSION QPSI(NTRM),VPSI(NTRM),SPSI(NTRM)
      DIMENSION RJ(NTRM),BP(NTRM),QP(NTRM)
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AEE    = 1.60219D-19
C
      CALL GSOPEN
      CALL TREQIN
C
      NTRMAX=51
      RR    = 3.D0
      RA    = 1.D0
      RB    = RA*1.1D0
      RKAP  = 1.6D0
      RDLT  = 0.25D0
      BB    = 3.D0
      RIP   = 3.D0
      PN0   = 0.5D0
      PT0   = 6.D0
      PROFJ1= 2
      PROFJ2= 2
C
      DRHO=1.D0/(NTRMAX-1)
      DO NTR=1,NTRMAX
         RHOG(NTR)=DRHO*(NTR-1)
         RHOM(NTR)=DRHO*(DBLE(NTR)-0.5D0)
         RJ(NTR)= (1.D0-RHOG(NTR)**PROFJ1)**PROFJ2
      ENDDO
C
      BP(1)=0.5D0*DRHO*RJ(1)
      DO NTR=2,NTRMAX
         BP(NTR)=(RHOM(NTR-1)*BP(NTR-1)+RHOG(NTR-1)*DRHO*RJ(NTR-1))
     &           /RHOM(NTR)
      ENDDO
      BPA= RMU0*RIP*1.D6/(2.D0*PI*RA)
      FACT=BPA/BP(NTRMAX)
      DO NTR=1,NTRMAX
         BP(NTR)=FACT*BP(NTR)
         QP(NTR)=RHOM(NTR)*RA*BB/(RR*BP(NTR))
      ENDDO
C
      P0=2*PN0*1.D20*PT0*1.D3*AEE/1.D6
      PSIT(1)=0.D0
      PSIP(1)=0.D0
      PPSI(1)=P0
      HJPSI(1)=RJ(1)
      DO NTR=2,NTRMAX
         PSIT(NTR)=PI*RHOG(NTR)**2*BB
         PSIP(NTR)=(PSIT(NTR)-PSIT(NTR-1))/QP(NTR)+PSIP(NTR-1)
         PPSI(NTR)=P0*(1.D0-RHOG(NTR)**2)
         HJPSI(NTR)=RJ(NTR)
         VTPSI(NTR)=0.D0
         TPSI(NTR)=PT0*(1.D0-RHOG(NTR)**2)
      ENDDO
C
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P5E12.4)') NTR,RHOG(NTR),RHOM(NTR),
     &                            RJ(NTR),BP(NTR),QP(NTR)
      ENDDO
C
C      DO NTR=1,NTRMAX
C         WRITE(6,'(1P4E12.4)') PSIP(NTR),PSIT(NTR),PPSI(NTR),HJPSI(NTR)
C      ENDDO
C
      CALL TREQEX(RR,RA,RB,RKAP,RDLT,BB,RIP,
     &            NTRMAX,PSIP,PPSI,HJPSI,VTPSI,TPSI,
     &            QPSI,VPSI,SPSI,IERR)
C
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P5E12.4)') NTR,PSIP(NTR),QP(NTR),QPSI(NTR),
     &                            VPSI(NTR),SPSI(NTR)
      ENDDO
C
 9000 CALL GSCLOS
      STOP
      END
