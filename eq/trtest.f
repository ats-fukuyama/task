C     $Id$
C
      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      PARAMETER (NTRM=51)
      DIMENSION RHOTR(NTRM)
      DIMENSION PRHO(NTRM),HJRHO(NTRM)
      DIMENSION VTRHO(NTRM),TRHO(NTRM)
      DIMENSION QRHO(NTRM),TTRHO(NTRM),DVRHO(NTRM)
      DIMENSION ABRHO(NTRM),ARRHO(NTRM)
C
C      PI     = 2.D0*ASIN(1.D0)
C      RMU0   = 4.D0*PI*1.D-7
      AEE    = 1.60219D-19
C
      CALL GSOPEN
C
      NTRMAX= 50
      RR    = 3.D0
      RA    = 1.D0
      RB    = RA*1.1D0
      RKAP  = 1.6D0
      RDLT  = 0.25D0
      BB    = 3.D0
      RIP   = 3.D0
      PN0   = 0.05D0
      PT0   = 0.6D0
      PROFJ1= 2
      PROFJ2= 2.5D0
C
      DRHO=1.D0/NTRMAX
      DO NTR=1,NTRMAX
         RHOTR(NTR)= DRHO*(DBLE(NTR)-0.5D0)
         HJRHO(NTR)= (1.D0-RHOTR(NTR)**PROFJ1)**PROFJ2
      ENDDO
C
      CALL TREQIN(RR,RA,RB,RKAP,RDLT,BB,RIP,
     &            NTRMAX,RHOTR,HJRHO,QRHO,IERR)
C
C      WRITE(6,'(A,I5,1P3E12.4)') 
C     &     ('NTR,RHOTR,HJRHO,QRHO=',
C     &       NTR,RHOTR(NTR),HJRHO(NTR),QRHO(NTR),
C     &                           NTR=1,NTRMAX)
C
      P0=2*PN0*1.D20*PT0*1.D3*AEE/1.D6
      DO NTR=1,NTRMAX
         PRHO(NTR)=P0*(1.D0-RHOTR(NTR)**2)**1.5D0
         HJRHO(NTR)= (1.D0-RHOTR(NTR)**PROFJ1)**PROFJ2
         VTRHO(NTR)=0.D0
         TRHO(NTR)=PT0*(1.D0-RHOTR(NTR)**2)
      ENDDO
C
      CALL TREQEX(NTRMAX,PRHO,HJRHO,VTRHO,TRHO,
     &            QRHO,TTRHO,DVRHO,ABRHO,ARRHO,IERR)
C
      CALL EQGOUT(1)
C
      WRITE(6,*) 'NTR ','RHO         ','q           ','F=BR        ',
     &                  'dV/drho     ','<Vrho^2/R^2>','<1/R^2>     '
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P6E12.4)') 
     &         NTR,RHOTR(NTR),QRHO(NTR),TTRHO(NTR),
     &             DVRHO(NTR),ABRHO(NTR),ARRHO(NTR)
      ENDDO
C
      CALL GSCLOS
      STOP
      END
