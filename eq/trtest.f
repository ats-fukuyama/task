C     $Id$
C
      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      PARAMETER (NTRM=51)
      DIMENSION RHOTR(NTRM)
      DIMENSION PRHO(NTRM),HJRHO(NTRM)
      DIMENSION VTRHO(NTRM),TRHO(NTRM)
      DIMENSION QRHO(NTRM),TTRHO(NTRM),DVRHO(NTRM),DSRHO(NTRM)
      DIMENSION ARHRRHO(NTRM),AIR2RHO(NTRM),ARH1RHO(NTRM)
      DIMENSION ARH2RHO(NTRM),ABB2RHO(NTRM),AIB2RHO(NTRM)
      DIMENSION ARHBRHO(NTRM),EPSRHO(NTRM),RMJRHO(NTRM)
      DIMENSION RMNRHO(NTRM),RKAPRHO(NTRM)
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
      RKAP  = 1.6D0
      RDLT  = 0.25D0
      BB    = 3.D0
      PN0   = 0.05D0
      PT0   = 0.6D0
      PROFJ1= 2
C      PROFJ2= 2.5D0
      PROFJ2= 4.0D0
C
      RIP1  = 3.D0
C
      DRHO=1.D0/NTRMAX
      DO NTR=1,NTRMAX
         RHOTR(NTR)= DRHO*(DBLE(NTR)-0.5D0)
      ENDDO
C
      CALL TREQIN(RR,RA,RKAP,RDLT,BB,IERR)
C
      DRHO=1.D0/NTRMAX
      P0=2*PN0*1.D20*PT0*1.D3*AEE/1.D6
      DO NTR=1,NTRMAX
         RHOTR(NTR)= DRHO*(DBLE(NTR)-0.5D0)
         PRHO(NTR)=P0*(1.D0-RHOTR(NTR)**2)**1.5D0
         HJRHO(NTR)= 1.D6*(1.D0-RHOTR(NTR)**PROFJ1)**PROFJ2
         VTRHO(NTR)=0.D0
         TRHO(NTR)=PT0*(1.D0-RHOTR(NTR)**2)
      ENDDO
C
      ICONT=0
      CALL TREQEX(NTRMAX,RHOTR,PRHO,HJRHO,VTRHO,TRHO,RIP1,ICONT,
     &            RSA,DPSIPDRHOA,IERR)
C
      CALL EQVIEW
C
      CALL TREQGET(NTRMAX,RHOTR,
     &             QRHO,TTRHO,DVRHO,DSRHO,
     &             ARHRRHO,AIR2RHO,ARH1RHO,ARH2RHO,
     &             ABB2RHO,AIB2RHO,ARHBRHO,
     &             EPSRHO,RMJRHO,RMNRHO,RKAPRHO,IERR)
C
      WRITE(6,*) 'RSA,DPSIPDRHOA=',RSA,DPSIPDRHOA
C
      CALL EQGOUT(1)
C
      WRITE(6,*) 'NTR ','RHO         ','q           ','2*Pi*B*R    ',
     &                  'dV/drho     ','dS/drho     ','<r/R>       '
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P6E12.4)') 
     &         NTR,RHOTR(NTR),QRHO(NTR),TTRHO(NTR),
     &             DVRHO(NTR),DSRHO(NTR),EPSRHO(NTR)
      ENDDO
      WRITE(6,*) 'NTR ','RHO         ','<Vrho^2/R^2>','<R_0^2/R^2> ',
     &                  '<Vrho>      ','<Vrho^2>    ','<B^2/B_0^2> '
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P6E12.4)') 
     &         NTR,RHOTR(NTR),ARHRRHO(NTR),AIR2RHO(NTR),
     &             ARH1RHO(NTR),ARH2RHO(NTR),ABB2RHO(NTR)
      ENDDO
      WRITE(6,*) 'NTR ','RHO         ','<B_0^2/B^2> ','<Vrho^2/B^2>',
     &                  'Rmaj        ','Rmin        ','Rkap        '
      DO NTR=1,NTRMAX
         WRITE(6,'(I5,1P6E12.4)') 
     &         NTR,RHOTR(NTR),AIB2RHO(NTR),ARHBRHO(NTR),
     &             RMJRHO(NTR),RMNRHO(NTR),RKAPRHO(NTR)
      ENDDO
C
      CALL GSCLOS
      STOP
      END
