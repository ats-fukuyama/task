C     $Id$
C
      INCLUDE '../eq/eqcomc.inc'
C      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      PARAMETER (NTRM=51)
      DIMENSION RHOTR(NTRM)
      DIMENSION PRHO(NTRM),HJRHO(NTRM)
      DIMENSION VTRHO(NTRM),TRHO(NTRM)
      DIMENSION QRHO(NTRM),TTRHO(NTRM),DVRHO(NTRM),DSRHO(NTRM)
      DIMENSION ARHRRHO(NTRM),AIR2RHO(NTRM),ARH1RHO(NTRM)
      DIMENSION ARH2RHO(NTRM),ABB2RHO(NTRM),AIB2RHO(NTRM)
      DIMENSION ARHBRHO(NTRM),EPSRHO(NTRM),RMJRHO(NTRM)
      DIMENSION RMNRHO(NTRM),RKAPRHO(NTRM)
      CHARACTER KIN*1
C      CHARACTER KNAMEQ*80
C
C      PI     = 2.D0*ASIN(1.D0)
C      RMU0   = 4.D0*PI*1.D-7
       AEE    = 1.60219D-19
C
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      NTRMAX= 50
      RR    = 3.D0
      RA    = 1.D0
      RKAP  = 1.6D0
      RDLT  = 0.25D0
      BB    = 3.D0
      PN0EQ = 0.5D0
      PT0 = 0.8D0
      PROFJ1= 2
C      PROFJ2= 2.5D0
      PROFJ2= 4.0D0
C
C      RIP  = 3.D0
      RIP  = 1.D0
C
      DRHO=1.D0/NTRMAX
      P0=2*PN0EQ*1.D20*PT0*1.D3*AEE/1.D6
      DO NTR=1,NTRMAX
         RHOTR(NTR)= DRHO*(DBLE(NTR)-0.5D0)
         PRHO(NTR)=P0*(1.D0-RHOTR(NTR)**2)**1.5D0
         HJRHO(NTR)= 1.D6*(1.D0-RHOTR(NTR)**PROFJ1)**PROFJ2
         VTRHO(NTR)=0.D0
         TRHO(NTR)=PT0*(1.D0-RHOTR(NTR)**2)
      ENDDO
      ICONT=0
C
 1    WRITE(6,*) '## TRTEST MENU: R/RUN L/LOAD&RUN Q/QUIT'
      READ(5,'(A1)') KIN
      CALL GUCPTL(KIN)
      IF(KIN.EQ.'Q') THEN
         GOTO 9000
      ELSEIF(KIN.EQ.'R'.OR.KIN.EQ.'L') THEN
         GOTO 100
      ELSE
         GOTO 1
      ENDIF
C
 100  CONTINUE
C     *** FOR INITIAL PROFILE ************
      CALL TREQIN(RR,RA,RKAP,RDLT,BB,IERR)
C     ************************************
C
      CALL PARM_READ(RR,RA,RKAP,RDLT,BB,RIP,PRHO,HJRHO,VTRHO,TRHO)
C     *** FOR PREVAILING GEOMETRIC QUANTITIES THROUGH COMMON BLOCKS ***
      CALL TREQIN(RR,RA,RKAP,RDLT,BB,IERR)
C     *****************************************************************
      CALL EQPARM(2,'NPRINT=2',IERR)
      CALL EQPARM(2,'NRGMAX=129',IERR)
      CALL EQPARM(2,'NZGMAX=129',IERR)
      CALL EQPARM(2,'NPSMAX=50',IERR)
      CALL EQPARM(2,'NLPMAX=20',IERR)
C      CALL EQPARM(2,'EPSEQ=1.D-5',IERR)
C
      IF(KIN.EQ.'L') THEN
         KNAMEQ='testeq'
         CALL EQLOAD(9,KNAMEQ,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX TRTEST: EQLOAD ERROR: IERR=',IERR
         CALL EQPARM(2,'NPRINT=2',IERR)
         ICONT=1
         RIP1=0.D0
         CALL EQGOUT(1)
      ENDIF
C
      CALL TREQEX(NTRMAX,RHOTR,PRHO,HJRHO,VTRHO,TRHO,RIP1,ICONT,
     &            RSA,DPSIPDRHOA,IERR)
C
      CALL EQVIEW
C
      IF(IERR.EQ.0) THEN
C
      CALL TREQGET(NTRMAX,RHOTR,
     &             QRHO,TTRHO,DVRHO,DSRHO,
     &             ARHRRHO,AIR2RHO,ARH1RHO,ARH2RHO,
     &             ABB2RHO,AIB2RHO,ARHBRHO,
     &             EPSRHO,RMJRHO,RMNRHO,RKAPRHO,IERR)
C
      WRITE(6,*) 'RSA,DPSIPDRHOA=',RSA,DPSIPDRHOA
C
      ENDIF
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
 9000 CALL GSCLOS
      STOP
      END
C
C     *** READING PARAMETER FILE ***
C
      SUBROUTINE PARM_READ
     &     (RR1,RA1,RKAP1,RDLT1,BB1,RIP1,PRHO1,HJRHO1,VTRHO1,TRHO1)
C
      IMPLICIT NONE
      INTEGER IDOPEN, IST1, IST2, KL, NTR, NTRM
      PARAMETER(NTRM=51)
      REAL*8 PRHO(NTRM),HJRHO(NTRM),VTRHO(NTRM),TRHO(NTRM)
      REAL*8 PRHO1(NTRM),HJRHO1(NTRM),VTRHO1(NTRM),TRHO1(NTRM)
      REAL*8 RR,RA,RKAP,RDLT,BB,RIP
      REAL*8 RR1,RA1,RKAP1,RDLT1,BB1,RIP1
      LOGICAL LEX
      CHARACTER*80 LINE
      NAMELIST /TRDATA/ RR,RA,RKAP,RDLT,BB,RIP,
     &                  PRHO,HJRHO,VTRHO,TRHO
C
      IDOPEN=25
      LINE='trdata'
C
      INQUIRE(FILE=LINE,EXIST=LEX,ERR=9100)
      IF(.NOT.LEX) THEN
         WRITE(6,'(A)') '## NO INPUT FILE EXISTS.'
         RETURN
      ENDIF
      OPEN(IDOPEN,FILE=LINE,IOSTAT=IST1,STATUS='OLD',ERR=9100)
      READ(IDOPEN,TRDATA,IOSTAT=IST2,ERR=9200,END=9900)
      CLOSE(IDOPEN)
C
      RR1=RR
      RA1=RA
      RKAP1=RKAP
      RDLT1=RDLT
      BB1=BB
      RIP1=RIP
      DO NTR=1,NTRM
         PRHO1(NTR)=PRHO(NTR)
         HJRHO1(NTR)=HJRHO(NTR)
         VTRHO1(NTR)=VTRHO(NTR)
         TRHO1(NTR)=TRHO(NTR)
      ENDDO
C
      CALL KTRIM(LINE,KL)
      WRITE(6,'(A,A,A)') 
     &     '## FILE (',LINE(1:KL),') IS ASSIGNED FOR PARM INPUT'
      GOTO 9900
C
 9100 WRITE(6,'(A,I6)') 'XX: FAILED TO OPEN PARM FILE : IOSTAT = ', IST1
      STOP
 9200 WRITE(6,'(A,I6)') 'XX: FAILED TO READ PARM FILE : IOSTAT = ', IST2
      STOP
 9900 CONTINUE
      RETURN
      END
