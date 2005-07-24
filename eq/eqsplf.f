C
C     ***** INTERPOLATE FUNCTIONS *****
C
      FUNCTION FNPSIN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      FTL=FTSA*RHON*RHON
      CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIN: SPL1DF ERROR : IERR=',IERR
      FNPSIN=1.D0-PSIL/PSI0
      RETURN
      END
C
      FUNCTION FNPSS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      FTL=FTSA*RHON*RHON
      CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIN: SPL1DF ERROR : IERR=',IERR
      FNPSS=PSIL
      RETURN
      END
C
      FUNCTION FNFTS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNFTS: SPL1DF ERROR : IERR=',IERR
      FNFTS=FTL
      RETURN
      END
C
      FUNCTION FNPPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,PPL,RHOT,UPPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPPS: SPL1DF ERROR : IERR=',IERR
      FNPPS=PPL
      RETURN
      END
C
      FUNCTION FNTTS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IF(PSIN.LT.1.D0) THEN
         PSIL=PSI0*(1.D0-PSIN)
         CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
         RHOTL=SQRT(FTL/FTSA)
         CALL SPL1DF(RHOTL,TTL,RHOT,UTTS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX FNTTS: SPL1DF ERROR : IERR=',IERR
      ELSE
         TTL=RR*BB
      ENDIF
      FNTTS=TTL
      RETURN
      END
C
      FUNCTION FNQPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,QPL,RHOT,UQPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNQPS: SPL1DF ERROR : IERR=',IERR
      FNQPS=QPL
      RETURN
      END
C
      FUNCTION FNG5S(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      FNG5S=0.D0
      RETURN
      END
C
      FUNCTION FNVPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,VPL,RHOT,UVPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNVPS: SPL1DF ERROR : IERR=',IERR
      FNVPS=VPL
      RETURN
      END
C
      FUNCTION FNSPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(PHOTL,SPL,RHOT,USPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNSPS: SPL1DF ERROR : IERR=',IERR
      FNSPS=SPL
      RETURN
      END
C
      FUNCTION FNRLEN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,RLENL,RHOT,URLEN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRLEN: SPL1DF ERROR : IERR=',IERR
      FNRLEN=RLENL
      RETURN
      END
C
      FUNCTION FNRRMN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,RRMINL,RHOT,URRMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMN: SPL1DF ERROR : IERR=',IERR
      FNRRMN=RRMINL
      RETURN
      END
C
      FUNCTION FNRRMX(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,RRMAXL,RHOT,URRMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMX: SPL1DF ERROR : IERR=',IERR
      FNRRMX=RRMAXL
      RETURN
      END
C
      FUNCTION FNBBMN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(RHOTL,BBMINL,RHOT,UBBMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMN: SPL1DF ERROR : IERR=',IERR
      FNBBMN=BBMINL
      RETURN
      END
C
      FUNCTION FNBBMX(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=PSI0*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      RHOTL=SQRT(FTL/FTSA)
      CALL SPL1DF(PSIL,BBMAXL,PSS,UBBMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMX: SPL1DF ERROR : IERR=',IERR
      FNBBMX=BBMAXL
      RETURN
      END
