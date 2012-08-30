C     $Id$
C
C     ***** INTERPOLATE FUNCTIONS *****
C
C     FNPSIN: evaluate normalized poloidal flux funcion (PSIN=psip/psipa)
C             corresponding to RHON (i.e. PSIN(RHON))
C
      FUNCTION FNPSIN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSITL=PSITA*RHON*RHON
      CALL SPL1DF(PSITL,PSIPL,PSIT,UPSIP,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIN: SPL1DF ERROR : IERR=',IERR
      FNPSIN=PSIPL/PSIPA
      RETURN
      END
C
C     FNRHON: evaluate rho corresponding to normalized psi, PSIN
C             (i.e. rho(PSIN))
C
      FUNCTION FNRHON(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IF(PSIN.LT.0.D0) THEN
         PSIPL=PSIPA*PSIN
         CALL SPL1DF(PSIPL,PSITL,PSIP,UPSIT,NRMAX,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX FNRHON: SPL1DF ERROR : IERR=',IERR
         IF(PSITL.LE.0.D0) THEN
            FNRHON=0.D0
         ELSE
            FNRHON=SQRT(PSITL/PSITA)
         ENDIF
      ELSE
         FNRHON=SQRT(PSIN)
      ENDIF
      RETURN
      END
C
C     FNPSIPT: evaluate PSIP corresponding to toroidal flux function
C            (PSIT) (i.e. PSIP(PSIT))
C
      FUNCTION FNPSIPT(PSITL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(PSITL,PSIPL,PSIT,UPSIP,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIPT: SPL1DF ERROR : IERR=',IERR
      FNPSIPT=PSIPL
      RETURN
      END
C
C     FNPSIP: evaluate PSI corresponding to RHON (i.e. psi(RHON))
C
      FUNCTION FNPSIP(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSITL=PSITA*RHON*RHON
      CALL SPL1DF(PSITL,PSIPL,PSIT,UPSIP,NRMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FNPSIP: SPL1DF ERROR : IERR=',IERR
         WRITE(6,'(A,1P3E12.4)') 'PSIT:',PSITL,PSIT(1),PSIT(NRMAX)
         WRITE(6,'(A,1P2E12.4)') 'PSITAB:',PSITA,PSITB
      ENDIF
      FNPSIP=PSIPL
      RETURN
      END
C
C     FNPSIP: evaluate PSIP-derivative corresponding to RHON (i.e. dpsi/dRHON)
C
      FUNCTION FNDPSIP(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSITL=PSITA*RHON*RHON
      CALL SPL1DD(PSITL,PSIPL,DPSIPL,PSIT,UPSIP,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIN: SPL1DF ERROR : IERR=',IERR
      FNDPSIP=2.D0*PSITA*RHON*DPSIPL
      RETURN
      END
C
C     FNPSIT: evaluate PSIT corresponding to RHON (i.e. PSIT(RHON))
C
      FUNCTION FNPSIT(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      FNPSIT=PSITA*RHON*RHON
      RETURN
      END
C
C     ************************************************************************
C
C     FNPPS: evaluate PPS corresponding to RHON (i.e. PPS(RHON))
C
      FUNCTION FNPPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),PPL,PSIP,UPPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPPS: SPL1DF ERROR : IERR=',IERR
      FNPPS=PPL
      RETURN
      END
C
C     FNTTS: evaluate TTS corresponding to RHON (i.e. TTS(RHON))
C
      FUNCTION FNTTS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IF(RHON.LT.1.D0) THEN
         CALL SPL1DF(FNPSIP(RHON),TTL,PSIP,UTTS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX FNTTS: SPL1DF ERROR : IERR=',IERR
         FNTTS=TTL
      ELSE
         FNTTS=2.D0*PI*BB*RR
      ENDIF
      RETURN
      END
C
C     FNQPS: evaluate QPS corresponding to RHON (i.e. QPS(RHON))
C
      FUNCTION FNQPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIPL=FNPSIP(RHON)
      CALL SPL1DF(PSIPL,QPL,PSIP,UQPS,NRMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FNQPS: SPL1DF ERROR : IERR=',IERR
         WRITE(6,'(A,1PE12.4)') 'RHON:',RHON
         WRITE(6,'(A,1P3E12.4)') 'PSIP:',PSIPL,PSIP(1),PSIP(NRMAX)
      ENDIF
      FNQPS=QPL
      RETURN
      END
C
C     FNVPS: evaluate VPS corresponding to RHON (i.e. VPS(RHON))
C
      FUNCTION FNVPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),VPL,PSIP,UVPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNVPS: SPL1DF ERROR : IERR=',IERR
      FNVPS=VPL
      RETURN
      END
C
C     FNSPS: evaluate SPS corresponding to RHON (i.e. SPS(RHON))
C
      FUNCTION FNSPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),SPL,PSIP,USPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNSPS: SPL1DF ERROR : IERR=',IERR
      FNSPS=SPL
      RETURN
      END
C
C     FNRLEN: evaluate RLEN corresponding to RHON (i.e. RLEN(RHON))
C
      FUNCTION FNRLEN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),RLENL,PSIP,URLEN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRLEN: SPL1DF ERROR : IERR=',IERR
      FNRLEN=RLENL
      RETURN
      END
C
C     FNRRMN: evaluate RRMIN corresponding to RHON (i.e. RRMIN(RHON))
C
      FUNCTION FNRRMN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),RRMINL,PSIP,URRMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMN: SPL1DF ERROR : IERR=',IERR
      FNRRMN=RRMINL
      RETURN
      END
C
C     FNRRMX: evaluate RRMAX corresponding to RHON (i.e. RRMAX(RHON))
C
      FUNCTION FNRRMX(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),RRMAXL,PSIP,URRMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMX: SPL1DF ERROR : IERR=',IERR
      FNRRMX=RRMAXL
      RETURN
      END
C
C     FNZZMN: evaluate ZZMIN corresponding to RHON (i.e. ZZMIN(RHON))
C
      FUNCTION FNZZMN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),ZZMINL,PSIP,UZZMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMN: SPL1DF ERROR : IERR=',IERR
      FNZZMN=ZZMINL
      RETURN
      END
C
C     FNZZMX: evaluate ZZMAX corresponding to RHON (i.e. ZZMAX(RHON))
C
      FUNCTION FNZZMX(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),ZZMAXL,PSIP,UZZMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMX: SPL1DF ERROR : IERR=',IERR
      FNZZMX=ZZMAXL
      RETURN
      END
C
C     FNBBMN: evaluate BBMIN corresponding to RHON (i.e. BBMIN(RHON))
C
      FUNCTION FNBBMN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),BBMINL,PSIP,UBBMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMN: SPL1DF ERROR : IERR=',IERR
      FNBBMN=BBMINL
      RETURN
      END
C
C     FNBBMX: evaluate BBMAX corresponding to RHON (i.e. BBMAX(RHON))
C
      FUNCTION FNBBMX(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),BBMAXL,PSIP,UBBMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMX: SPL1DF ERROR : IERR=',IERR
      FNBBMX=BBMAXL
      RETURN
      END
C
C     FNAVRR2: evaluate AVERR2 corresponding to RHON (i.e. AVERR2(RHON))
C
      FUNCTION FNAVRR2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVERR2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVRR2: SPL1DF ERROR : IERR=',IERR
      FNAVRR2=DAT
      RETURN
      END
C
C     FNAVIR2: evaluate AVEIR2 corresponding to RHON (i.e. AVEIR2(RHON))
C
      FUNCTION FNAVIR2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEIR2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVIR2: SPL1DF ERROR : IERR=',IERR
      FNAVIR2=DAT
      RETURN
      END
C
C     FNAVBB2: evaluate AVEBB2 corresponding to RHON (i.e. AVEBB2(RHON))
C
      FUNCTION FNAVBB2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEBB2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVBB2: SPL1DF ERROR : IERR=',IERR
      FNAVBB2=DAT
      RETURN
      END
C
C     FNAVIB2: evaluate AVEIB2 corresponding to RHON (i.e. AVEIB2(RHON))
C
      FUNCTION FNAVIB2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEIB2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVIB2: SPL1DF ERROR : IERR=',IERR
      FNAVIB2=DAT
      RETURN
      END
C
C     FNAVBB: evaluate AVEBB corresponding to RHON (i.e. AVEBB(RHON))
C
      FUNCTION FNAVBB(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEBB,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVBB: SPL1DF ERROR : IERR=',IERR
      FNAVBB=DAT
      RETURN
      END
C
C     FNAVGV: evaluate AVEGV corresponding to RHON (i.e. AVEGV(RHON))
C
      FUNCTION FNAVGV(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGV,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGV: SPL1DF ERROR : IERR=',IERR
      FNAVGV=DAT
      RETURN
      END
C
C     FNAVGV2: evaluate AVEGV2 corresponding to RHON (i.e. AVEGV2(RHON))
C
      FUNCTION FNAVGV2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGV2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGV2: SPL1DF ERROR : IERR=',IERR
      FNAVGV2=DAT
      RETURN
      END
C
C     FNAVGR2: evaluate AVEGVR2 corresponding to RHON (i.e. AVEGVR2(RHON))
C
      FUNCTION FNAVGVR2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGVR2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGVR2: SPL1DF ERROR : IERR=',IERR
      FNAVGVR2=DAT
      RETURN
      END
C
C     FNAVGP2: evaluate AVEGP2 corresponding to RHON (i.e. AVEGP2(RHON))
C
      FUNCTION FNAVGP2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGP2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGP2: SPL1DF ERROR : IERR=',IERR
      FNAVGP2=DAT
      RETURN
      END
C
C     FNRRPS: evaluate RRPSI corresponding to RHON (i.e. RRPSI(RHON))
C
      FUNCTION FNRRPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,URRPSI,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRPS: SPL1DF ERROR : IERR=',IERR
      FNRRPS=DAT
      RETURN
      END
C
C     FNRSPS: evaluate RSPSI corresponding to RHON (i.e. RSPSI(RHON))
C
      FUNCTION FNRSPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,URSPSI,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRSPS: SPL1DF ERROR : IERR=',IERR
      FNRSPS=DAT
      RETURN
      END
C
C     FNELPPS: evaluate ELIPPSI corresponding to RHON (i.e. ELIPPSI(RHON))
C
      FUNCTION FNELPPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UELIPPSI,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNELPPS: SPL1DF ERROR : IERR=',IERR
      FNELPPS=DAT
      RETURN
      END
C
C     FNTRGPS: evaluate TRIGPSI corresponding to RHON (i.e. TRIGPSI(RHON))
C
      FUNCTION FNTRGPS(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UTRIGPSI,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNTRGPS: SPL1DF ERROR : IERR=',IERR
      FNTRGPS=DAT
      RETURN
      END
C
C     FNDVDPSP: evaluate DVDPSIP corresponding to RHON (i.e. DVDPSIP(RHON))
C
      FUNCTION FNDVDPSP(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UDVDPSIP,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNDVDPSP: SPL1DF ERROR : IERR=',IERR
      FNDVDPSP=DAT
      RETURN
      END
C
C     FNDVDPST: evaluate DVDPSIT corresponding to RHON (i.e. DVDPSIT(RHON))
C
      FUNCTION FNDVDPST(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UDVDPSIT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNDVDPST: SPL1DF ERROR : IERR=',IERR
      FNDVDPST=DAT
      RETURN
      END
C
C     FNAVGR: evaluate AVEGR corresponding to RHON (i.e. AVEGR(RHON))
C
      FUNCTION FNAVGR(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGR,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGR: SPL1DF ERROR : IERR=',IERR
      FNAVGR=DAT
      RETURN
      END
C
C     FNAVGR2: evaluate AVEGR2 corresponding to RHON (i.e. AVEGR2(RHON))
C
      FUNCTION FNAVGR2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGR2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGR2: SPL1DF ERROR : IERR=',IERR
      FNAVGR2=DAT
      RETURN
      END
C
C     FNAVGRR2: evaluate AVEGVRR2 corresponding to RHON (i.e. AVEGVRR2(RHON))
C
      FUNCTION FNAVGRR2(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEGRR2,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVGRR2: SPL1DF ERROR : IERR=',IERR
      FNAVGRR2=DAT
      RETURN
      END
C
C     FNAVIR2: evaluate AVEIR corresponding to RHON (i.e. AVEIR(RHON))
C
      FUNCTION FNAVIR(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(FNPSIP(RHON),DAT,PSIP,UAVEIR,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNAVIR: SPL1DF ERROR : IERR=',IERR
      FNAVIR=DAT
      RETURN
      END
