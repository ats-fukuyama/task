C     $Id$
C
      PARAMETER (NRM=201,NTHM=128)
      PARAMETER (NSUM=1025)
C      PARAMETER (NRM=50,NTHM=32)
C      PARAMETER (NSUM=41)
      PARAMETER (NNM=200)
C
      PARAMETER (NTHMP=NTHM+1)
C
      COMMON /EQWMU1/ RPS(NTHM,NRM),ZPS(NTHM,NRM)
      COMMON /EQWMU2/ DRPSI(NTHM,NRM),DZPSI(NTHM,NRM)
      COMMON /EQWMU3/ DRCHI(NTHM,NRM),DZCHI(NTHM,NRM)
      COMMON /EQWMU4/ PSS(NRM),PPS(NRM),TTS(NRM),FTS(NRM),FTSA,FTSB
      COMMON /EQWMU5/ QPS(NRM),VPS(NRM),SPS(NRM),RLEN(NRM)
      COMMON /EQWMU6/ RRMIN(NRM),RRMAX(NRM),BBMIN(NRM),BBMAX(NRM)
      COMMON /EQWMU7/ RGMIN,RGMAX,ZGMIN,ZGMAX
      COMMON /EQWMU8/ RSU(NSUM),ZSU(NSUM)
      COMMON /EQWMU9/ RSW(NSUM),ZSW(NSUM)
      COMMON /EQTRU1/ AVBR(NRM),AVRR(NRM)
      COMMON /EQTRU2/ AVR1(NRM),AVR2(NRM)
C
      COMMON /EQRSP2/ UPPPS(4,NPSM),UTTPS(4,NPSM)
      COMMON /EQRSP3/ UPPS(4,NRM),UTTS(4,NRM),UFTS(4,NRM),UFTT(4,NRM)
      COMMON /EQRSP4/ UQPS(4,NRM),UVPS(4,NRM),USPS(4,NRM),URLEN(4,NRM)
      COMMON /EQRSP5/ URRMIN(4,NRM),URRMAX(4,NRM)
      COMMON /EQRSP6/ UBBMIN(4,NRM),UBBMAX(4,NRM)
      COMMON /EQRSP7/ UAVBR(4,NRM),UAVRR(4,NRM)
      COMMON /EQRSP8/ UAVR1(4,NRM),UAVR2(4,NRM)
