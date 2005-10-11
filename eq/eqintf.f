C     $Id$
C
C     ****** PARAMETER EXCHANGE INTERFACE ******
C
      SUBROUTINE EQ_SET_PLPARM
C
      INCLUDE '../eq/eqcomm.inc'
      INCLUDE '../pl/plcom0.inc'
      DIMENSION DPARM(22),IPARM(4),DPARMS(13,NSM)
      CHARACTER KPARM(6)*80
C
      CALL PL_GET_PARM(DPARM,IPARM,KPARM,DPARMS)
C
      DPARM( 1) = RR
      DPARM( 2) = RA
      DPARM( 3) = RB
      DPARM( 4) = RKAP
      DPARM( 5) = RDLT
C
      DPARM( 6) = BB
      DPARM( 9) = RIP
C
      KPARM( 1) = KNAMEQ
C
      CALL PL_SET_PARM(DPARM,IPARM,KPARM,DPARMS)
      RETURN
      END
C
C     ****** PARAMETER EXCHANGE INTERFACE ******
C
      SUBROUTINE EQ_GET_PLPARM
C
      INCLUDE '../eq/eqcomm.inc'
      INCLUDE '../pl/plcom0.inc'
      DIMENSION DPARM(22),IPARM(4),DPARMS(13,NSM)
      CHARACTER KPARM(6)*80
C
      CALL PL_GET_PARM(DPARM,IPARM,KPARM,DPARMS)
C
      RR     = DPARM( 1)
      RA     = DPARM( 2)
      RB     = DPARM( 3)
      RKAP   = DPARM( 4)
      RDLT   = DPARM( 5)
C
      BB     = DPARM( 6)
      RIP    = DPARM( 9)
C
      KNAMEQ = KPARM( 1)
C
      RETURN
      END
C
C     ***** GET PSIN AND MAGNETIC FIELD *****
C
C           PSIN=0 ON MAGNETIC AXIS
C           PSIN=1 ON PLASMA BOUNDARY
C
C           PSI=PSI0 ON MAGNETIC AXIS
C           PSI=0     ON PLASMA BOUNDARY
C
      SUBROUTINE GETRZ(RP,ZP,PHIP,BR,BZ,BT,RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL2DD(RP,ZP,PSI,PSIR,PSIZ,
     &            RG,ZG,UPSIRZ,NRGM,NRGMAX,NZGMAX,IERR)
C
      PSIN=1.D0-PSI/PSI0
      IF(PSIN.LE.0.D0) THEN
         PSIN=0.D0
         BT=FNTTS(0.D0)/RR
         BR=0.D0
         BZ=0.D0
      ELSE
         BT=FNTTS(PSIN)/RP
         BR=-PSIZ/RP
         BZ= PSIR/RP
      ENDIF
      RHON=FNRHON(PSIN)
C
      RETURN
      END
C
C     ***** GET PARAMETERS *****
C
      SUBROUTINE EQGETB(BB1,RR1,RIP1,RA1,RKAP1,RDEL1,RB1)
C
      INCLUDE '../eq/eqcomq.inc'
C
      BB1  =BB
      RR1  =RR
      RIP1 =RIP
      RA1  =RA
      RKAP1=RKAP
      RDEL1=RDLT
      RB1  =RB
      RETURN
      END
C
C     ***** GET PRESSURE AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETPP(PSIN,PP)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PP=FNPPS(PSIN)
      RETURN
      END
C
C     ***** GET SAFETY FACTOR AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETQP(PSIN,QP)
C
      INCLUDE '../eq/eqcomq.inc'
C
      QP=FNQPS(PSIN)
      RETURN
      END
C
C     ***** GET MINIMUM R AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETRMN(PSIN,RRMINL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RRMINL=FNRRMN(PSIN)
      RETURN
      END
C
C     ***** GET MAXIMUM R AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETRMX(PSIN,RRMAXL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RRMAXL=FNRRMX(PSIN)
      RETURN
      END
C
C     ***** GET MAGNETIC AXIS *****
C
      SUBROUTINE GETAXS(RAXIS1,ZAXIS1)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RAXIS1=RAXIS
      ZAXIS1=ZAXIS
      RETURN
      END
C
C     ***** GET PLASMA BOUNDARY POSITION *****
C
      SUBROUTINE GETRSU(RSU1,ZSU1,N,NSUMAX1)
C
      INCLUDE '../eq/eqcomq.inc'
      DIMENSION RSU1(N),ZSU1(N)
C
      NSUMAX1=NSUMAX
      DO NSU=1,MIN(N,NSUMAX)
         RSU1(NSU)=RSU(NSU)
         ZSU1(NSU)=ZSU(NSU)
      ENDDO
      RETURN
      END
