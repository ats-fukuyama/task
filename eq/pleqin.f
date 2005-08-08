C     $Id$
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
     &            RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
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
