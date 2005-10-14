C     $Id$
C
C     ****** PARAMETER EXCHANGE INTERFACE ******
C
      SUBROUTINE WM_SET_PLPARM
C
      INCLUDE '../wm/wmcomm.inc'
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
      SUBROUTINE WM_GET_PLPARM
C
      INCLUDE '../wm/wmcomm.inc'
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
