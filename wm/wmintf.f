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
      DPARM( 7) = Q0
      DPARM( 8) = QA
      DPARM( 9) = RIP
      DPARM(10) = PROFJ
C
      DPARM(11) = PROFN1
      DPARM(12) = PROFN2
      DPARM(13) = PROFT1
      DPARM(14) = PROFT2
      DPARM(15) = PROFU1
      DPARM(16) = PROFU2
C
      DPARM(17) = RHOMIN
      DPARM(18) = QMIN
      DPARM(19) = RHOITB
      DPARM(20) = RHOEDG
      DPARM(21) = RHOGMN
      DPARM(22) = RHOGMX
C
      IPARM( 1) = NSMAX
      IPARM( 2) = MODELG
      IPARM( 3) = MODELN
      IPARM( 4) = MODELQ
C
      KPARM( 1) = KNAMEQ
      KPARM( 2) = KNAMWR
      KPARM( 3) = KNAMWM
      KPARM( 4) = KNAMFP
      KPARM( 5) = KNAMFO
      KPARM( 6) = KNAMPF
C
      DO NS=1,NSMAX
         DPARMS( 1,NS) = PA(NS)
         DPARMS( 2,NS) = PZ(NS)
         DPARMS( 3,NS) = PN(NS)
         DPARMS( 4,NS) = PNS(NS)
         DPARMS( 5,NS) = PZCL(NS)
         DPARMS( 6,NS) = PTPR(NS)
         DPARMS( 7,NS) = PTPP(NS)
         DPARMS( 8,NS) = PTS(NS)
         DPARMS( 9,NS) = PU(NS)
         DPARMS(10,NS) = PUS(NS)
         DPARMS(11,NS) = PNITB(NS)
         DPARMS(12,NS) = PTITB(NS)
         DPARMS(13,NS) = PUITB(NS)
      ENDDO
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
      Q0     = DPARM( 7)
      QA     = DPARM( 8)
      RIP    = DPARM( 9)
      PROFJ  = DPARM(10)
C
      PROFN1 = DPARM(11)
      PROFN2 = DPARM(12)
      PROFT1 = DPARM(13)
      PROFT2 = DPARM(14)
      PROFU1 = DPARM(15)
      PROFU2 = DPARM(16)
C
      RHOMIN = DPARM(17)
      QMIN   = DPARM(18)
      RHOITB = DPARM(19)
      RHOEDG = DPARM(20)
      RHOGMN = DPARM(21)
      RHOGMX = DPARM(22)
C
      NSMAX  = IPARM( 1)
      MODELG = IPARM( 2)
      MODELN = IPARM( 3)
      MODELQ = IPARM( 4)
C
      KNAMEQ = KPARM( 1)
      KNAMWR = KPARM( 2)
      KNAMWM = KPARM( 3)
      KNAMFP = KPARM( 4)
      KNAMFO = KPARM( 5)
      KNAMPF = KPARM( 6)
C
      DO NS=1,NSMAX
         PA(NS)    = DPARMS( 1,NS)
         PZ(NS)    = DPARMS( 2,NS)
         PN(NS)    = DPARMS( 3,NS)
         PNS(NS)   = DPARMS( 4,NS)
         PZCL(NS)  = DPARMS( 5,NS)
         PTPR(NS)  = DPARMS( 6,NS)
         PTPP(NS)  = DPARMS( 7,NS)
         PTS(NS)   = DPARMS( 8,NS)
         PU(NS)    = DPARMS( 9,NS)
         PUS(NS)   = DPARMS(10,NS)
         PNITB(NS) = DPARMS(11,NS)
         PTITB(NS) = DPARMS(12,NS)
         PUITB(NS) = DPARMS(13,NS)
      ENDDO
C
      RETURN
      END
