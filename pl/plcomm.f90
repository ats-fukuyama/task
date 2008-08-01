!     $Id$

      module plcomm

      use bpsd_kinds
      use bpsd_constants

      implicit none

      public pl_allocate_ns

      integer:: NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW

      real(rkind):: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
      real(rkind):: PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
      real(rkind):: RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX

      real(rkind),dimension(:),allocatable:: & ! size=NSM
           PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PNITB,PTITB,PUITB

      character(len=80):: KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR

      contains

        subroutine pl_allocate_ns
          implicit none
          integer,save:: NSMAX_save
          integer,save:: init=0
          integer NS

          if(init.eq.0) then
             init=1
          else
             if(NSMAX.eq.NSMAX_save) return
             call pl_deallocate_ns
          endif

          allocate(PA(NSMAX))
          allocate(PZ(NSMAX))
          allocate(PZ0(NSMAX))

          allocate(PN(NSMAX))
          allocate(PNS(NSMAX))
          allocate(PTPR(NSMAX))
          allocate(PTPP(NSMAX))
          allocate(PTS(NSMAX))
          allocate(PU(NSMAX))
          allocate(PUS(NSMAX))
          allocate(PNITB(NSMAX))
          allocate(PTITB(NSMAX))
          allocate(PUITB(NSMAX))

          NSMAX_save=NSMAX

          IF(NSMAX.GE.1) THEN
             NS=1
             PA(NS)   = AME/AMP
             PZ(NS)   =-1.0D0
             PZ0(NS)  =-1.0D0
             PN(NS)   = 1.0D0
             PNS(NS)  = 0.0D0
             PTPR(NS) = 5.0D0
             PTPP(NS) = 5.0D0
             PTS(NS)  = 0.05D0
             PU(NS)   = 0.D0
             PUS(NS)  = 0.D0
             PNITB(NS)= 0.D0
             PTITB(NS)= 0.D0
             PUITB(NS)= 0.D0
          ENDIF

          IF(NSMAX.GE.2) THEN
             NS=2
             PA(NS)   = 1.0D0
             PZ(NS)   = 1.0D0
             PZ0(NS)  = 1.00D0
             PN(NS)   = 1.0D0
             PNS(NS)  = 0.0D0
             PTPR(NS) = 5.0D0
             PTPP(NS) = 5.0D0
             PTS(NS)  = 0.05D0
             PU(NS)   = 0.D0
             PUS(NS)  = 0.D0
             PNITB(NS)= 0.D0
             PTITB(NS)= 0.D0
             PUITB(NS)= 0.D0
          ENDIF

          DO NS=3,NSMAX
             PA(NS)   = 1.0D0
             PZ(NS)   = 1.0D0
             PZ0(NS)  = 0.0D0
             PN(NS)   = 0.0D0
             PNS(NS)  = 0.0D0
             PTPR(NS) = 5.0D0
             PTPP(NS) = 5.0D0
             PTS(NS)  = 0.0D0
             PU(NS)   = 0.D0
             PUS(NS)  = 0.D0
             PNITB(NS)= 0.D0
             PTITB(NS)= 0.D0
             PUITB(NS)= 0.D0
          ENDDO
          return
        end subroutine pl_allocate_ns

        subroutine pl_deallocate_ns
          implicit none

          deallocate(PA)
          deallocate(PZ)
          deallocate(PZ0)

          deallocate(PN)
          deallocate(PNS)
          deallocate(PTPR)
          deallocate(PTPP)
          deallocate(PTS)
          deallocate(PU)
          deallocate(PUS)
          deallocate(PNITB)
          deallocate(PTITB)
          deallocate(PUITB)
          return

        end subroutine pl_deallocate_ns

      end module plcomm
