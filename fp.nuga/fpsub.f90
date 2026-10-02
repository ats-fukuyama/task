MODULE fpsub

  PRIVATE
  PUBLIC fpbave_dpp
  PUBLIC fpbave_dth

CONTAINS

  SUBROUTINE fpbave_dpp(DPP,NR,NSA,ID)
    USE fpcomm,ONLY: NTHMAX,NPSTART,NPENDWG,NRSTART,NRENDWM,NSAMAX, &
                     RLAMDA,ITL,ITU
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: &
         DPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX)
    INTEGER,INTENT(IN):: NR,NSA,ID  ! ID=0 for DPP,FPP, ID=1 for DPT
    INTEGER:: NP,NTH

    DO NP=NPSTART,NPENDWG
       IF(ID.EQ.0) THEN
          DO NTH=ITL(NR)+1,NTHMAX/2
             DPP(NTH,NP,NR,NSA)    =0.5D0*(DPP(NTH,NP,NR,NSA) &
                                          +DPP(NTHMAX-NTH+1,NP,NR,NSA))
             DPP(NTHMAX-NTH+1,NP,NR,NSA)  =DPP(NTH,NP,NR,NSA)
          END DO
       ELSE
          DO NTH=ITL(NR)+1,NTHMAX/2
             DPP(NTH,NP,NR,NSA)    =0.5D0*(DPP(NTH,NP,NR,NSA) &
                                          -DPP(NTHMAX-NTH+1,NP,NR,NSA))
             DPP(NTHMAX-NTH+1,NP,NR,NSA) =-DPP(NTH,NP,NR,NSA)
          END DO
       END IF
       IF(ID.EQ.0) THEN
          DPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                      *( DPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                        +DPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                        +DPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                        +DPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR)) 
          DPP(ITU(NR),NP,NR,NSA)=DPP(ITL(NR),NP,NR,NSA)
       ELSE
          DPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                      *( DPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                        +DPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                        -DPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                        -DPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR)) 
          DPP(ITU(NR),NP,NR,NSA)=-DPP(ITL(NR),NP,NR,NSA)
       END IF
    END DO
    RETURN
  END SUBROUTINE fpbave_dpp

  SUBROUTINE fpbave_dth(DTH,NR,NSA,ID)
    USE fpcomm,ONLY: NTHMAX,NPSTARTW,NPENDWM,NRSTART,NRENDWM,NSAMAX, &
                     RLAMDA,ITL,ITU
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: &
         DTH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX)
    INTEGER,INTENT(IN):: NR,NSA,ID  ! ID=0 for DTT, ID=1 for DTP,FTH
    INTEGER:: NP,NTH

    DO NP=NPSTARTW,NPENDWM
       DO NTH=ITL(NR)+1,NTHMAX/2
          IF(ID.EQ.0) THEN
             DTH(NTH,NP,NR,NSA)    =0.5D0*(DTH(NTH,         NP,NR,NSA) &
                                          +DTH(NTHMAX-NTH+2,NP,NR,NSA))
             DTH(NTHMAX-NTH+2,NP,NR,NSA)  =DTH(NTH,NP,NR,NSA)
          ELSE
             DTH(NTH,NP,NR,NSA)    =0.5D0*(DTH(NTH,         NP,NR,NSA) &
                                          -DTH(NTHMAX-NTH+2,NP,NR,NSA))
             DTH(NTHMAX-NTH+2,NP,NR,NSA) =-DTH(NTH,NP,NR,NSA)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE fpbave_dth
END MODULE fpsub
