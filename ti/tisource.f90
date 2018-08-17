! MODULE tisource

MODULE tisource

  PUBLIC ti_source

CONTAINS

  SUBROUTINE ti_source(NR)
    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind):: RHON
    INTEGER:: NSA

    RHON=RM(NR)/RA

    DO NSA=1,NSAMAX
       PNB(NSA,NR)=0.D0
       SNB(NSA,NR)=0.D0
       AJNB(NSA,NR)=0.D0
       PEC(NSA,NR)=0.D0
       AJEC(NSA,NR)=0.D0
       PLH(NSA,NR)=0.D0
       AJLH(NSA,NR)=0.D0
       PIC(NSA,NR)=0.D0
       AJIC(NSA,NR)=0.D0
       PNF(NSA,NR)=0.D0
       SNF(NSA,NR)=0.D0
       SPEL(NSA,NR)=0.D0
       SPSC(NSA,NR)=0.D0
       POH(NSA,NR)=0.D0
       PRB(NSA,NR)=0.D0
       PRC(NSA,NR)=0.D0
       PRL(NSA,NR)=0.D0
       PCX(NSA,NR)=0.D0
       PIE(NSA,NR)=0.D0
       SCX(NSA,NR)=0.D0
       SIE(NSA,NR)=0.D0
    END DO

!    CALL TIPNB(NR)
!    CALL TIPEC(NR)
!    CALL TIPLH(NR)
!    CALL TIPIC(NR)
!    CALL TIPNF(NR)
!    CALL TIPEL(NR)
!    CALL TIPOH(NR)
!    CALL TIPRB(NR)
!    CALL TIPRC(NR)
!    CALL TIPRL(NR)
!    CALL TISCX(NR)
!    CALL TISIE(NR)

    DO NSA=1,NSAMAX
       PSIN(NSA,NR)=PNB(NSA,NR)+PEC(NSA,NR)+PLH(NSA,NR)+PIC(NSA,NR) &
                   +PNF(NSA,NR)+POH(NSA,NR)+PRB(NSA,NR)+PRC(NSA,NR) &
                   +PRL(NSA,NR)+PCX(NSA,NR)+PIE(NSA,NR)
       VSIN(NSA,NR)=0.D0
       SSIN(NSA,NR)=SNB(NSA,NR)+SNF(NSA,NR)+SPEL(NSA,NR)+SPSC(NSA,NR) &
                   +SCX(NSA,NR)+SIE(NSA,NR)
       AJIN(NSA,NR)=AJNB(NSA,NR)+AJEC(NSA,NR)+AJLH(NSA,NR)+AJIC(NSA,NR) &
                   +AJBS(NSA,NR)
    END DO

    RETURN
  END SUBROUTINE ti_source
END MODULE tisource
