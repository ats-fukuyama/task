! wrexec.f90

MODULE wrexec

  PRIVATE
  PUBLIC wr_exec

CONTAINS

  SUBROUTINE wr_exec(nstat,ierr)

    USE wrcomm
    USE wrexecr
    USE wrexecb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: nstat,ierr
    INTEGER:: NS,NSS
    REAL(rkind):: ZA1,ZA2

!     --- eliminate disp factor for same z/a species ---

    DO NS=1,NSMAX
       NSDP(NS)=1
       DO NSS=1,NS-1
          ZA1=PZ(NS)/PA(NS)
          ZA2=PZ(NSS)/PA(NSS)
          IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
       ENDDO
    ENDDO

!     --- exec ray/beam tracing ---

    IF(mode_beam.EQ.0) THEN
       CALL wr_exec_rays(ierr)
       nstat=1
    ELSE
       CALL wr_exec_beams(ierr)
       nstat=2
    END IF
  END SUBROUTINE wr_exec
END MODULE wrexec
