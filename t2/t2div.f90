!
!
!
MODULE T2DIV

  PRIVATE
  PUBLIC T2_DIV

CONTAINS

  SUBROUTINE T2_DIV

    USE T2NGRA, ONLY: T2_NGRA
    USE T2VGRA, ONLY: T2VGRA_EXECUTE
    USE T2INTG, ONLY: T2INTG_EXECUTE
    IMPLICIT NONE
    REAL(4):: e0time_0,e0time_1

    CALL CPU_TIME(e0time_0)
    CALL T2_NGRA
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- Node graph generated:       cpu=', &
                           e0time_1-e0time_0,' [s]'
    CALL CPU_TIME(e0time_0)
    CALL T2VGRA_EXECUTE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- Variable graph generated:   cpu=', &
                           e0time_1-e0time_0,' [s]'
    CALL CPU_TIME(e0time_0)
    print*,'T2DIV.end'
    STOP

    CALL T2INTG_EXECUTE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- Integral table calculated:  cpu=', &
                           e0time_1-e0time_0,' [s]'
    RETURN

  END SUBROUTINE T2_DIV
  
END MODULE T2DIV
