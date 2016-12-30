MODULE xxcomm_parm

  USE bpsd_kinds
  USE bpsd_constants

  INTEGER:: nxmax
  REAL(rkind):: x0,dx,am

END MODULE xxcomm_parm

MODULE xxcomm

  USE xxcomm_parm

  REAL(rkind),ALLOCATABLE,DIMENSION(:):: x,y

CONTAINS

  SUBROUTINE xx_allocate

    INTEGER,SAVE:: nxmax_save = 0

    IF(nxmax == nxmax_save) RETURN

    CALL xx_deallocate

    ALLOCATE(x(nxmax))
    ALLOCATE(y(nxmax))

    nxmax_save=nxmax

    RETURN
  END SUBROUTINE xx_allocate

  SUBROUTINE xx_deallocate

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(y)) DEALLOCATE(y)

    RETURN
  END SUBROUTINE xx_deallocate
END MODULE xxcomm
