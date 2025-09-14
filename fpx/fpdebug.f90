! fpdebug.f90

MODULE fpdebug

  PRIVATE
  PUBLIC fp_debug_drr

CONTAINS
  
  SUBROUTINE fp_debug_drr
    USE fpcomm
    IMPLICIT NONE
    

    WRITE(70,'(3I8)') NTHMAX/2,NPMAX/2,NRMAX/2
    WRITE(70,'(6ES12.4)') &
         DRR(NTHMAX/2-1,NPMAX/2,NRMAX/2,2), &
         DRR(NTHMAX/2,  NPMAX/2,NRMAX/2,2), &
         DRR(NTHMAX/2+1,NPMAX/2,NRMAX/2,2), &
         DRR(NTHMAX/2-1,NPMAX/2-1,NRMAX/2,2), &
         DRR(NTHMAX/2,  NPMAX/2-1,NRMAX/2,2), &
         DRR(NTHMAX/2+1,NPMAX/2-1,NRMAX/2,2)

    RETURN
  END SUBROUTINE fp_debug_drr
END MODULE fpdebug
