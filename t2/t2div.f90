!
!
!
MODULE T2DIV

  PRIVATE
  PUBLIC T2_DIV

CONTAINS

  SUBROUTINE T2_DIV

    USE T2NGRA, ONLY: T2_NGRA
    USE T2VGRA, ONLY: T2_VGRA
    USE T2INTG, ONLY: T2_INTG
    IMPLICIT NONE
  
    CALL T2_NGRA
    WRITE(6,*)'-- Node graph generated'
    CALL T2_VGRA
    WRITE(6,*)'-- Variable graph generated'
    CALL T2_INTG
    WRITE(6,*)'-- Integral table calculated'
    RETURN

  END SUBROUTINE T2_DIV
  
END MODULE T2DIV
