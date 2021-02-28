MODULE libqsort

! Recursive Fortran 95 quicksort routine
! sorts numbers and associated array into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

! According the order of A, B is sorted.  

  USE task_kinds
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qsort_ic,qsort_lc,qsort_di,qsort_dl,qsort_ii,qsort_li

CONTAINS

RECURSIVE SUBROUTINE qsort_ic(A, B)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  COMPLEX(dp), INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: iq

  IF(SIZE(A) > 1) then
     CALL partition_ic(A, iq, B)
     CALL qsort_ic(A(:iq-1),B(:iq-1))
     CALL qsort_ic(A(iq:),B(iq:))
  ENDIF
END SUBROUTINE qsort_ic

SUBROUTINE partition_ic(A, marker, B)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER :: atemp
  INTEGER :: x      ! pivot point
  COMPLEX(dp), INTENT(IN OUT), DIMENSION(:) :: B
  COMPLEX(dp):: btemp
  INTEGER, INTENT(out) :: marker
  INTEGER :: i, j

  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END do

END SUBROUTINE partition_ic

RECURSIVE SUBROUTINE qsort_lc(A, B)
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: A
  COMPLEX(dp), INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER(long) :: iq

  IF(SIZE(A) > 1) THEN
     CALL partition_lc(A, iq, B)
     CALL qsort_lc(A(:iq-1),B(:iq-1))
     CALL qsort_lc(A(iq:),B(iq:))
  ENDIF
END SUBROUTINE qsort_lc

SUBROUTINE partition_lc(A, marker, B)
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER(long) :: atemp
  INTEGER(long) :: x      ! pivot point
  COMPLEX(dp), INTENT(IN OUT), DIMENSION(:) :: B
  COMPLEX(dp):: btemp
  INTEGER(long), INTENT(OUT) :: marker
  INTEGER(long) :: i, j
  
  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END DO

END SUBROUTINE partition_lc

RECURSIVE SUBROUTINE qsort_di(A, B)
  REAL(dp), INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: iq

  IF(SIZE(A) > 1) then
     CALL partition_di(A, iq, B)
     CALL qsort_di(A(:iq-1),B(:iq-1))
     CALL qsort_di(A(iq:),B(iq:))
  ENDIF
END SUBROUTINE qsort_di

SUBROUTINE partition_di(A, marker, B)
  REAL(dp), INTENT(IN OUT), DIMENSION(:) :: A
  REAL(dp) :: atemp
  REAL(dp) :: x      ! pivot point
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: btemp
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: i, j

  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END DO

END SUBROUTINE partition_di

RECURSIVE SUBROUTINE qsort_dl(A, B)
  REAL(dp), INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER(long) :: iq

  IF(SIZE(A) > 1) THEN
     CALL partition_dl(A, iq, B)
     CALL qsort_dl(A(:iq-1),B(:iq-1))
     CALL qsort_dl(A(iq:),B(iq:))
  ENDIF
END SUBROUTINE qsort_dl

SUBROUTINE partition_dl(A, marker, B)
  REAL(dp), INTENT(IN OUT), DIMENSION(:) :: A
  REAL(dp) :: temp
  REAL(dp) :: x      ! pivot point
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER(long) :: ctemp
  INTEGER(long), INTENT(OUT) :: marker
  INTEGER(long) :: i, j

  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ctemp = B(i)
        B(i) = B(j)
        B(j) = ctemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END DO

END SUBROUTINE partition_dl

RECURSIVE SUBROUTINE qsort_li(A, B)
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: iq

  IF(SIZE(A) > 1) THEN
     CALL partition_li(A, iq, B)
     CALL qsort_li(A(:iq-1),B(:iq-1))
     CALL qsort_li(A(iq:),B(iq:))
  ENDIF
  RETURN
END SUBROUTINE qsort_li

SUBROUTINE partition_li(A, marker, B)
  INTEGER(long), INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER(long) :: atemp
  INTEGER(long) :: x      ! pivot point
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: btemp
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: i, j

  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END DO
  RETURN
END SUBROUTINE partition_li

RECURSIVE SUBROUTINE qsort_ii(A, B)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: iq

  IF(SIZE(A) > 1) THEN
     CALL partition_ii(A, iq, B)
     CALL qsort_ii(A(:iq-1),B(:iq-1))
     CALL qsort_ii(A(iq:),B(iq:))
  ENDIF
  RETURN
END SUBROUTINE qsort_ii

SUBROUTINE partition_ii(A, marker, B)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER :: atemp
  INTEGER :: x      ! pivot point
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER :: btemp
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: i, j

  x = A(1)
  i= 0
  j= SIZE(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     ELSEIF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     ENDIF
  END DO
  RETURN
END SUBROUTINE partition_ii

END MODULE libqsort
