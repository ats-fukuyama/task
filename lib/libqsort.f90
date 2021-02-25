MODULE libqsort

! Recursive Fortran 95 quicksort routine
! sorts numbers and associated array into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

! According the order of A, B is sorted.  

  USE task_kinds,ONLY: dp,ikind
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qsort_ic,qsort_lc,qsort_di,qsort_dl,qsort_ii,qsort_li

CONTAINS

RECURSIVE SUBROUTINE qsort_ic(A, B)
  integer, intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call partition_ic(A, iq, B)
     call qsort_ic(A(:iq-1),B(:iq-1))
     call qsort_ic(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_ic

SUBROUTINE partition_ic(A, marker, B)
  integer, intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  complex(dp):: ctemp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ctemp = B(i)
        B(i) = B(j)
        B(j) = ctemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

END SUBROUTINE partition_ic

RECURSIVE SUBROUTINE qsort_lc(A, B)
  INTEGER(ikind), intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  INTEGER(ikind) :: iq

  if(size(A) > 1) then
     call partition_lc(A, iq, B)
     call qsort_lc(A(:iq-1),B(:iq-1))
     call qsort_lc(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_lc

SUBROUTINE PARTITION_LC(A, marker, B)
  INTEGER(ikind), intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  INTEGER(ikind), intent(out) :: marker
  INTEGER(ikind) :: i, j
  INTEGER(ikind) :: temp
  complex(dp):: ctemp
  INTEGER(ikind) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ctemp = B(i)
        B(i) = B(j)
        B(j) = ctemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

END SUBROUTINE partition_lc

RECURSIVE SUBROUTINE qsort_di(A, B)
  real(dp), intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call partition_di(A, iq, B)
     call qsort_di(A(:iq-1),B(:iq-1))
     call qsort_di(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_di

SUBROUTINE partition_di(A, marker, B)
  real(dp), intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  real(dp) :: temp
  integer :: ctemp
  real(dp) :: x      ! pivot point

  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ctemp = B(i)
        B(i) = B(j)
        B(j) = ctemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

END SUBROUTINE partition_di

RECURSIVE SUBROUTINE qsort_dl(A, B)
  real(dp), intent(in out), dimension(:) :: A
  INTEGER(ikind), intent(in out), dimension(:) :: B
  INTEGER(ikind) :: iq

  if(size(A) > 1) then
     call partition_dl(A, iq, B)
     call qsort_dl(A(:iq-1),B(:iq-1))
     call qsort_dl(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_dl

SUBROUTINE partition_dl(A, marker, B)
  real(dp), intent(in out), dimension(:) :: A
  INTEGER(ikind), intent(in out), dimension(:) :: B
  INTEGER(ikind), intent(out) :: marker
  INTEGER(ikind) :: i, j
  real(dp) :: temp
  INTEGER(ikind) :: ctemp
  real(dp) :: x      ! pivot point

  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ctemp = B(i)
        B(i) = B(j)
        B(j) = ctemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

END SUBROUTINE partition_dl

RECURSIVE SUBROUTINE qsort_li(A, B)
  INTEGER(ikind), intent(in out), dimension(:) :: A
  INTEGER, intent(in out), dimension(:) :: B
  INTEGER :: iq

  if(size(A) > 1) then
     call partition_li(A, iq, B)
     call qsort_li(A(:iq-1),B(:iq-1))
     call qsort_li(A(iq:),B(iq:))
  endif
  RETURN
END SUBROUTINE qsort_li

SUBROUTINE partition_li(A, marker, B)
  INTEGER(ikind), intent(in out), dimension(:) :: A
  INTEGER, intent(in out), dimension(:) :: B
  INTEGER, intent(out) :: marker
  INTEGER :: i, j
  INTEGER(ikind) :: atemp
  INTEGER :: btemp
  INTEGER(ikind) :: x      ! pivot point

  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
  RETURN
END SUBROUTINE partition_li

RECURSIVE SUBROUTINE qsort_ii(A, B)
  INTEGER, intent(in out), dimension(:) :: A
  INTEGER, intent(in out), dimension(:) :: B
  INTEGER :: iq

  if(size(A) > 1) then
     call partition_ii(A, iq, B)
     call qsort_ii(A(:iq-1),B(:iq-1))
     call qsort_ii(A(iq:),B(iq:))
  endif
  RETURN
END SUBROUTINE qsort_ii

SUBROUTINE partition_ii(A, marker, B)
  INTEGER(ikind), intent(in out), dimension(:) :: A
  INTEGER, intent(in out), dimension(:) :: B
  INTEGER, intent(out) :: marker
  INTEGER :: i, j
  INTEGER :: atemp
  INTEGER :: btemp
  INTEGER :: x      ! pivot point

  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        atemp = A(i)
        A(i) = A(j)
        A(j) = atemp
        btemp = B(i)
        B(i) = B(j)
        B(j) = btemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
  RETURN
END SUBROUTINE partition_ii

END MODULE libqsort
