MODULE libqsort

! Recursive Fortran 95 quicksort routine
! sorts numbers and associated array into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

  USE task_kinds,ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qsort_ic,qsort_lc,qsort_di,qsort_dl

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
  integer(dp), intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  integer(dp) :: iq

  if(size(A) > 1) then
     call partition_lc(A, iq, B)
     call qsort_lc(A(:iq-1),B(:iq-1))
     call qsort_lc(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_lc

SUBROUTINE PARTITION_LC(A, marker, B)
  integer(dp), intent(in out), dimension(:) :: A
  complex(dp), intent(in out), dimension(:) :: B
  integer(dp), intent(out) :: marker
  integer(dp) :: i, j
  integer(dp) :: temp
  complex(dp):: ctemp
  integer(dp) :: x      ! pivot point
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
  integer(dp), intent(in out), dimension(:) :: B
  integer(dp) :: iq

  if(size(A) > 1) then
     call partition_dl(A, iq, B)
     call qsort_dl(A(:iq-1),B(:iq-1))
     call qsort_dl(A(iq:),B(iq:))
  endif
END SUBROUTINE qsort_dl

SUBROUTINE partition_dl(A, marker, B)
  real(dp), intent(in out), dimension(:) :: A
  integer(dp), intent(in out), dimension(:) :: B
  integer(dp), intent(out) :: marker
  integer(dp) :: i, j
  real(dp) :: temp
  integer(dp) :: ctemp
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

END MODULE libqsort
