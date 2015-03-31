module libqsort

! Recursive Fortran 95 quicksort routine
! sorts numbers and associated array into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

implicit none
public :: qsort_ic,qsort_lc
private :: partition_ic,partition_lc

contains

recursive subroutine qsort_ic(A, B)
  integer, intent(in out), dimension(:) :: A
  complex(8), intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call partition_ic(A, iq, B)
     call qsort_ic(A(:iq-1),B(:iq-1))
     call qsort_ic(A(iq:),B(iq:))
  endif
end subroutine qsort_ic

subroutine partition_ic(A, marker, B)
  integer, intent(in out), dimension(:) :: A
  complex(8), intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  complex(8):: ctemp
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

end subroutine partition_ic

recursive subroutine qsort_lc(A, B)
  integer(8), intent(in out), dimension(:) :: A
  complex(8), intent(in out), dimension(:) :: B
  integer(8) :: iq

  if(size(A) > 1) then
     call partition_lc(A, iq, B)
     call qsort_lc(A(:iq-1),B(:iq-1))
     call qsort_lc(A(iq:),B(iq:))
  endif
end subroutine qsort_lc

subroutine partition_lc(A, marker, B)
  integer(8), intent(in out), dimension(:) :: A
  complex(8), intent(in out), dimension(:) :: B
  integer(8), intent(out) :: marker
  integer(8) :: i, j
  integer(8) :: temp
  complex(8):: ctemp
  integer(8) :: x      ! pivot point
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

end subroutine partition_lc

end module libqsort
