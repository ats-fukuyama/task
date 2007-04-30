!     $Id$
module libfem_mod
  implicit none
  private
  real(8), dimension(1:8,1:8), save :: table_hh
  real(8), dimension(1:8,1:4), save :: table_lh
  real(8), dimension(1:8,1:8,1:8), save :: table_hhh
  real(8), dimension(1:8,1:4,1:8), save :: table_hlh
  real(8), dimension(1:4,1:8,1:8), save :: table_lhh
  real(8), dimension(1:4,1:4,1:8), save :: table_llh
  integer :: table_initialize_flag = 0
  public :: fem_integrate

contains

  subroutine table_initialize

!  Integral value table
!  table_hh : table of double product integrand
!  table_hhh : table of triple product integrand

!  Relation
!     form  : L = h1 + g1 + h2 + g2,  dL = dh1 + dg1 + dh2 + dg2
!     label :     1    2    3    4          5     6     7     8
!  where   h1 = h_i, g1 =-g_i / h, h2 = h_{i+1}, g2 = g_{i+1} / h
!          dh1 = h'_i * h, dg1 = g'_i, dh2 =-h'_{i+1} * h, dg2 = g'_{i+1}
!          h = r_{i+1} - r_i

    implicit none
    integer:: k,i,j

    table_hh(1:8,1:8)         = 0.d0
    table_hhh(1:8,1:8,1:8)     = 0.d0

    !  *** double product hh ***

    table_hh(1,1) = 13.d0 / 35.d0
    table_hh(1,2) = 11.d0 / 210.d0
    table_hh(1,3) =  9.d0 / 70.d0
    table_hh(1,4) =-13.d0 / 420.d0
    table_hh(1,5) =- 1.d0 / 2.d0
    table_hh(1,6) =  1.d0 / 10.d0
    table_hh(1,7) =  1.d0 / 2.d0
    table_hh(1,8) =- 1.d0 / 10.d0

    table_hh(2,1) = 11.d0 / 210.d0
    table_hh(2,2) =  1.d0 / 105.d0
    table_hh(2,3) = 13.d0 / 420.d0
    table_hh(2,4) =- 1.d0 / 140.d0
    table_hh(2,5) =- 1.d0 / 10.d0
    table_hh(2,6) =  0.d0
    table_hh(2,7) =  1.d0 / 10.d0
    table_hh(2,8) =- 1.d0 / 60.d0

    table_hh(3,1) =  9.d0 / 70.d0
    table_hh(3,2) = 13.d0 / 420.d0
    table_hh(3,3) = 13.d0 / 35.d0
    table_hh(3,4) =-11.d0 / 210.d0
    table_hh(3,5) =- 1.d0 / 2.d0
    table_hh(3,6) =- 1.d0 / 10.d0
    table_hh(3,7) =  1.d0 / 2.d0
    table_hh(3,8) =  1.d0 / 10.d0

    table_hh(4,1) =-13.d0 / 420.d0
    table_hh(4,2) =- 1.d0 / 140.d0
    table_hh(4,3) =-11.d0 / 210.d0
    table_hh(4,4) =  1.d0 / 105.d0
    table_hh(4,5) =  1.d0 / 10.d0
    table_hh(4,6) =  1.d0 / 60.d0
    table_hh(4,7) =- 1.d0 / 10.d0
    table_hh(4,8) =  0.d0

    table_hh(5,1) =- 1.d0 / 2.d0
    table_hh(5,2) =- 1.d0 / 10.d0
    table_hh(5,3) =- 1.d0 / 2.d0
    table_hh(5,4) =  1.d0 / 10.d0
    table_hh(5,5) =  6.d0 / 5.d0
    table_hh(5,6) =  1.d0 / 10.d0
    table_hh(5,7) =- 6.d0 / 5.d0
    table_hh(5,8) =  1.d0 / 10.d0

    table_hh(6,1) =  1.d0 / 10.d0
    table_hh(6,2) =  0.d0
    table_hh(6,3) =- 1.d0 / 10.d0
    table_hh(6,4) =  1.d0 / 60.d0
    table_hh(6,5) =  1.d0 / 10.d0
    table_hh(6,6) =  2.d0 / 15.d0
    table_hh(6,7) =- 1.d0 / 10.d0
    table_hh(6,8) =- 1.d0 / 30.d0

    table_hh(7,1) =  1.d0 / 2.d0
    table_hh(7,2) =  1.d0 / 10.d0
    table_hh(7,3) =  1.d0 / 2.d0
    table_hh(7,4) =- 1.d0 / 10.d0
    table_hh(7,5) =- 6.d0 / 5.d0
    table_hh(7,6) =- 1.d0 / 10.d0
    table_hh(7,7) =  6.d0 / 5.d0
    table_hh(7,8) =- 1.d0 / 10.d0

    table_hh(8,1) =- 1.d0 / 10.d0
    table_hh(8,2) =- 1.d0 / 60.d0
    table_hh(8,3) =  1.d0 / 10.d0
    table_hh(8,4) =  0.d0
    table_hh(8,5) =  1.d0 / 10.d0
    table_hh(8,6) =- 1.d0 / 30.d0
    table_hh(8,7) =- 1.d0 / 10.d0
    table_hh(8,8) =  2.d0 / 15.d0

    !  *** double product hl ***

    table_lh(1,1) =  7.d0 / 20.d0
    table_lh(1,2) =  1.d0 / 20.d0
    table_lh(1,3) =  3.d0 / 20.d0
    table_lh(1,4) =- 1.d0 / 30.d0
    table_lh(1,5) =- 1.d0 / 2.d0
    table_lh(1,6) =  1.d0 / 12.d0
    table_lh(1,7) =  1.d0 / 2.d0
    table_lh(1,8) =- 1.d0 / 12.d0

    table_lh(2,1) =  3.d0 / 20.d0
    table_lh(2,2) =  1.d0 / 30.d0
    table_lh(2,3) =  7.d0 / 20.d0
    table_lh(2,4) =- 1.d0 / 20.d0
    table_lh(2,5) =- 1.d0 / 2.d0
    table_lh(2,6) =- 1.d0 / 12.d0
    table_lh(2,7) =  1.d0 / 2.d0
    table_lh(2,8) =  1.d0 / 12.d0

    table_lh(3,1) =- 1.d0 / 2.d0
    table_lh(3,2) =- 1.d0 / 12.d0
    table_lh(3,3) =- 1.d0 / 2.d0
    table_lh(3,4) =  1.d0 / 12.d0
    table_lh(3,5) =  1.d0
    table_lh(3,6) =  0.d0
    table_lh(3,7) =  1.d0
    table_lh(3,8) =  0.d0

    table_lh(4,1) =  1.d0 / 2.d0
    table_lh(4,2) =  1.d0 / 12.d0
    table_lh(4,3) =  1.d0 / 2.d0
    table_lh(4,4) =- 1.d0 / 12.d0
    table_lh(4,5) =- 1.d0
    table_lh(4,6) =  0.d0
    table_lh(4,7) =  1.d0
    table_lh(4,8) =  0.d0

    !  *** triple product hhh ***

    table_hhh(1,1,1) = 43.d0 / 140.d0
    table_hhh(1,1,2) = 97.d0 / 2520.d0
    table_hhh(1,1,3) =  9.d0 / 140.d0
    table_hhh(1,1,4) =-43.d0 / 2520.d0
    table_hhh(1,1,5) =- 1.d0 / 3.d0
    table_hhh(1,1,6) =  5.d0 / 42.d0
    table_hhh(1,1,7) =  1.d0 / 3.d0
    table_hhh(1,1,8) =-17.d0 / 210.d0

    table_hhh(1,2,1) = 97.d0 / 2520.d0
    table_hhh(1,2,2) =  2.d0 / 315.d0
    table_hhh(1,2,3) =  1.d0 / 72.d0
    table_hhh(1,2,4) =- 1.d0 / 280.d0
    table_hhh(1,2,5) =- 5.d0 / 84.d0
    table_hhh(1,2,6) =  1.d0 / 168.d0
    table_hhh(1,2,7) =  5.d0 / 84.d0
    table_hhh(1,2,8) =-11.d0 / 840.d0

    table_hhh(1,3,1) =  9.d0 / 140.d0
    table_hhh(1,3,2) =  1.d0 / 72.d0
    table_hhh(1,3,3) =  9.d0 / 140.d0
    table_hhh(1,3,4) =- 1.d0 / 72.d0
    table_hhh(1,3,5) =- 1.d0 / 6.d0
    table_hhh(1,3,6) =- 2.d0 / 105.d0
    table_hhh(1,3,7) =  1.d0 / 6.d0
    table_hhh(1,3,8) =- 2.d0 / 105.d0

    table_hhh(1,4,1) =-43.d0 / 2520.d0
    table_hhh(1,4,2) =- 1.d0 / 280.d0
    table_hhh(1,4,3) =- 1.d0 / 72.d0
    table_hhh(1,4,4) =  1.d0 / 315.d0
    table_hhh(1,4,5) = 17.d0 / 420.d0
    table_hhh(1,4,6) =  1.d0 / 280.d0
    table_hhh(1,4,7) =-17.d0 / 420.d0
    table_hhh(1,4,8) =  1.d0 / 168.d0

    table_hhh(1,5,1) =- 1.d0 / 3.d0
    table_hhh(1,5,2) =- 5.d0 / 84.d0
    table_hhh(1,5,3) =- 1.d0 / 6.d0
    table_hhh(1,5,4) = 17.d0 / 420.d0
    table_hhh(1,5,5) =  3.d0 / 5.d0
    table_hhh(1,5,6) =- 1.d0 / 70.d0
    table_hhh(1,5,7) =- 3.d0 / 5.d0
    table_hhh(1,5,8) =  4.d0 / 35.d0

    table_hhh(1,6,1) =  5.d0 / 42.d0
    table_hhh(1,6,2) =  1.d0 / 168.d0
    table_hhh(1,6,3) =- 2.d0 / 105.d0
    table_hhh(1,6,4) =  1.d0 / 280.d0
    table_hhh(1,6,5) =- 1.d0 / 70.d0
    table_hhh(1,6,6) = 43.d0 / 420.d0
    table_hhh(1,6,7) =  1.d0 / 70.d0
    table_hhh(1,6,8) =- 1.d0 / 60.d0

    table_hhh(1,7,1) =  1.d0 / 3.d0
    table_hhh(1,7,2) =  5.d0 / 84.d0
    table_hhh(1,7,3) =  1.d0 / 6.d0
    table_hhh(1,7,4) =-17.d0 / 420.d0
    table_hhh(1,7,5) =- 3.d0 / 5.d0
    table_hhh(1,7,6) =  1.d0 / 70.d0
    table_hhh(1,7,7) =  3.d0 / 5.d0
    table_hhh(1,7,8) =- 4.d0 / 35.d0

    table_hhh(1,8,1) =-17.d0 / 210.d0
    table_hhh(1,8,2) =-11.d0 / 840.d0
    table_hhh(1,8,3) =- 2.d0 / 105.d0
    table_hhh(1,8,4) =  1.d0 / 168.d0
    table_hhh(1,8,5) =  4.d0 / 35.d0
    table_hhh(1,8,6) =- 1.d0 / 60.d0
    table_hhh(1,8,7) =- 4.d0 / 35.d0
    table_hhh(1,8,8) = 13.d0 / 420.d0

!   -----

    table_hhh(2,1,1) = 97.d0 / 2520.d0
    table_hhh(2,1,2) =  2.d0 / 315.d0
    table_hhh(2,1,3) =  1.d0 / 72.d0
    table_hhh(2,1,4) =- 1.d0 / 280.d0
    table_hhh(2,1,5) =- 5.d0 / 84.d0
    table_hhh(2,1,6) =  1.d0 / 168.d0
    table_hhh(2,1,7) =  5.d0 / 84.d0
    table_hhh(2,1,8) =-11.d0 / 840.d0

    table_hhh(2,2,1) =  2.d0 / 315.d0
    table_hhh(2,2,2) =  1.d0 / 840.d0
    table_hhh(2,2,3) =  1.d0 / 315.d0
    table_hhh(2,2,4) =- 1.d0 / 1260.d0
    table_hhh(2,2,5) =- 1.d0 / 84.d0
    table_hhh(2,2,6) =  0.d0 
    table_hhh(2,2,7) =  1.d0 / 84.d0
    table_hhh(2,2,8) =- 1.d0 / 420.d0

    table_hhh(2,3,1) =  1.d0 / 72.d0
    table_hhh(2,3,2) =  1.d0 / 315.d0
    table_hhh(2,3,3) = 43.d0 / 2520.d0
    table_hhh(2,3,4) =- 1.d0 / 280.d0
    table_hhh(2,3,5) =-17.d0 / 420.d0
    table_hhh(2,3,6) =- 1.d0 / 168.d0
    table_hhh(2,3,7) = 17.d0 / 420.d0
    table_hhh(2,3,8) =- 1.d0 / 280.d0

    table_hhh(2,4,1) = -1.d0 / 280.d0
    table_hhh(2,4,2) =- 1.d0 / 1260.d0
    table_hhh(2,4,3) = -1.d0 / 280.d0
    table_hhh(2,4,4) =  1.d0 / 1260.d0
    table_hhh(2,4,5) =  1.d0 / 105.d0
    table_hhh(2,4,6) =  1.d0 / 840.d0
    table_hhh(2,4,7) =- 1.d0 / 105.d0
    table_hhh(2,4,8) =  1.d0 / 840.d0

    table_hhh(2,5,1) =- 5.d0 / 84.d0
    table_hhh(2,5,2) =- 1.d0 / 84.d0
    table_hhh(2,5,3) =-17.d0 / 420.d0
    table_hhh(2,5,4) =  1.d0 / 105.d0
    table_hhh(2,5,5) =  9.d0 / 70.d0
    table_hhh(2,5,6) =  1.d0 / 140.d0
    table_hhh(2,5,7) =- 9.d0 / 70.d0
    table_hhh(2,5,8) =  3.d0 / 140.d0

    table_hhh(2,6,1) =  1.d0 / 168.d0
    table_hhh(2,6,2) =  0.d0
    table_hhh(2,6,3) =- 1.d0 / 168.d0
    table_hhh(2,6,4) =  1.d0 / 840.d0
    table_hhh(2,6,5) =  1.d0 / 140.d0 
    table_hhh(2,6,6) =  1.d0 / 120.d0
    table_hhh(2,6,7) =- 1.d0 / 140.d0
    table_hhh(2,6,8) =- 1.d0 / 840.d0

    table_hhh(2,7,1) =  5.d0 / 84.d0
    table_hhh(2,7,2) =  1.d0 / 84.d0
    table_hhh(2,7,3) = 17.d0 / 420.d0
    table_hhh(2,7,4) =- 1.d0 / 105.d0
    table_hhh(2,7,5) =- 9.d0 / 70.d0
    table_hhh(2,7,6) =- 1.d0 / 140.d0
    table_hhh(2,7,7) =  9.d0 / 70.d0
    table_hhh(2,7,8) =- 3.d0 / 140.d0

    table_hhh(2,8,1) =-11.d0 / 840.d0
    table_hhh(2,8,2) =- 1.d0 / 420.d0
    table_hhh(2,8,3) =- 1.d0 / 280.d0
    table_hhh(2,8,4) =  1.d0 / 840.d0
    table_hhh(2,8,5) =  3.d0 / 140.d0
    table_hhh(2,8,6) =- 1.d0 / 840.d0
    table_hhh(2,8,7) =- 3.d0 / 140.d0
    table_hhh(2,8,8) =  1.d0 / 168.d0

!   -----

    table_hhh(3,1,1) =  9.d0 / 140.d0
    table_hhh(3,1,2) =  1.d0 / 72.d0
    table_hhh(3,1,3) =  9.d0 / 140.d0
    table_hhh(3,1,4) =- 1.d0 / 72.d0
    table_hhh(3,1,5) =- 1.d0 / 6.d0
    table_hhh(3,1,6) =- 2.d0 / 105.d0
    table_hhh(3,1,7) =  1.d0 / 6.d0
    table_hhh(3,1,8) =- 2.d0 / 105.d0

    table_hhh(3,2,1) =  1.d0 / 72.d0
    table_hhh(3,2,2) =  1.d0 / 315.d0
    table_hhh(3,2,3) = 43.d0 / 2520.d0
    table_hhh(3,2,4) =- 1.d0 / 280.d0
    table_hhh(3,2,5) =-17.d0 / 420.d0
    table_hhh(3,2,6) =- 1.d0 / 168.d0
    table_hhh(3,2,7) = 17.d0 / 420.d0
    table_hhh(3,2,8) =- 1.d0 / 280.d0

    table_hhh(3,3,1) =  9.d0 / 140.d0
    table_hhh(3,3,2) = 43.d0 / 2520.d0
    table_hhh(3,3,3) = 43.d0 / 140.d0
    table_hhh(3,3,4) =-97.d0 / 2520.d0
    table_hhh(3,3,5) =- 1.d0 / 3.d0
    table_hhh(3,3,6) =-17.d0 / 210.d0
    table_hhh(3,3,7) =  1.d0 / 3.d0
    table_hhh(3,3,8) =  5.d0 / 42.d0

    table_hhh(3,4,1) =- 1.d0 / 72.d0
    table_hhh(3,4,2) =- 1.d0 / 280.d0
    table_hhh(3,4,3) =-97.d0 / 2520.d0
    table_hhh(3,4,4) =  2.d0 / 315.d0
    table_hhh(3,4,5) =  5.d0 / 84.d0
    table_hhh(3,4,6) = 11.d0 / 840.d0
    table_hhh(3,4,7) =- 5.d0 / 84.d0
    table_hhh(3,4,8) =- 1.d0 / 168.d0

    table_hhh(3,5,1) =- 1.d0 / 6.d0
    table_hhh(3,5,2) =-17.d0 / 420.d0
    table_hhh(3,5,3) =- 1.d0 / 3.d0
    table_hhh(3,5,4) =  5.d0 / 84.d0
    table_hhh(3,5,5) =  3.d0 / 5.d0
    table_hhh(3,5,6) =  4.d0 / 35.d0
    table_hhh(3,5,7) =- 3.d0 / 5.d0
    table_hhh(3,5,8) =- 1.d0 / 70.d0

    table_hhh(3,6,1) =- 2.d0 / 105.d0
    table_hhh(3,6,2) =- 1.d0 / 168.d0
    table_hhh(3,6,3) =-17.d0 / 210.d0
    table_hhh(3,6,4) = 11.d0 / 840.d0
    table_hhh(3,6,5) =  4.d0 / 35.d0
    table_hhh(3,6,6) = 13.d0 / 420.d0
    table_hhh(3,6,7) =- 4.d0 / 35.d0
    table_hhh(3,6,8) =- 1.d0 / 60.d0

    table_hhh(3,7,1) =  1.d0 / 6.d0
    table_hhh(3,7,2) = 17.d0 / 420.d0
    table_hhh(3,7,3) =  1.d0 / 3.d0
    table_hhh(3,7,4) =- 5.d0 / 84.d0
    table_hhh(3,7,5) =- 3.d0 / 5.d0
    table_hhh(3,7,6) =- 4.d0 / 35.d0
    table_hhh(3,7,7) =  3.d0 / 5.d0
    table_hhh(3,7,8) =  1.d0 / 70.d0

    table_hhh(3,8,1) =- 2.d0 / 105.d0
    table_hhh(3,8,2) =- 1.d0 / 280.d0
    table_hhh(3,8,3) =  5.d0 / 42.d0
    table_hhh(3,8,4) =- 1.d0 / 168.d0
    table_hhh(3,8,5) =- 1.d0 / 70.d0
    table_hhh(3,8,6) =- 1.d0 / 60.d0
    table_hhh(3,8,7) =  1.d0 / 70.d0
    table_hhh(3,8,8) = 43.d0 / 420.d0

!   -----

    table_hhh(4,1,1) =-43.d0 / 2520.d0
    table_hhh(4,1,2) =- 1.d0 / 280.d0
    table_hhh(4,1,3) =- 1.d0 / 72.d0
    table_hhh(4,1,4) =  1.d0 / 315.d0
    table_hhh(4,1,5) = 17.d0 / 420.d0
    table_hhh(4,1,6) =  1.d0 / 280.d0
    table_hhh(4,1,7) =-17.d0 / 420.d0
    table_hhh(4,1,8) =  1.d0 / 168.d0

    table_hhh(4,2,1) = -1.d0 / 280.d0
    table_hhh(4,2,2) =- 1.d0 / 1260.d0
    table_hhh(4,2,3) = -1.d0 / 280.d0
    table_hhh(4,2,4) =  1.d0 / 1260.d0
    table_hhh(4,2,5) =  1.d0 / 105.d0
    table_hhh(4,2,6) =  1.d0 / 840.d0
    table_hhh(4,2,7) =- 1.d0 / 105.d0
    table_hhh(4,2,8) =  1.d0 / 840.d0

    table_hhh(4,3,1) =- 1.d0 / 72.d0
    table_hhh(4,3,2) =- 1.d0 / 280.d0
    table_hhh(4,3,3) =-97.d0 / 2520.d0
    table_hhh(4,3,4) =  2.d0 / 315.d0
    table_hhh(4,3,5) =  5.d0 / 84.d0
    table_hhh(4,3,6) = 11.d0 / 840.d0
    table_hhh(4,3,7) =- 5.d0 / 84.d0
    table_hhh(4,3,8) =- 1.d0 / 168.d0

    table_hhh(4,4,1) =  1.d0 / 315.d0
    table_hhh(4,4,2) =  1.d0 / 1260.d0
    table_hhh(4,4,3) =  2.d0 / 315.d0
    table_hhh(4,4,4) =- 1.d0 / 840.d0
    table_hhh(4,4,5) =- 1.d0 / 84.d0
    table_hhh(4,4,6) =- 1.d0 / 420.d0
    table_hhh(4,4,7) =  1.d0 / 84.d0
    table_hhh(4,4,8) =  0.d0

    table_hhh(4,5,1) = 17.d0 / 420.d0
    table_hhh(4,5,2) =  1.d0 / 105.d0
    table_hhh(4,5,3) =  5.d0 / 84.d0
    table_hhh(4,5,4) =- 1.d0 / 84.d0
    table_hhh(4,5,5) =- 9.d0 / 70.d0
    table_hhh(4,5,6) =- 3.d0 / 140.d0
    table_hhh(4,5,7) =  9.d0 / 70.d0
    table_hhh(4,5,8) =- 1.d0 / 140.d0

    table_hhh(4,6,1) =  1.d0 / 280.d0
    table_hhh(4,6,2) =  1.d0 / 840.d0
    table_hhh(4,6,3) = 11.d0 / 840.d0
    table_hhh(4,6,4) =- 1.d0 / 420.d0
    table_hhh(4,6,5) =- 3.d0 / 140.d0
    table_hhh(4,6,6) =- 1.d0 / 168.d0
    table_hhh(4,6,7) =  3.d0 / 140.d0
    table_hhh(4,6,8) =  1.d0 / 840.d0

    table_hhh(4,7,1) =-17.d0 / 420.d0
    table_hhh(4,7,2) =- 1.d0 / 105.d0
    table_hhh(4,7,3) =- 5.d0 / 84.d0
    table_hhh(4,7,4) =  1.d0 / 84.d0
    table_hhh(4,7,5) =  9.d0 / 70.d0
    table_hhh(4,7,6) =  3.d0 / 140.d0
    table_hhh(4,7,7) =- 9.d0 / 70.d0
    table_hhh(4,7,8) =  1.d0 / 140.d0

    table_hhh(4,8,1) =  1.d0 / 168.d0
    table_hhh(4,8,2) =  1.d0 / 840.d0
    table_hhh(4,8,3) =- 1.d0 / 168.d0
    table_hhh(4,8,4) =  0.d0
    table_hhh(4,8,5) =- 1.d0 / 140.d0
    table_hhh(4,8,6) =  1.d0 / 840.d0
    table_hhh(4,8,7) =  1.d0 / 140.d0
    table_hhh(4,8,8) =- 1.d0 / 120.d0

!   -----

    table_hhh(5,1,1) =- 1.d0 / 3.d0
    table_hhh(5,1,2) =- 5.d0 / 84.d0
    table_hhh(5,1,3) =- 1.d0 / 6.d0
    table_hhh(5,1,4) = 17.d0 / 420.d0
    table_hhh(5,1,5) =  3.d0 / 5.d0
    table_hhh(5,1,6) =- 1.d0 / 70.d0
    table_hhh(5,1,7) =- 3.d0 / 5.d0
    table_hhh(5,1,8) =  4.d0 / 35.d0

    table_hhh(5,2,1) =- 5.d0 / 84.d0
    table_hhh(5,2,2) =- 1.d0 / 84.d0 
    table_hhh(5,2,3) =-17.d0 / 420.d0
    table_hhh(5,2,4) =  1.d0 / 105.d0
    table_hhh(5,2,5) =  9.d0 / 70.d0
    table_hhh(5,2,6) =  1.d0 / 140.d0
    table_hhh(5,2,7) =- 9.d0 / 70.d0
    table_hhh(5,2,8) =  3.d0 / 140.d0

    table_hhh(5,3,1) =- 1.d0 / 6.d0
    table_hhh(5,3,2) =-17.d0 / 420.d0
    table_hhh(5,3,3) =- 1.d0 / 3.d0
    table_hhh(5,3,4) =  5.d0 / 84.d0
    table_hhh(5,3,5) =  3.d0 / 5.d0
    table_hhh(5,3,6) =  4.d0 / 35.d0
    table_hhh(5,3,7) =- 3.d0 / 5.d0
    table_hhh(5,3,8) =- 1.d0 / 70.d0

    table_hhh(5,4,1) = 17.d0 / 420.d0
    table_hhh(5,4,2) =  1.d0 / 105.d0
    table_hhh(5,4,3) =  5.d0 / 84.d0
    table_hhh(5,4,4) =- 1.d0 / 84.d0 
    table_hhh(5,4,5) =- 9.d0 / 70.d0
    table_hhh(5,4,6) =- 3.d0 / 140.d0
    table_hhh(5,4,7) =  9.d0 / 70.d0
    table_hhh(5,4,8) =- 1.d0 / 140.d0

    table_hhh(5,5,1) =  3.d0 / 5.d0
    table_hhh(5,5,2) =  9.d0 / 70.d0
    table_hhh(5,5,3) =  3.d0 / 5.d0
    table_hhh(5,5,4) =- 9.d0 / 70.d0
    table_hhh(5,5,5) =-54.d0 / 35.d0
    table_hhh(5,5,6) =- 6.d0 / 35.d0
    table_hhh(5,5,7) = 54.d0 / 35.d0
    table_hhh(5,5,8) =- 6.d0 / 35.d0

    table_hhh(5,6,1) =- 1.d0 / 70.d0
    table_hhh(5,6,2) =  1.d0 / 140.d0
    table_hhh(5,6,3) =  4.d0 / 35.d0
    table_hhh(5,6,4) =- 3.d0 / 140.d0
    table_hhh(5,6,5) =- 6.d0 / 35.d0
    table_hhh(5,6,6) =- 3.d0 / 35.d0
    table_hhh(5,6,7) =  6.d0 / 35.d0
    table_hhh(5,6,8) =  1.d0 / 70.d0

    table_hhh(5,7,1) =- 3.d0 / 5.d0
    table_hhh(5,7,2) =- 9.d0 / 70.d0
    table_hhh(5,7,3) =- 3.d0 / 5.d0
    table_hhh(5,7,4) =  9.d0 / 70.d0
    table_hhh(5,7,5) = 54.d0 / 35.d0
    table_hhh(5,7,6) =  6.d0 / 35.d0
    table_hhh(5,7,7) =-54.d0 / 35.d0
    table_hhh(5,7,8) =  6.d0 / 35.d0

    table_hhh(5,8,1) =  4.d0 / 35.d0
    table_hhh(5,8,2) =  3.d0 / 140.d0
    table_hhh(5,8,3) =- 1.d0 / 70.d0
    table_hhh(5,8,4) =- 1.d0 / 140.d0
    table_hhh(5,8,5) =- 6.d0 / 35.d0
    table_hhh(5,8,6) =  1.d0 / 70.d0
    table_hhh(5,8,7) =  6.d0 / 35.d0
    table_hhh(5,8,8) =- 3.d0 / 35.d0

!   -----

    table_hhh(6,1,1) =  5.d0 / 42.d0
    table_hhh(6,1,2) =  1.d0 / 168.d0
    table_hhh(6,1,3) =- 2.d0 / 105.d0
    table_hhh(6,1,4) =  1.d0 / 280.d0
    table_hhh(6,1,5) =- 1.d0 / 70.d0
    table_hhh(6,1,6) = 43.d0 / 420.d0
    table_hhh(6,1,7) =  1.d0 / 70.d0
    table_hhh(6,1,8) =- 1.d0 / 60.d0

    table_hhh(6,2,1) =  1.d0 / 168.d0
    table_hhh(6,2,2) =  0.d0         
    table_hhh(6,2,3) =- 1.d0 / 168.d0
    table_hhh(6,2,4) =  1.d0 / 840.d0
    table_hhh(6,2,5) =  1.d0 / 140.d0
    table_hhh(6,2,6) =  1.d0 / 120.d0
    table_hhh(6,2,7) =- 1.d0 / 140.d0
    table_hhh(6,2,8) =- 1.d0 / 840.d0

    table_hhh(6,3,1) =- 2.d0 / 105.d0
    table_hhh(6,3,2) =- 1.d0 / 168.d0
    table_hhh(6,3,3) =-17.d0 / 210.d0
    table_hhh(6,3,4) = 11.d0 / 840.d0
    table_hhh(6,3,5) =  4.d0 / 35.d0
    table_hhh(6,3,6) = 13.d0 / 420.d0
    table_hhh(6,3,7) =- 4.d0 / 35.d0
    table_hhh(6,3,8) =- 1.d0 / 60.d0

    table_hhh(6,4,1) =  1.d0 / 280.d0
    table_hhh(6,4,2) =  1.d0 / 840.d0
    table_hhh(6,4,3) = 11.d0 / 840.d0
    table_hhh(6,4,4) =- 1.d0 / 420.d0
    table_hhh(6,4,5) =- 3.d0 / 140.d0
    table_hhh(6,4,6) =- 1.d0 / 168.d0
    table_hhh(6,4,7) =  3.d0 / 140.d0
    table_hhh(6,4,8) =  1.d0 / 840.d0

    table_hhh(6,5,1) =- 1.d0 / 70.d0
    table_hhh(6,5,2) =  1.d0 / 140.d0
    table_hhh(6,5,3) =  4.d0 / 35.d0
    table_hhh(6,5,4) =- 3.d0 / 140.d0
    table_hhh(6,5,5) =- 6.d0 / 35.d0
    table_hhh(6,5,6) =- 3.d0 / 35.d0
    table_hhh(6,5,7) =  6.d0 / 35.d0
    table_hhh(6,5,8) =  1.d0 / 70.d0

    table_hhh(6,6,1) = 43.d0 / 420.d0
    table_hhh(6,6,2) =  1.d0 / 120.d0
    table_hhh(6,6,3) = 13.d0 / 420.d0
    table_hhh(6,6,4) =- 1.d0 / 168.d0
    table_hhh(6,6,5) =- 3.d0 / 35.d0
    table_hhh(6,6,6) =  2.d0 / 35.d0
    table_hhh(6,6,7) =  3.d0 / 35.d0
    table_hhh(6,6,8) =- 1.d0 / 105.d0

    table_hhh(6,7,1) =  1.d0 / 70.d0
    table_hhh(6,7,2) =- 1.d0 / 140.d0
    table_hhh(6,7,3) =- 4.d0 / 35.d0
    table_hhh(6,7,4) =  3.d0 / 140.d0
    table_hhh(6,7,5) =  6.d0 / 35.d0
    table_hhh(6,7,6) =- 3.d0 / 35.d0
    table_hhh(6,7,7) =- 6.d0 / 35.d0
    table_hhh(6,7,8) =  1.d0 / 70.d0

    table_hhh(6,8,1) =- 1.d0 / 60.d0
    table_hhh(6,8,2) =- 1.d0 / 840.d0
    table_hhh(6,8,3) =- 1.d0 / 60.d0
    table_hhh(6,8,4) =  1.d0 / 840.d0
    table_hhh(6,8,5) =  1.d0 / 70.d0
    table_hhh(6,8,6) =- 1.d0 / 105.d0
    table_hhh(6,8,7) =- 1.d0 / 70.d0
    table_hhh(6,8,8) =- 1.d0 / 105.d0

!   -----

    table_hhh(7,1,1) =  1.d0 / 3.d0
    table_hhh(7,1,2) =  5.d0 / 84.d0
    table_hhh(7,1,3) =  1.d0 / 6.d0
    table_hhh(7,1,4) =-17.d0 / 420.d0
    table_hhh(7,1,5) =- 3.d0 / 5.d0
    table_hhh(7,1,6) =  1.d0 / 70.d0
    table_hhh(7,1,7) =  3.d0 / 5.d0
    table_hhh(7,1,8) =- 4.d0 / 35.d0

    table_hhh(7,2,1) =  5.d0 / 84.d0
    table_hhh(7,2,2) =  1.d0 / 84.d0 
    table_hhh(7,2,3) = 17.d0 / 420.d0
    table_hhh(7,2,4) =- 1.d0 / 105.d0
    table_hhh(7,2,5) =- 9.d0 / 70.d0
    table_hhh(7,2,6) =- 1.d0 / 140.d0
    table_hhh(7,2,7) =  9.d0 / 70.d0
    table_hhh(7,2,8) =- 3.d0 / 140.d0

    table_hhh(7,3,1) =  1.d0 / 6.d0
    table_hhh(7,3,2) = 17.d0 / 420.d0
    table_hhh(7,3,3) =  1.d0 / 3.d0
    table_hhh(7,3,4) =- 5.d0 / 84.d0
    table_hhh(7,3,5) =- 3.d0 / 5.d0
    table_hhh(7,3,6) =- 4.d0 / 35.d0
    table_hhh(7,3,7) =  3.d0 / 5.d0
    table_hhh(7,3,8) =  1.d0 / 70.d0

    table_hhh(7,4,1) =-17.d0 / 420.d0
    table_hhh(7,4,2) =- 1.d0 / 105.d0
    table_hhh(7,4,3) =- 5.d0 / 84.d0
    table_hhh(7,4,4) =  1.d0 / 84.d0 
    table_hhh(7,4,5) =  9.d0 / 70.d0
    table_hhh(7,4,6) =  3.d0 / 140.d0
    table_hhh(7,4,7) =- 9.d0 / 70.d0
    table_hhh(7,4,8) =  1.d0 / 140.d0

    table_hhh(7,5,1) =- 3.d0 / 5.d0
    table_hhh(7,5,2) =- 9.d0 / 70.d0
    table_hhh(7,5,3) =- 3.d0 / 5.d0
    table_hhh(7,5,4) =  9.d0 / 70.d0
    table_hhh(7,5,5) = 54.d0 / 35.d0
    table_hhh(7,5,6) =  6.d0 / 35.d0
    table_hhh(7,5,7) =-54.d0 / 35.d0
    table_hhh(7,5,8) =  6.d0 / 35.d0

    table_hhh(7,6,1) =  1.d0 / 70.d0
    table_hhh(7,6,2) =- 1.d0 / 140.d0
    table_hhh(7,6,3) =- 4.d0 / 35.d0
    table_hhh(7,6,4) =  3.d0 / 140.d0
    table_hhh(7,6,5) =  6.d0 / 35.d0
    table_hhh(7,6,6) =  3.d0 / 35.d0
    table_hhh(7,6,7) =- 6.d0 / 35.d0
    table_hhh(7,6,8) =- 1.d0 / 70.d0

    table_hhh(7,7,1) =  3.d0 / 5.d0
    table_hhh(7,7,2) =  9.d0 / 70.d0
    table_hhh(7,7,3) =  3.d0 / 5.d0
    table_hhh(7,7,4) =- 9.d0 / 70.d0
    table_hhh(7,7,5) =-54.d0 / 35.d0
    table_hhh(7,7,6) =- 6.d0 / 35.d0
    table_hhh(7,7,7) = 54.d0 / 35.d0
    table_hhh(7,7,8) =- 6.d0 / 35.d0

    table_hhh(7,8,1) =- 4.d0 / 35.d0
    table_hhh(7,8,2) =- 3.d0 / 140.d0
    table_hhh(7,8,3) =  1.d0 / 70.d0
    table_hhh(7,8,4) =  1.d0 / 140.d0
    table_hhh(7,8,5) =  6.d0 / 35.d0
    table_hhh(7,8,6) =- 1.d0 / 70.d0
    table_hhh(7,8,7) =- 6.d0 / 35.d0
    table_hhh(7,8,8) =  3.d0 / 35.d0

!   -----

    table_hhh(8,1,1) =-17.d0 / 210.d0
    table_hhh(8,1,2) =-11.d0 / 840.d0
    table_hhh(8,1,3) =- 2.d0 / 105.d0
    table_hhh(8,1,4) =  1.d0 / 168.d0
    table_hhh(8,1,5) =  4.d0 / 35.d0
    table_hhh(8,1,6) =- 1.d0 / 60.d0
    table_hhh(8,1,7) =- 4.d0 / 35.d0
    table_hhh(8,1,8) = 13.d0 / 420.d0

    table_hhh(8,2,1) =-11.d0 / 840.d0
    table_hhh(8,2,2) =- 1.d0 / 420.d0
    table_hhh(8,2,3) =- 1.d0 / 280.d0
    table_hhh(8,2,4) =  1.d0 / 840.d0
    table_hhh(8,2,5) =  3.d0 / 140.d0
    table_hhh(8,2,6) =- 1.d0 / 840.d0
    table_hhh(8,2,7) =- 3.d0 / 140.d0
    table_hhh(8,2,8) =  1.d0 / 168.d0

    table_hhh(8,3,1) =- 2.d0 / 105.d0
    table_hhh(8,3,2) =- 1.d0 / 280.d0
    table_hhh(8,3,3) =  5.d0 / 42.d0
    table_hhh(8,3,4) =- 1.d0 / 168.d0
    table_hhh(8,3,5) =- 1.d0 / 70.d0
    table_hhh(8,3,6) =- 1.d0 / 60.d0
    table_hhh(8,3,7) =  1.d0 / 70.d0
    table_hhh(8,3,8) = 43.d0 / 420.d0

    table_hhh(8,4,1) =  1.d0 / 168.d0
    table_hhh(8,4,2) =  1.d0 / 840.d0
    table_hhh(8,4,3) =- 1.d0 / 168.d0
    table_hhh(8,4,4) =  0.d0         
    table_hhh(8,4,5) =- 1.d0 / 140.d0
    table_hhh(8,4,6) =  1.d0 / 840.d0
    table_hhh(8,4,7) =  1.d0 / 140.d0
    table_hhh(8,4,8) =- 1.d0 / 120.d0

    table_hhh(8,5,1) =  4.d0 / 35.d0
    table_hhh(8,5,2) =  3.d0 / 140.d0
    table_hhh(8,5,3) =- 1.d0 / 70.d0
    table_hhh(8,5,4) =- 1.d0 / 140.d0
    table_hhh(8,5,5) =- 6.d0 / 35.d0
    table_hhh(8,5,6) =  1.d0 / 70.d0
    table_hhh(8,5,7) =  6.d0 / 35.d0
    table_hhh(8,5,8) =- 3.d0 / 35.d0

    table_hhh(8,6,1) =- 1.d0 / 60.d0
    table_hhh(8,6,2) =- 1.d0 / 840.d0
    table_hhh(8,6,3) =- 1.d0 / 60.d0
    table_hhh(8,6,4) =  1.d0 / 840.d0
    table_hhh(8,6,5) =  1.d0 / 70.d0
    table_hhh(8,6,6) =- 1.d0 / 105.d0
    table_hhh(8,6,7) =- 1.d0 / 70.d0
    table_hhh(8,6,8) =- 1.d0 / 105.d0

    table_hhh(8,7,1) =- 4.d0 / 35.d0
    table_hhh(8,7,2) =- 3.d0 / 140.d0
    table_hhh(8,7,3) =  1.d0 / 70.d0
    table_hhh(8,7,4) =  1.d0 / 140.d0
    table_hhh(8,7,5) =  6.d0 / 35.d0
    table_hhh(8,7,6) =- 3.d0 / 35.d0
    table_hhh(8,7,7) =- 6.d0 / 35.d0
    table_hhh(8,7,8) =  1.d0 / 70.d0

    table_hhh(8,8,1) = 13.d0 / 420.d0
    table_hhh(8,8,2) =  1.d0 / 168.d0
    table_hhh(8,8,3) = 43.d0 / 420.d0
    table_hhh(8,8,4) =- 1.d0 / 120.d0
    table_hhh(8,8,5) =- 3.d0 / 35.d0
    table_hhh(8,8,6) =- 1.d0 / 105.d0
    table_hhh(8,8,7) =  3.d0 / 35.d0
    table_hhh(8,8,8) =  2.d0 / 35.d0

    !  *** triple product lhh ***

    table_lhh(1,1,1) = 43.d0 / 140.d0
    table_lhh(1,1,2) = 97.d0 / 2520.d0
    table_lhh(1,1,3) =  9.d0 / 140.d0
    table_lhh(1,1,4) =-43.d0 / 2520.d0
    table_lhh(1,1,5) =- 1.d0 / 3.d0
    table_lhh(1,1,6) =  5.d0 / 42.d0
    table_lhh(1,1,7) =  1.d0 / 3.d0
    table_lhh(1,1,8) =-17.d0 / 210.d0

    table_lhh(1,2,1) = 97.d0 / 2520.d0
    table_lhh(1,2,2) =  2.d0 / 315.d0
    table_lhh(1,2,3) =  1.d0 / 72.d0
    table_lhh(1,2,4) =- 1.d0 / 280.d0
    table_lhh(1,2,5) =- 5.d0 / 84.d0
    table_lhh(1,2,6) =  1.d0 / 168.d0
    table_lhh(1,2,7) =  5.d0 / 84.d0
    table_lhh(1,2,8) =-11.d0 / 840.d0

    table_lhh(1,3,1) =  9.d0 / 140.d0
    table_lhh(1,3,2) =  1.d0 / 72.d0
    table_lhh(1,3,3) =  9.d0 / 140.d0
    table_lhh(1,3,4) =- 1.d0 / 72.d0
    table_lhh(1,3,5) =- 1.d0 / 6.d0
    table_lhh(1,3,6) =- 2.d0 / 105.d0
    table_lhh(1,3,7) =  1.d0 / 6.d0
    table_lhh(1,3,8) =- 2.d0 / 105.d0

    table_lhh(1,4,1) =-43.d0 / 2520.d0
    table_lhh(1,4,2) =- 1.d0 / 280.d0
    table_lhh(1,4,3) =- 1.d0 / 72.d0
    table_lhh(1,4,4) =  1.d0 / 315.d0
    table_lhh(1,4,5) = 17.d0 / 420.d0
    table_lhh(1,4,6) =  1.d0 / 280.d0
    table_lhh(1,4,7) =-17.d0 / 420.d0
    table_lhh(1,4,8) =  1.d0 / 168.d0

    table_lhh(1,5,1) =- 1.d0 / 3.d0
    table_lhh(1,5,2) =- 5.d0 / 84.d0
    table_lhh(1,5,3) =- 1.d0 / 6.d0
    table_lhh(1,5,4) = 17.d0 / 420.d0
    table_lhh(1,5,5) =  3.d0 / 5.d0
    table_lhh(1,5,6) =- 1.d0 / 70.d0
    table_lhh(1,5,7) =- 3.d0 / 5.d0
    table_lhh(1,5,8) =  4.d0 / 35.d0

    table_lhh(1,6,1) =  5.d0 / 42.d0
    table_lhh(1,6,2) =  1.d0 / 168.d0
    table_lhh(1,6,3) =- 2.d0 / 105.d0
    table_lhh(1,6,4) =  1.d0 / 280.d0
    table_lhh(1,6,5) =- 1.d0 / 70.d0
    table_lhh(1,6,6) = 43.d0 / 420.d0
    table_lhh(1,6,7) =  1.d0 / 70.d0
    table_lhh(1,6,8) =- 1.d0 / 60.d0

    table_lhh(1,7,1) =  1.d0 / 3.d0
    table_lhh(1,7,2) =  5.d0 / 84.d0
    table_lhh(1,7,3) =  1.d0 / 6.d0
    table_lhh(1,7,4) =-17.d0 / 420.d0
    table_lhh(1,7,5) =- 3.d0 / 5.d0
    table_lhh(1,7,6) =  1.d0 / 70.d0
    table_lhh(1,7,7) =  3.d0 / 5.d0
    table_lhh(1,7,8) =- 4.d0 / 35.d0

    table_lhh(1,8,1) =-17.d0 / 210.d0
    table_lhh(1,8,2) =-11.d0 / 840.d0
    table_lhh(1,8,3) =- 2.d0 / 105.d0
    table_lhh(1,8,4) =  1.d0 / 168.d0
    table_lhh(1,8,5) =  4.d0 / 35.d0
    table_lhh(1,8,6) =- 1.d0 / 60.d0
    table_lhh(1,8,7) =- 4.d0 / 35.d0
    table_lhh(1,8,8) = 13.d0 / 420.d0

!   -----

    table_lhh(2,1,1) = 97.d0 / 2520.d0
    table_lhh(2,1,2) =  2.d0 / 315.d0
    table_lhh(2,1,3) =  1.d0 / 72.d0
    table_lhh(2,1,4) =- 1.d0 / 280.d0
    table_lhh(2,1,5) =- 5.d0 / 84.d0
    table_lhh(2,1,6) =  1.d0 / 168.d0
    table_lhh(2,1,7) =  5.d0 / 84.d0
    table_lhh(2,1,8) =-11.d0 / 840.d0

    table_lhh(2,2,1) =  2.d0 / 315.d0
    table_lhh(2,2,2) =  1.d0 / 840.d0
    table_lhh(2,2,3) =  1.d0 / 315.d0
    table_lhh(2,2,4) =- 1.d0 / 1260.d0
    table_lhh(2,2,5) =- 1.d0 / 84.d0
    table_lhh(2,2,6) =  0.d0 
    table_lhh(2,2,7) =  1.d0 / 84.d0
    table_lhh(2,2,8) =- 1.d0 / 420.d0

    table_lhh(2,3,1) =  1.d0 / 72.d0
    table_lhh(2,3,2) =  1.d0 / 315.d0
    table_lhh(2,3,3) = 43.d0 / 2520.d0
    table_lhh(2,3,4) =- 1.d0 / 280.d0
    table_lhh(2,3,5) =-17.d0 / 420.d0
    table_lhh(2,3,6) =- 1.d0 / 168.d0
    table_lhh(2,3,7) = 17.d0 / 420.d0
    table_lhh(2,3,8) =- 1.d0 / 280.d0

    table_lhh(2,4,1) = -1.d0 / 280.d0
    table_lhh(2,4,2) =- 1.d0 / 1260.d0
    table_lhh(2,4,3) = -1.d0 / 280.d0
    table_lhh(2,4,4) =  1.d0 / 1260.d0
    table_lhh(2,4,5) =  1.d0 / 105.d0
    table_lhh(2,4,6) =  1.d0 / 840.d0
    table_lhh(2,4,7) =- 1.d0 / 105.d0
    table_lhh(2,4,8) =  1.d0 / 840.d0

    table_lhh(2,5,1) =- 5.d0 / 84.d0
    table_lhh(2,5,2) =- 1.d0 / 84.d0
    table_lhh(2,5,3) =-17.d0 / 420.d0
    table_lhh(2,5,4) =  1.d0 / 105.d0
    table_lhh(2,5,5) =  9.d0 / 70.d0
    table_lhh(2,5,6) =  1.d0 / 140.d0
    table_lhh(2,5,7) =- 9.d0 / 70.d0
    table_lhh(2,5,8) =  3.d0 / 140.d0

    table_lhh(2,6,1) =  1.d0 / 168.d0
    table_lhh(2,6,2) =  0.d0
    table_lhh(2,6,3) =- 1.d0 / 168.d0
    table_lhh(2,6,4) =  1.d0 / 840.d0
    table_lhh(2,6,5) =  1.d0 / 140.d0 
    table_lhh(2,6,6) =  1.d0 / 120.d0
    table_lhh(2,6,7) =- 1.d0 / 140.d0
    table_lhh(2,6,8) =- 1.d0 / 840.d0

    table_lhh(2,7,1) =  5.d0 / 84.d0
    table_lhh(2,7,2) =  1.d0 / 84.d0
    table_lhh(2,7,3) = 17.d0 / 420.d0
    table_lhh(2,7,4) =- 1.d0 / 105.d0
    table_lhh(2,7,5) =- 9.d0 / 70.d0
    table_lhh(2,7,6) =- 1.d0 / 140.d0
    table_lhh(2,7,7) =  9.d0 / 70.d0
    table_lhh(2,7,8) =- 3.d0 / 140.d0

    table_lhh(2,8,1) =-11.d0 / 840.d0
    table_lhh(2,8,2) =- 1.d0 / 420.d0
    table_lhh(2,8,3) =- 1.d0 / 280.d0
    table_lhh(2,8,4) =  1.d0 / 840.d0
    table_lhh(2,8,5) =  3.d0 / 140.d0
    table_lhh(2,8,6) =- 1.d0 / 840.d0
    table_lhh(2,8,7) =- 3.d0 / 140.d0
    table_lhh(2,8,8) =  1.d0 / 168.d0

!   -----

    table_lhh(3,1,1) =  9.d0 / 140.d0
    table_lhh(3,1,2) =  1.d0 / 72.d0
    table_lhh(3,1,3) =  9.d0 / 140.d0
    table_lhh(3,1,4) =- 1.d0 / 72.d0
    table_lhh(3,1,5) =- 1.d0 / 6.d0
    table_lhh(3,1,6) =- 2.d0 / 105.d0
    table_lhh(3,1,7) =  1.d0 / 6.d0
    table_lhh(3,1,8) =- 2.d0 / 105.d0

    table_lhh(3,2,1) =  1.d0 / 72.d0
    table_lhh(3,2,2) =  1.d0 / 315.d0
    table_lhh(3,2,3) = 43.d0 / 2520.d0
    table_lhh(3,2,4) =- 1.d0 / 280.d0
    table_lhh(3,2,5) =-17.d0 / 420.d0
    table_lhh(3,2,6) =- 1.d0 / 168.d0
    table_lhh(3,2,7) = 17.d0 / 420.d0
    table_lhh(3,2,8) =- 1.d0 / 280.d0

    table_lhh(3,3,1) =  9.d0 / 140.d0
    table_lhh(3,3,2) = 43.d0 / 2520.d0
    table_lhh(3,3,3) = 43.d0 / 140.d0
    table_lhh(3,3,4) =-97.d0 / 2520.d0
    table_lhh(3,3,5) =- 1.d0 / 3.d0
    table_lhh(3,3,6) =-17.d0 / 210.d0
    table_lhh(3,3,7) =  1.d0 / 3.d0
    table_lhh(3,3,8) =  5.d0 / 42.d0

    table_lhh(3,4,1) =- 1.d0 / 72.d0
    table_lhh(3,4,2) =- 1.d0 / 280.d0
    table_lhh(3,4,3) =-97.d0 / 2520.d0
    table_lhh(3,4,4) =  2.d0 / 315.d0
    table_lhh(3,4,5) =  5.d0 / 84.d0
    table_lhh(3,4,6) = 11.d0 / 840.d0
    table_lhh(3,4,7) =- 5.d0 / 84.d0
    table_lhh(3,4,8) =- 1.d0 / 168.d0

    table_lhh(3,5,1) =- 1.d0 / 6.d0
    table_lhh(3,5,2) =-17.d0 / 420.d0
    table_lhh(3,5,3) =- 1.d0 / 3.d0
    table_lhh(3,5,4) =  5.d0 / 84.d0
    table_lhh(3,5,5) =  3.d0 / 5.d0
    table_lhh(3,5,6) =  4.d0 / 35.d0
    table_lhh(3,5,7) =- 3.d0 / 5.d0
    table_lhh(3,5,8) =- 1.d0 / 70.d0

    table_lhh(3,6,1) =- 2.d0 / 105.d0
    table_lhh(3,6,2) =- 1.d0 / 168.d0
    table_lhh(3,6,3) =-17.d0 / 210.d0
    table_lhh(3,6,4) = 11.d0 / 840.d0
    table_lhh(3,6,5) =  4.d0 / 35.d0
    table_lhh(3,6,6) = 13.d0 / 420.d0
    table_lhh(3,6,7) =- 4.d0 / 35.d0
    table_lhh(3,6,8) =- 1.d0 / 60.d0

    table_lhh(3,7,1) =  1.d0 / 6.d0
    table_lhh(3,7,2) = 17.d0 / 420.d0
    table_lhh(3,7,3) =  1.d0 / 3.d0
    table_lhh(3,7,4) =- 5.d0 / 84.d0
    table_lhh(3,7,5) =- 3.d0 / 5.d0
    table_lhh(3,7,6) =- 4.d0 / 35.d0
    table_lhh(3,7,7) =  3.d0 / 5.d0
    table_lhh(3,7,8) =  1.d0 / 70.d0

    table_lhh(3,8,1) =- 2.d0 / 105.d0
    table_lhh(3,8,2) =- 1.d0 / 280.d0
    table_lhh(3,8,3) =  5.d0 / 42.d0
    table_lhh(3,8,4) =- 1.d0 / 168.d0
    table_lhh(3,8,5) =- 1.d0 / 70.d0
    table_lhh(3,8,6) =- 1.d0 / 60.d0
    table_lhh(3,8,7) =  1.d0 / 70.d0
    table_lhh(3,8,8) = 43.d0 / 420.d0

!   -----

    table_lhh(4,1,1) =-43.d0 / 2520.d0
    table_lhh(4,1,2) =- 1.d0 / 280.d0
    table_lhh(4,1,3) =- 1.d0 / 72.d0
    table_lhh(4,1,4) =  1.d0 / 315.d0
    table_lhh(4,1,5) = 17.d0 / 420.d0
    table_lhh(4,1,6) =  1.d0 / 280.d0
    table_lhh(4,1,7) =-17.d0 / 420.d0
    table_lhh(4,1,8) =  1.d0 / 168.d0

    table_lhh(4,2,1) = -1.d0 / 280.d0
    table_lhh(4,2,2) =- 1.d0 / 1260.d0
    table_lhh(4,2,3) = -1.d0 / 280.d0
    table_lhh(4,2,4) =  1.d0 / 1260.d0
    table_lhh(4,2,5) =  1.d0 / 105.d0
    table_lhh(4,2,6) =  1.d0 / 840.d0
    table_lhh(4,2,7) =- 1.d0 / 105.d0
    table_lhh(4,2,8) =  1.d0 / 840.d0

    table_lhh(4,3,1) =- 1.d0 / 72.d0
    table_lhh(4,3,2) =- 1.d0 / 280.d0
    table_lhh(4,3,3) =-97.d0 / 2520.d0
    table_lhh(4,3,4) =  2.d0 / 315.d0
    table_lhh(4,3,5) =  5.d0 / 84.d0
    table_lhh(4,3,6) = 11.d0 / 840.d0
    table_lhh(4,3,7) =- 5.d0 / 84.d0
    table_lhh(4,3,8) =- 1.d0 / 168.d0

    table_lhh(4,4,1) =  1.d0 / 315.d0
    table_lhh(4,4,2) =  1.d0 / 1260.d0
    table_lhh(4,4,3) =  2.d0 / 315.d0
    table_lhh(4,4,4) =- 1.d0 / 840.d0
    table_lhh(4,4,5) =- 1.d0 / 84.d0
    table_lhh(4,4,6) =- 1.d0 / 420.d0
    table_lhh(4,4,7) =  1.d0 / 84.d0
    table_lhh(4,4,8) =  0.d0

    table_lhh(4,5,1) = 17.d0 / 420.d0
    table_lhh(4,5,2) =  1.d0 / 105.d0
    table_lhh(4,5,3) =  5.d0 / 84.d0
    table_lhh(4,5,4) =- 1.d0 / 84.d0
    table_lhh(4,5,5) =- 9.d0 / 70.d0
    table_lhh(4,5,6) =- 3.d0 / 140.d0
    table_lhh(4,5,7) =  9.d0 / 70.d0
    table_lhh(4,5,8) =- 1.d0 / 140.d0

    table_lhh(4,6,1) =  1.d0 / 280.d0
    table_lhh(4,6,2) =  1.d0 / 840.d0
    table_lhh(4,6,3) = 11.d0 / 840.d0
    table_lhh(4,6,4) =- 1.d0 / 420.d0
    table_lhh(4,6,5) =- 3.d0 / 140.d0
    table_lhh(4,6,6) =- 1.d0 / 168.d0
    table_lhh(4,6,7) =  3.d0 / 140.d0
    table_lhh(4,6,8) =  1.d0 / 840.d0

    table_lhh(4,7,1) =-17.d0 / 420.d0
    table_lhh(4,7,2) =- 1.d0 / 105.d0
    table_lhh(4,7,3) =- 5.d0 / 84.d0
    table_lhh(4,7,4) =  1.d0 / 84.d0
    table_lhh(4,7,5) =  9.d0 / 70.d0
    table_lhh(4,7,6) =  3.d0 / 140.d0
    table_lhh(4,7,7) =- 9.d0 / 70.d0
    table_lhh(4,7,8) =  1.d0 / 140.d0

    table_lhh(4,8,1) =  1.d0 / 168.d0
    table_lhh(4,8,2) =  1.d0 / 840.d0
    table_lhh(4,8,3) =- 1.d0 / 168.d0
    table_lhh(4,8,4) =  0.d0
    table_lhh(4,8,5) =- 1.d0 / 140.d0
    table_lhh(4,8,6) =  1.d0 / 840.d0
    table_lhh(4,8,7) =  1.d0 / 140.d0
    table_lhh(4,8,8) =- 1.d0 / 120.d0

    do k=1,8
       do j=1,8
          do i=1,4
             table_hlh(j,i,k)=table_lhh(i,j,k)
          enddo
       enddo
    enddo

    !  *** triple product llh ***

    table_llh(1,1,1) = 43.d0 / 140.d0
    table_llh(1,1,2) = 97.d0 / 2520.d0
    table_llh(1,1,3) =  9.d0 / 140.d0
    table_llh(1,1,4) =-43.d0 / 2520.d0
    table_llh(1,1,5) =- 1.d0 / 3.d0
    table_llh(1,1,6) =  5.d0 / 42.d0
    table_llh(1,1,7) =  1.d0 / 3.d0
    table_llh(1,1,8) =-17.d0 / 210.d0

    table_llh(1,2,1) = 97.d0 / 2520.d0
    table_llh(1,2,2) =  2.d0 / 315.d0
    table_llh(1,2,3) =  1.d0 / 72.d0
    table_llh(1,2,4) =- 1.d0 / 280.d0
    table_llh(1,2,5) =- 5.d0 / 84.d0
    table_llh(1,2,6) =  1.d0 / 168.d0
    table_llh(1,2,7) =  5.d0 / 84.d0
    table_llh(1,2,8) =-11.d0 / 840.d0

    table_llh(1,3,1) =  9.d0 / 140.d0
    table_llh(1,3,2) =  1.d0 / 72.d0
    table_llh(1,3,3) =  9.d0 / 140.d0
    table_llh(1,3,4) =- 1.d0 / 72.d0
    table_llh(1,3,5) =- 1.d0 / 6.d0
    table_llh(1,3,6) =- 2.d0 / 105.d0
    table_llh(1,3,7) =  1.d0 / 6.d0
    table_llh(1,3,8) =- 2.d0 / 105.d0

    table_llh(1,4,1) =-43.d0 / 2520.d0
    table_llh(1,4,2) =- 1.d0 / 280.d0
    table_llh(1,4,3) =- 1.d0 / 72.d0
    table_llh(1,4,4) =  1.d0 / 315.d0
    table_llh(1,4,5) = 17.d0 / 420.d0
    table_llh(1,4,6) =  1.d0 / 280.d0
    table_llh(1,4,7) =-17.d0 / 420.d0
    table_llh(1,4,8) =  1.d0 / 168.d0

!   -----

    table_llh(2,1,1) = 97.d0 / 2520.d0
    table_llh(2,1,2) =  2.d0 / 315.d0
    table_llh(2,1,3) =  1.d0 / 72.d0
    table_llh(2,1,4) =- 1.d0 / 280.d0
    table_llh(2,1,5) =- 5.d0 / 84.d0
    table_llh(2,1,6) =  1.d0 / 168.d0
    table_llh(2,1,7) =  5.d0 / 84.d0
    table_llh(2,1,8) =-11.d0 / 840.d0

    table_llh(2,2,1) =  2.d0 / 315.d0
    table_llh(2,2,2) =  1.d0 / 840.d0
    table_llh(2,2,3) =  1.d0 / 315.d0
    table_llh(2,2,4) =- 1.d0 / 1260.d0
    table_llh(2,2,5) =- 1.d0 / 84.d0
    table_llh(2,2,6) =  0.d0 
    table_llh(2,2,7) =  1.d0 / 84.d0
    table_llh(2,2,8) =- 1.d0 / 420.d0

    table_llh(2,3,1) =  1.d0 / 72.d0
    table_llh(2,3,2) =  1.d0 / 315.d0
    table_llh(2,3,3) = 43.d0 / 2520.d0
    table_llh(2,3,4) =- 1.d0 / 280.d0
    table_llh(2,3,5) =-17.d0 / 420.d0
    table_llh(2,3,6) =- 1.d0 / 168.d0
    table_llh(2,3,7) = 17.d0 / 420.d0
    table_llh(2,3,8) =- 1.d0 / 280.d0

    table_llh(2,4,1) = -1.d0 / 280.d0
    table_llh(2,4,2) =- 1.d0 / 1260.d0
    table_llh(2,4,3) = -1.d0 / 280.d0
    table_llh(2,4,4) =  1.d0 / 1260.d0
    table_llh(2,4,5) =  1.d0 / 105.d0
    table_llh(2,4,6) =  1.d0 / 840.d0
    table_llh(2,4,7) =- 1.d0 / 105.d0
    table_llh(2,4,8) =  1.d0 / 840.d0

!   -----

    table_llh(3,1,1) =  9.d0 / 140.d0
    table_llh(3,1,2) =  1.d0 / 72.d0
    table_llh(3,1,3) =  9.d0 / 140.d0
    table_llh(3,1,4) =- 1.d0 / 72.d0
    table_llh(3,1,5) =- 1.d0 / 6.d0
    table_llh(3,1,6) =- 2.d0 / 105.d0
    table_llh(3,1,7) =  1.d0 / 6.d0
    table_llh(3,1,8) =- 2.d0 / 105.d0

    table_llh(3,2,1) =  1.d0 / 72.d0
    table_llh(3,2,2) =  1.d0 / 315.d0
    table_llh(3,2,3) = 43.d0 / 2520.d0
    table_llh(3,2,4) =- 1.d0 / 280.d0
    table_llh(3,2,5) =-17.d0 / 420.d0
    table_llh(3,2,6) =- 1.d0 / 168.d0
    table_llh(3,2,7) = 17.d0 / 420.d0
    table_llh(3,2,8) =- 1.d0 / 280.d0

    table_llh(3,3,1) =  9.d0 / 140.d0
    table_llh(3,3,2) = 43.d0 / 2520.d0
    table_llh(3,3,3) = 43.d0 / 140.d0
    table_llh(3,3,4) =-97.d0 / 2520.d0
    table_llh(3,3,5) =- 1.d0 / 3.d0
    table_llh(3,3,6) =-17.d0 / 210.d0
    table_llh(3,3,7) =  1.d0 / 3.d0
    table_llh(3,3,8) =  5.d0 / 42.d0

    table_llh(3,4,1) =- 1.d0 / 72.d0
    table_llh(3,4,2) =- 1.d0 / 280.d0
    table_llh(3,4,3) =-97.d0 / 2520.d0
    table_llh(3,4,4) =  2.d0 / 315.d0
    table_llh(3,4,5) =  5.d0 / 84.d0
    table_llh(3,4,6) = 11.d0 / 840.d0
    table_llh(3,4,7) =- 5.d0 / 84.d0
    table_llh(3,4,8) =- 1.d0 / 168.d0

!   -----

    table_llh(4,1,1) =-43.d0 / 2520.d0
    table_llh(4,1,2) =- 1.d0 / 280.d0
    table_llh(4,1,3) =- 1.d0 / 72.d0
    table_llh(4,1,4) =  1.d0 / 315.d0
    table_llh(4,1,5) = 17.d0 / 420.d0
    table_llh(4,1,6) =  1.d0 / 280.d0
    table_llh(4,1,7) =-17.d0 / 420.d0
    table_llh(4,1,8) =  1.d0 / 168.d0

    table_llh(4,2,1) = -1.d0 / 280.d0
    table_llh(4,2,2) =- 1.d0 / 1260.d0
    table_llh(4,2,3) = -1.d0 / 280.d0
    table_llh(4,2,4) =  1.d0 / 1260.d0
    table_llh(4,2,5) =  1.d0 / 105.d0
    table_llh(4,2,6) =  1.d0 / 840.d0
    table_llh(4,2,7) =- 1.d0 / 105.d0
    table_llh(4,2,8) =  1.d0 / 840.d0

    table_llh(4,3,1) =- 1.d0 / 72.d0
    table_llh(4,3,2) =- 1.d0 / 280.d0
    table_llh(4,3,3) =-97.d0 / 2520.d0
    table_llh(4,3,4) =  2.d0 / 315.d0
    table_llh(4,3,5) =  5.d0 / 84.d0
    table_llh(4,3,6) = 11.d0 / 840.d0
    table_llh(4,3,7) =- 5.d0 / 84.d0
    table_llh(4,3,8) =- 1.d0 / 168.d0

    table_llh(4,4,1) =  1.d0 / 315.d0
    table_llh(4,4,2) =  1.d0 / 1260.d0
    table_llh(4,4,3) =  2.d0 / 315.d0
    table_llh(4,4,4) =- 1.d0 / 840.d0
    table_llh(4,4,5) =- 1.d0 / 84.d0
    table_llh(4,4,6) =- 1.d0 / 420.d0
    table_llh(4,4,7) =  1.d0 / 84.d0
    table_llh(4,4,8) =  0.d0


 end subroutine table_initialize

 subroutine fem_integrate(id,x,a,da,b,db)

!-------------------------------------------------------
!
!   Calculate "\int_{psi_i}^{psi_{i+1}} function(psi) dpsi"
!      dpsi : mesh interval
!      a    : coefficient vector
!      b    : coefficient vector
!      u    : variable vector
!      w    : weighting vector
!
!   function(psi) is classified as follows:
!      id = 1  : u * w
!      id = 4  : u'* w
!      id = 8  : u * w'
!      id = 11 : u'* w'
!      id = 2  : a * u * w
!      id = 6  : a'* u * w
!      id = 5  : a * u'* w
!      id = 7  : a'* u'* w
!      id = 9  : a * u * w'
!      id = 13 : a'* u * w'
!      id = 12 : a * u'* w'
!      id = 14 : a'* u'* w'
!
!      id = -1 : a * w
!
!   where ' means the derivative of psi
!
!  < input >
!     id        : mode select
!     a(1:2)  : coefficient vector, optional
!     da(1:2) : derivative a(nnode), optional
!     b(1:2)  : coefficient vector, optional
!     db(1:2) : derivative a(nnode), optional
!  < output >
!     x(4,4)    : matrix of integrated values
!
!-------------------------------------------------------
    integer, intent(in) :: id
    real(8), intent(in), dimension(1:2), optional :: a, da, b, db
    integer :: i, j
    real(8) :: x(4,4), a1, a2, da1, da2, b1, b2, db1, db2

    if(table_initialize_flag.eq.0) then
       call table_initialize
       table_initialize_flag=1
    endif
    
    if(present(a).and.present(da)) then
       a1  = a(1)  ; a2  = a(2)
       da1 = da(1) ; da2 = da(2)
       if(present(b).and.present(db)) then
          b1  = b(1)  ; b2  = b(2)
          db1 = db(1) ; db2 = db(2)
       end if
    end if

    select case(id)
    case(-1)
       x(1,1) = table_hh(1,1) *  a1
       x(2,1) = table_hh(2,1) * da1
       x(3,1) = table_hh(3,1) *  a2
       x(4,1) = table_hh(4,1) * da2
       x(1,2) = table_hh(1,2) *  a1
       x(2,2) = table_hh(2,2) * da1
       x(3,2) = table_hh(3,2) *  a2
       x(4,2) = table_hh(4,2) * da2
       x(1,3) = table_hh(1,3) *  a1
       x(2,3) = table_hh(2,3) * da1
       x(3,3) = table_hh(3,3) *  a2
       x(4,3) = table_hh(4,3) * da2
       x(1,4) = table_hh(1,4) *  a1
       x(2,4) = table_hh(2,4) * da1
       x(3,4) = table_hh(3,4) *  a2
       x(4,4) = table_hh(4,4) * da2
    case(1)
       do i = 1, 4
          do j = 1, 4
             x(j,i) = table_hh(j,i)
          end do
       end do
    case(4)
       do i = 1, 4
          do j = 1, 4
             x(j,i) = table_hh(j+4,i)
          end do
       end do
    case(8)
       do i = 1, 4
          do j = 1, 4
             x(j,i) = table_hh(j,i+4)
          end do
       end do
    case(11)
       do i = 1, 4
          do j = 1, 4
             x(j,i) = table_hh(j+4,i+4)
          end do
       end do
    case(2)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(1,j,i) &
                  & + da1 * table_hhh(2,j,i) &
                  & +  a2 * table_hhh(3,j,i) &
                  & + da2 * table_hhh(4,j,i)
          end do
       end do
    case(6)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(5,j,i) &
                  & + da1 * table_hhh(6,j,i) &
                  & +  a2 * table_hhh(7,j,i) &
                  & + da2 * table_hhh(8,j,i)
          end do
       end do
    case(5)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(1,j+4,i) &
                  & + da1 * table_hhh(2,j+4,i) &
                  & +  a2 * table_hhh(3,j+4,i) &
                  & + da2 * table_hhh(4,j+4,i)
          end do
       end do
    case(7)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(5,j+4,i) &
                  & + da1 * table_hhh(6,j+4,i) &
                  & +  a2 * table_hhh(7,j+4,i) &
                  & + da2 * table_hhh(8,j+4,i)
          end do
       end do
    case(9)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(1,j,i+4) &
                  & + da1 * table_hhh(2,j,i+4) &
                  & +  a2 * table_hhh(3,j,i+4) &
                  & + da2 * table_hhh(4,j,i+4)
          end do
       end do
    case(13)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(5,j,i+4) &
                  & + da1 * table_hhh(6,j,i+4) &
                  & +  a2 * table_hhh(7,j,i+4) &
                  & + da2 * table_hhh(8,j,i+4) 
          end do
       end do
    case(12)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(1,j+4,i+4) &
                  & + da1 * table_hhh(2,j+4,i+4) &
                  & +  a2 * table_hhh(3,j+4,i+4) &
                  & + da2 * table_hhh(4,j+4,i+4)
          end do
       end do
    case(14)
       do i = 1, 4
          do j = 1, 4
             x(j,i) =  a1 * table_hhh(5,j+4,i+4) &
                  & + da1 * table_hhh(6,j+4,i+4) &
                  & +  a2 * table_hhh(7,j+4,i+4) &
                  & + da2 * table_hhh(8,j+4,i+4) 
          end do
       end do
    case default
       stop 'XX Unknown ID in fem_integrate'
    end select

  end subroutine fem_integrate

end module libfem_mod
