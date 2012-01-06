MODULE libfem
  IMPLICIT NONE
  PRIVATE
  real(8), dimension(1:4,1:4), save, public :: table_ll
  real(8), dimension(1:4,1:4), save, public :: table_lg
  real(8), dimension(1:8,1:8), save, public :: table_hh
  real(8), dimension(1:8,1:8), save, public :: table_hg
  real(8), dimension(1:8,1:8), save, public :: table_gg
  real(8), dimension(1:4,1:4,1:4), save, public :: table_lll
  real(8), dimension(1:4,1:4,1:4), save, public :: table_lgl
  real(8), dimension(1:4,1:4,1:4), save, public :: table_llg
  real(8), dimension(1:4,1:4,1:4), save, public :: table_lgg
  real(8), dimension(1:8,1:8,1:8), save, public :: table_hhh
  real(8), dimension(1:8,1:8,1:8), save, public :: table_hgh
  real(8), dimension(1:8,1:8,1:8), save, public :: table_hhg
  real(8), dimension(1:8,1:8,1:8), save, public :: table_hgg
  real(8), dimension(1:8,1:6), save, public :: table_hq
  real(8), dimension(1:6,1:6), save, public :: table_qq
  real(8), dimension(1:8,1:6,1:8), save, public :: table_hqh
  real(8), dimension(1:8,1:8,1:6), save, public :: table_hhq
  real(8), dimension(1:8,1:6,1:6), save, public :: table_hqq
  real(8), dimension(1:4,1:6), save, public :: table_lq
  real(8), dimension(1:6,1:4), save, public :: table_ql
  real(8), dimension(1:4,1:6,1:4), save, public :: table_lql
  real(8), dimension(1:4,1:4,1:6), save, public :: table_llq
  real(8), dimension(1:4,1:6,1:6), save, public :: table_lqq
  integer, save, public :: table_initialize_flag = 0
  public :: table_initialize, fem_integrate

CONTAINS

  subroutine table_initialize

!  Integral value table

!  Relation
!     l     : l1=1-x             dl1=-1
!           : l2=x               dl2= 1
!     h     : h1=(1-x)^2(1+2x)   dh1=-6x(1-x)
!           : h2= x(1-x)^2       dh2=(1-x)(1-3x)
!           : h3= x^2(3-2x)      dh3= 6x(1-x)
!           : h4=-x^2(1-x)       dh4=x(-2+3x)
!     g     : g1=1               dg1=0
!           : g2=x-1/2           dg2=1
!           : g3=(x-1/2)^2/2     dg3=(x-1/2)
!           : g4=(x-1/2)^3/6     dg4=(x-1/2)^2/2
!     q     : q1=(1-2x)(1-x)     dq1=-3+4x
!           : q2=4x(1-x)         dq2=4-8x
!           : q3=x(2x-1)         dq3=4x-1
!

    implicit none
    integer:: k,i,j

    !  *** double product ll ***

      table_ll(1,1)= 1.d0/3.d0
      table_ll(1,2)= 1.d0/6.d0
      table_ll(1,3)= -1.d0/2.d0
      table_ll(1,4)= 1.d0/2.d0

      table_ll(2,1)= 1.d0/6.d0
      table_ll(2,2)= 1.d0/3.d0
      table_ll(2,3)= -1.d0/2.d0
      table_ll(2,4)= 1.d0/2.d0

      table_ll(3,1)= -1.d0/2.d0
      table_ll(3,2)= -1.d0/2.d0
      table_ll(3,3)= 1.d0
      table_ll(3,4)= -1.d0

      table_ll(4,1)= 1.d0/2.d0
      table_ll(4,2)= 1.d0/2.d0
      table_ll(4,3)= -1.d0
      table_ll(4,4)= 1.d0

    !  *** double product lg ***

      table_lg(1,1)= 1.d0/2.d0
      table_lg(1,2)= 0.d0
      table_lg(1,3)= -1.d0/12.d0
      table_lg(1,4)= 1.d0/2.d0

      table_lg(2,1)= 1.d0/2.d0
      table_lg(2,2)= 0.d0
      table_lg(2,3)= 1.d0/12.d0
      table_lg(2,4)= 1.d0/2.d0

      table_lg(3,1)= -1.d0
      table_lg(3,2)= 0.d0
      table_lg(3,3)= 0.d0
      table_lg(3,4)= -1.d0

      table_lg(4,1)= 1.d0
      table_lg(4,2)= 0.d0
      table_lg(4,3)= 0.d0
      table_lg(4,4)= 1.d0

    !  *** double product hh ***

      table_hh(1,1)= 13.d0/35.d0
      table_hh(1,2)= 11.d0/210.d0
      table_hh(1,3)= 9.d0/70.d0
      table_hh(1,4)= -13.d0/420.d0
      table_hh(1,5)= -1.d0/2.d0
      table_hh(1,6)= 1.d0/10.d0
      table_hh(1,7)= 1.d0/2.d0
      table_hh(1,8)= -1.d0/10.d0

      table_hh(2,1)= 11.d0/210.d0
      table_hh(2,2)= 1.d0/105.d0
      table_hh(2,3)= 13.d0/420.d0
      table_hh(2,4)= -1.d0/140.d0
      table_hh(2,5)= -1.d0/10.d0
      table_hh(2,6)= 0.d0
      table_hh(2,7)= 1.d0/10.d0
      table_hh(2,8)= -1.d0/60.d0

      table_hh(3,1)= 9.d0/70.d0
      table_hh(3,2)= 13.d0/420.d0
      table_hh(3,3)= 13.d0/35.d0
      table_hh(3,4)= -11.d0/210.d0
      table_hh(3,5)= -1.d0/2.d0
      table_hh(3,6)= -1.d0/10.d0
      table_hh(3,7)= 1.d0/2.d0
      table_hh(3,8)= 1.d0/10.d0

      table_hh(4,1)= -13.d0/420.d0
      table_hh(4,2)= -1.d0/140.d0
      table_hh(4,3)= -11.d0/210.d0
      table_hh(4,4)= 1.d0/105.d0
      table_hh(4,5)= 1.d0/10.d0
      table_hh(4,6)= 1.d0/60.d0
      table_hh(4,7)= -1.d0/10.d0
      table_hh(4,8)= 0.d0

      table_hh(5,1)= -1.d0/2.d0
      table_hh(5,2)= -1.d0/10.d0
      table_hh(5,3)= -1.d0/2.d0
      table_hh(5,4)= 1.d0/10.d0
      table_hh(5,5)= 6.d0/5.d0
      table_hh(5,6)= 1.d0/10.d0
      table_hh(5,7)= -6.d0/5.d0
      table_hh(5,8)= 1.d0/10.d0

      table_hh(6,1)= 1.d0/10.d0
      table_hh(6,2)= 0.d0
      table_hh(6,3)= -1.d0/10.d0
      table_hh(6,4)= 1.d0/60.d0
      table_hh(6,5)= 1.d0/10.d0
      table_hh(6,6)= 2.d0/15.d0
      table_hh(6,7)= -1.d0/10.d0
      table_hh(6,8)= -1.d0/30.d0

      table_hh(7,1)= 1.d0/2.d0
      table_hh(7,2)= 1.d0/10.d0
      table_hh(7,3)= 1.d0/2.d0
      table_hh(7,4)= -1.d0/10.d0
      table_hh(7,5)= -6.d0/5.d0
      table_hh(7,6)= -1.d0/10.d0
      table_hh(7,7)= 6.d0/5.d0
      table_hh(7,8)= -1.d0/10.d0

      table_hh(8,1)= -1.d0/10.d0
      table_hh(8,2)= -1.d0/60.d0
      table_hh(8,3)= 1.d0/10.d0
      table_hh(8,4)= 0.d0
      table_hh(8,5)= 1.d0/10.d0
      table_hh(8,6)= -1.d0/30.d0
      table_hh(8,7)= -1.d0/10.d0
      table_hh(8,8)= 2.d0/15.d0

    !  *** double product hg ***

      table_hg(1,1)= 1.d0/2.d0
      table_hg(1,2)= 0.d0
      table_hg(1,3)= -1.d0/10.d0
      table_hg(1,4)= 1.d0/2.d0
      table_hg(1,5)= 1.d0/48.d0
      table_hg(1,6)= -1.d0/10.d0
      table_hg(1,7)= -1.d0/420.d0
      table_hg(1,8)= 1.d0/48.d0

      table_hg(2,1)= 1.d0/12.d0
      table_hg(2,2)= 0.d0
      table_hg(2,3)= -1.d0/120.d0
      table_hg(2,4)= 1.d0/12.d0
      table_hg(2,5)= 1.d0/480.d0
      table_hg(2,6)= -1.d0/120.d0
      table_hg(2,7)= -1.d0/6720.d0
      table_hg(2,8)= 1.d0/480.d0

      table_hg(3,1)= 1.d0/2.d0
      table_hg(3,2)= 0.d0
      table_hg(3,3)= 1.d0/10.d0
      table_hg(3,4)= 1.d0/2.d0
      table_hg(3,5)= 1.d0/48.d0
      table_hg(3,6)= 1.d0/10.d0
      table_hg(3,7)= 1.d0/420.d0
      table_hg(3,8)= 1.d0/48.d0

      table_hg(4,1)= -1.d0/12.d0
      table_hg(4,2)= 0.d0
      table_hg(4,3)= -1.d0/120.d0
      table_hg(4,4)= -1.d0/12.d0
      table_hg(4,5)= -1.d0/480.d0
      table_hg(4,6)= -1.d0/120.d0
      table_hg(4,7)= -1.d0/6720.d0
      table_hg(4,8)= -1.d0/480.d0

      table_hg(5,1)= -1.d0
      table_hg(5,2)= 0.d0
      table_hg(5,3)= 0.d0
      table_hg(5,4)= -1.d0
      table_hg(5,5)= -1.d0/40.d0
      table_hg(5,6)= 0.d0
      table_hg(5,7)= 0.d0
      table_hg(5,8)= -1.d0/40.d0

      table_hg(6,1)= 0.d0
      table_hg(6,2)= 0.d0
      table_hg(6,3)= -1.d0/12.d0
      table_hg(6,4)= 0.d0
      table_hg(6,5)= 1.d0/120.d0
      table_hg(6,6)= -1.d0/12.d0
      table_hg(6,7)= -1.d0/480.d0
      table_hg(6,8)= 1.d0/120.d0

      table_hg(7,1)= 1.d0
      table_hg(7,2)= 0.d0
      table_hg(7,3)= 0.d0
      table_hg(7,4)= 1.d0
      table_hg(7,5)= 1.d0/40.d0
      table_hg(7,6)= 0.d0
      table_hg(7,7)= 0.d0
      table_hg(7,8)= 1.d0/40.d0

      table_hg(8,1)= 0.d0
      table_hg(8,2)= 0.d0
      table_hg(8,3)= 1.d0/12.d0
      table_hg(8,4)= 0.d0
      table_hg(8,5)= 1.d0/120.d0
      table_hg(8,6)= 1.d0/12.d0
      table_hg(8,7)= 1.d0/480.d0
      table_hg(8,8)= 1.d0/120.d0

    !  *** double product gg ***

      table_gg(1,1)= 1.d0
      table_gg(1,2)= 0.d0
      table_gg(1,3)= 0.d0
      table_gg(1,4)= 1.d0
      table_gg(1,5)= 1.d0/24.d0
      table_gg(1,6)= 0.d0
      table_gg(1,7)= 0.d0
      table_gg(1,8)= 1.d0/24.d0

      table_gg(2,1)= 0.d0
      table_gg(2,2)= 0.d0
      table_gg(2,3)= 0.d0
      table_gg(2,4)= 0.d0
      table_gg(2,5)= 0.d0
      table_gg(2,6)= 0.d0
      table_gg(2,7)= 0.d0
      table_gg(2,8)= 0.d0

      table_gg(3,1)= 0.d0
      table_gg(3,2)= 0.d0
      table_gg(3,3)= 1.d0/12.d0
      table_gg(3,4)= 0.d0
      table_gg(3,5)= 0.d0
      table_gg(3,6)= 1.d0/12.d0
      table_gg(3,7)= 1.d0/480.d0
      table_gg(3,8)= 0.d0

      table_gg(4,1)= 1.d0
      table_gg(4,2)= 0.d0
      table_gg(4,3)= 0.d0
      table_gg(4,4)= 1.d0
      table_gg(4,5)= 1.d0/24.d0
      table_gg(4,6)= 0.d0
      table_gg(4,7)= 0.d0
      table_gg(4,8)= 1.d0/24.d0

      table_gg(5,1)= 1.d0/24.d0
      table_gg(5,2)= 0.d0
      table_gg(5,3)= 0.d0
      table_gg(5,4)= 1.d0/24.d0
      table_gg(5,5)= 1.d0/320.d0
      table_gg(5,6)= 0.d0
      table_gg(5,7)= 0.d0
      table_gg(5,8)= 1.d0/320.d0

      table_gg(6,1)= 0.d0
      table_gg(6,2)= 0.d0
      table_gg(6,3)= 1.d0/12.d0
      table_gg(6,4)= 0.d0
      table_gg(6,5)= 0.d0
      table_gg(6,6)= 1.d0/12.d0
      table_gg(6,7)= 1.d0/480.d0
      table_gg(6,8)= 0.d0

      table_gg(7,1)= 0.d0
      table_gg(7,2)= 0.d0
      table_gg(7,3)= 1.d0/480.d0
      table_gg(7,4)= 0.d0
      table_gg(7,5)= 0.d0
      table_gg(7,6)= 1.d0/480.d0
      table_gg(7,7)= 1.d0/16128.d0
      table_gg(7,8)= 0.d0

      table_gg(8,1)= 1.d0/24.d0
      table_gg(8,2)= 0.d0
      table_gg(8,3)= 0.d0
      table_gg(8,4)= 1.d0/24.d0
      table_gg(8,5)= 1.d0/320.d0
      table_gg(8,6)= 0.d0
      table_gg(8,7)= 0.d0
      table_gg(8,8)= 1.d0/320.d0

    !  *** triple product lll ***

      table_lll(1,1,1)= 1.d0/4.d0
      table_lll(1,1,2)= 1.d0/12.d0
      table_lll(1,1,3)= -1.d0/3.d0
      table_lll(1,1,4)= 1.d0/3.d0

      table_lll(1,2,1)= 1.d0/12.d0
      table_lll(1,2,2)= 1.d0/12.d0
      table_lll(1,2,3)= -1.d0/6.d0
      table_lll(1,2,4)= 1.d0/6.d0

      table_lll(1,3,1)= -1.d0/3.d0
      table_lll(1,3,2)= -1.d0/6.d0
      table_lll(1,3,3)= 1.d0/2.d0
      table_lll(1,3,4)= -1.d0/2.d0

      table_lll(1,4,1)= 1.d0/3.d0
      table_lll(1,4,2)= 1.d0/6.d0
      table_lll(1,4,3)= -1.d0/2.d0
      table_lll(1,4,4)= 1.d0/2.d0


      table_lll(2,1,1)= 1.d0/12.d0
      table_lll(2,1,2)= 1.d0/12.d0
      table_lll(2,1,3)= -1.d0/6.d0
      table_lll(2,1,4)= 1.d0/6.d0

      table_lll(2,2,1)= 1.d0/12.d0
      table_lll(2,2,2)= 1.d0/4.d0
      table_lll(2,2,3)= -1.d0/3.d0
      table_lll(2,2,4)= 1.d0/3.d0

      table_lll(2,3,1)= -1.d0/6.d0
      table_lll(2,3,2)= -1.d0/3.d0
      table_lll(2,3,3)= 1.d0/2.d0
      table_lll(2,3,4)= -1.d0/2.d0

      table_lll(2,4,1)= 1.d0/6.d0
      table_lll(2,4,2)= 1.d0/3.d0
      table_lll(2,4,3)= -1.d0/2.d0
      table_lll(2,4,4)= 1.d0/2.d0


      table_lll(3,1,1)= -1.d0/3.d0
      table_lll(3,1,2)= -1.d0/6.d0
      table_lll(3,1,3)= 1.d0/2.d0
      table_lll(3,1,4)= -1.d0/2.d0

      table_lll(3,2,1)= -1.d0/6.d0
      table_lll(3,2,2)= -1.d0/3.d0
      table_lll(3,2,3)= 1.d0/2.d0
      table_lll(3,2,4)= -1.d0/2.d0

      table_lll(3,3,1)= 1.d0/2.d0
      table_lll(3,3,2)= 1.d0/2.d0
      table_lll(3,3,3)= -1.d0
      table_lll(3,3,4)= 1.d0

      table_lll(3,4,1)= -1.d0/2.d0
      table_lll(3,4,2)= -1.d0/2.d0
      table_lll(3,4,3)= 1.d0
      table_lll(3,4,4)= -1.d0


      table_lll(4,1,1)= 1.d0/3.d0
      table_lll(4,1,2)= 1.d0/6.d0
      table_lll(4,1,3)= -1.d0/2.d0
      table_lll(4,1,4)= 1.d0/2.d0

      table_lll(4,2,1)= 1.d0/6.d0
      table_lll(4,2,2)= 1.d0/3.d0
      table_lll(4,2,3)= -1.d0/2.d0
      table_lll(4,2,4)= 1.d0/2.d0

      table_lll(4,3,1)= -1.d0/2.d0
      table_lll(4,3,2)= -1.d0/2.d0
      table_lll(4,3,3)= 1.d0
      table_lll(4,3,4)= -1.d0

      table_lll(4,4,1)= 1.d0/2.d0
      table_lll(4,4,2)= 1.d0/2.d0
      table_lll(4,4,3)= -1.d0
      table_lll(4,4,4)= 1.d0

    !  *** triple product lgl ***

      table_lgl(1,1,1)= 1.d0/3.d0
      table_lgl(1,1,2)= 1.d0/6.d0
      table_lgl(1,1,3)= -1.d0/2.d0
      table_lgl(1,1,4)= 1.d0/2.d0

      table_lgl(1,2,1)= 0.d0
      table_lgl(1,2,2)= 0.d0
      table_lgl(1,2,3)= 0.d0
      table_lgl(1,2,4)= 0.d0

      table_lgl(1,3,1)= -1.d0/12.d0
      table_lgl(1,3,2)= 0.d0
      table_lgl(1,3,3)= 1.d0/12.d0
      table_lgl(1,3,4)= -1.d0/12.d0

      table_lgl(1,4,1)= 1.d0/3.d0
      table_lgl(1,4,2)= 1.d0/6.d0
      table_lgl(1,4,3)= -1.d0/2.d0
      table_lgl(1,4,4)= 1.d0/2.d0


      table_lgl(2,1,1)= 1.d0/6.d0
      table_lgl(2,1,2)= 1.d0/3.d0
      table_lgl(2,1,3)= -1.d0/2.d0
      table_lgl(2,1,4)= 1.d0/2.d0

      table_lgl(2,2,1)= 0.d0
      table_lgl(2,2,2)= 0.d0
      table_lgl(2,2,3)= 0.d0
      table_lgl(2,2,4)= 0.d0

      table_lgl(2,3,1)= 0.d0
      table_lgl(2,3,2)= 1.d0/12.d0
      table_lgl(2,3,3)= -1.d0/12.d0
      table_lgl(2,3,4)= 1.d0/12.d0

      table_lgl(2,4,1)= 1.d0/6.d0
      table_lgl(2,4,2)= 1.d0/3.d0
      table_lgl(2,4,3)= -1.d0/2.d0
      table_lgl(2,4,4)= 1.d0/2.d0


      table_lgl(3,1,1)= -1.d0/2.d0
      table_lgl(3,1,2)= -1.d0/2.d0
      table_lgl(3,1,3)= 1.d0
      table_lgl(3,1,4)= -1.d0

      table_lgl(3,2,1)= 0.d0
      table_lgl(3,2,2)= 0.d0
      table_lgl(3,2,3)= 0.d0
      table_lgl(3,2,4)= 0.d0

      table_lgl(3,3,1)= 1.d0/12.d0
      table_lgl(3,3,2)= -1.d0/12.d0
      table_lgl(3,3,3)= 0.d0
      table_lgl(3,3,4)= 0.d0

      table_lgl(3,4,1)= -1.d0/2.d0
      table_lgl(3,4,2)= -1.d0/2.d0
      table_lgl(3,4,3)= 1.d0
      table_lgl(3,4,4)= -1.d0


      table_lgl(4,1,1)= 1.d0/2.d0
      table_lgl(4,1,2)= 1.d0/2.d0
      table_lgl(4,1,3)= -1.d0
      table_lgl(4,1,4)= 1.d0

      table_lgl(4,2,1)= 0.d0
      table_lgl(4,2,2)= 0.d0
      table_lgl(4,2,3)= 0.d0
      table_lgl(4,2,4)= 0.d0

      table_lgl(4,3,1)= -1.d0/12.d0
      table_lgl(4,3,2)= 1.d0/12.d0
      table_lgl(4,3,3)= 0.d0
      table_lgl(4,3,4)= 0.d0

      table_lgl(4,4,1)= 1.d0/2.d0
      table_lgl(4,4,2)= 1.d0/2.d0
      table_lgl(4,4,3)= -1.d0
      table_lgl(4,4,4)= 1.d0

    !  *** triple product hhg ***

      do k=1,4
         do j=1,4
            do i=1,4
               table_llg(i,k,j)=table_lgl(i,j,k)
            enddo
         enddo
      enddo

    !  *** triple product lgg ***

      table_lgg(1,1,1)= 1.d0/2.d0
      table_lgg(1,1,2)= 0.d0
      table_lgg(1,1,3)= -1.d0/12.d0
      table_lgg(1,1,4)= 1.d0/2.d0

      table_lgg(1,2,1)= 0.d0
      table_lgg(1,2,2)= 0.d0
      table_lgg(1,2,3)= 0.d0
      table_lgg(1,2,4)= 0.d0

      table_lgg(1,3,1)= -1.d0/12.d0
      table_lgg(1,3,2)= 0.d0
      table_lgg(1,3,3)= 1.d0/24.d0
      table_lgg(1,3,4)= -1.d0/12.d0

      table_lgg(1,4,1)= 1.d0/2.d0
      table_lgg(1,4,2)= 0.d0
      table_lgg(1,4,3)= -1.d0/12.d0
      table_lgg(1,4,4)= 1.d0/2.d0


      table_lgg(2,1,1)= 1.d0/2.d0
      table_lgg(2,1,2)= 0.d0
      table_lgg(2,1,3)= 1.d0/12.d0
      table_lgg(2,1,4)= 1.d0/2.d0

      table_lgg(2,2,1)= 0.d0
      table_lgg(2,2,2)= 0.d0
      table_lgg(2,2,3)= 0.d0
      table_lgg(2,2,4)= 0.d0

      table_lgg(2,3,1)= 1.d0/12.d0
      table_lgg(2,3,2)= 0.d0
      table_lgg(2,3,3)= 1.d0/24.d0
      table_lgg(2,3,4)= 1.d0/12.d0

      table_lgg(2,4,1)= 1.d0/2.d0
      table_lgg(2,4,2)= 0.d0
      table_lgg(2,4,3)= 1.d0/12.d0
      table_lgg(2,4,4)= 1.d0/2.d0


      table_lgg(3,1,1)= -1.d0
      table_lgg(3,1,2)= 0.d0
      table_lgg(3,1,3)= 0.d0
      table_lgg(3,1,4)= -1.d0

      table_lgg(3,2,1)= 0.d0
      table_lgg(3,2,2)= 0.d0
      table_lgg(3,2,3)= 0.d0
      table_lgg(3,2,4)= 0.d0

      table_lgg(3,3,1)= 0.d0
      table_lgg(3,3,2)= 0.d0
      table_lgg(3,3,3)= -1.d0/12.d0
      table_lgg(3,3,4)= 0.d0

      table_lgg(3,4,1)= -1.d0
      table_lgg(3,4,2)= 0.d0
      table_lgg(3,4,3)= 0.d0
      table_lgg(3,4,4)= -1.d0


      table_lgg(4,1,1)= 1.d0
      table_lgg(4,1,2)= 0.d0
      table_lgg(4,1,3)= 0.d0
      table_lgg(4,1,4)= 1.d0

      table_lgg(4,2,1)= 0.d0
      table_lgg(4,2,2)= 0.d0
      table_lgg(4,2,3)= 0.d0
      table_lgg(4,2,4)= 0.d0

      table_lgg(4,3,1)= 0.d0
      table_lgg(4,3,2)= 0.d0
      table_lgg(4,3,3)= 1.d0/12.d0
      table_lgg(4,3,4)= 0.d0

      table_lgg(4,4,1)= 1.d0
      table_lgg(4,4,2)= 0.d0
      table_lgg(4,4,3)= 0.d0
      table_lgg(4,4,4)= 1.d0

    !  *** triple product hhh ***

      table_hhh(1,1,1)= 43.d0/140.d0
      table_hhh(1,1,2)= 97.d0/2520.d0
      table_hhh(1,1,3)= 9.d0/140.d0
      table_hhh(1,1,4)= -43.d0/2520.d0
      table_hhh(1,1,5)= -1.d0/3.d0
      table_hhh(1,1,6)= 5.d0/42.d0
      table_hhh(1,1,7)= 1.d0/3.d0
      table_hhh(1,1,8)= -17.d0/210.d0

      table_hhh(1,2,1)= 97.d0/2520.d0
      table_hhh(1,2,2)= 2.d0/315.d0
      table_hhh(1,2,3)= 1.d0/72.d0
      table_hhh(1,2,4)= -1.d0/280.d0
      table_hhh(1,2,5)= -5.d0/84.d0
      table_hhh(1,2,6)= 1.d0/168.d0
      table_hhh(1,2,7)= 5.d0/84.d0
      table_hhh(1,2,8)= -11.d0/840.d0

      table_hhh(1,3,1)= 9.d0/140.d0
      table_hhh(1,3,2)= 1.d0/72.d0
      table_hhh(1,3,3)= 9.d0/140.d0
      table_hhh(1,3,4)= -1.d0/72.d0
      table_hhh(1,3,5)= -1.d0/6.d0
      table_hhh(1,3,6)= -2.d0/105.d0
      table_hhh(1,3,7)= 1.d0/6.d0
      table_hhh(1,3,8)= -2.d0/105.d0

      table_hhh(1,4,1)= -43.d0/2520.d0
      table_hhh(1,4,2)= -1.d0/280.d0
      table_hhh(1,4,3)= -1.d0/72.d0
      table_hhh(1,4,4)= 1.d0/315.d0
      table_hhh(1,4,5)= 17.d0/420.d0
      table_hhh(1,4,6)= 1.d0/280.d0
      table_hhh(1,4,7)= -17.d0/420.d0
      table_hhh(1,4,8)= 1.d0/168.d0

      table_hhh(1,5,1)= -1.d0/3.d0
      table_hhh(1,5,2)= -5.d0/84.d0
      table_hhh(1,5,3)= -1.d0/6.d0
      table_hhh(1,5,4)= 17.d0/420.d0
      table_hhh(1,5,5)= 3.d0/5.d0
      table_hhh(1,5,6)= -1.d0/70.d0
      table_hhh(1,5,7)= -3.d0/5.d0
      table_hhh(1,5,8)= 4.d0/35.d0

      table_hhh(1,6,1)= 5.d0/42.d0
      table_hhh(1,6,2)= 1.d0/168.d0
      table_hhh(1,6,3)= -2.d0/105.d0
      table_hhh(1,6,4)= 1.d0/280.d0
      table_hhh(1,6,5)= -1.d0/70.d0
      table_hhh(1,6,6)= 43.d0/420.d0
      table_hhh(1,6,7)= 1.d0/70.d0
      table_hhh(1,6,8)= -1.d0/60.d0

      table_hhh(1,7,1)= 1.d0/3.d0
      table_hhh(1,7,2)= 5.d0/84.d0
      table_hhh(1,7,3)= 1.d0/6.d0
      table_hhh(1,7,4)= -17.d0/420.d0
      table_hhh(1,7,5)= -3.d0/5.d0
      table_hhh(1,7,6)= 1.d0/70.d0
      table_hhh(1,7,7)= 3.d0/5.d0
      table_hhh(1,7,8)= -4.d0/35.d0

      table_hhh(1,8,1)= -17.d0/210.d0
      table_hhh(1,8,2)= -11.d0/840.d0
      table_hhh(1,8,3)= -2.d0/105.d0
      table_hhh(1,8,4)= 1.d0/168.d0
      table_hhh(1,8,5)= 4.d0/35.d0
      table_hhh(1,8,6)= -1.d0/60.d0
      table_hhh(1,8,7)= -4.d0/35.d0
      table_hhh(1,8,8)= 13.d0/420.d0


      table_hhh(2,1,1)= 97.d0/2520.d0
      table_hhh(2,1,2)= 2.d0/315.d0
      table_hhh(2,1,3)= 1.d0/72.d0
      table_hhh(2,1,4)= -1.d0/280.d0
      table_hhh(2,1,5)= -5.d0/84.d0
      table_hhh(2,1,6)= 1.d0/168.d0
      table_hhh(2,1,7)= 5.d0/84.d0
      table_hhh(2,1,8)= -11.d0/840.d0

      table_hhh(2,2,1)= 2.d0/315.d0
      table_hhh(2,2,2)= 1.d0/840.d0
      table_hhh(2,2,3)= 1.d0/315.d0
      table_hhh(2,2,4)= -1.d0/1260.d0
      table_hhh(2,2,5)= -1.d0/84.d0
      table_hhh(2,2,6)= 0.d0
      table_hhh(2,2,7)= 1.d0/84.d0
      table_hhh(2,2,8)= -1.d0/420.d0

      table_hhh(2,3,1)= 1.d0/72.d0
      table_hhh(2,3,2)= 1.d0/315.d0
      table_hhh(2,3,3)= 43.d0/2520.d0
      table_hhh(2,3,4)= -1.d0/280.d0
      table_hhh(2,3,5)= -17.d0/420.d0
      table_hhh(2,3,6)= -1.d0/168.d0
      table_hhh(2,3,7)= 17.d0/420.d0
      table_hhh(2,3,8)= -1.d0/280.d0

      table_hhh(2,4,1)= -1.d0/280.d0
      table_hhh(2,4,2)= -1.d0/1260.d0
      table_hhh(2,4,3)= -1.d0/280.d0
      table_hhh(2,4,4)= 1.d0/1260.d0
      table_hhh(2,4,5)= 1.d0/105.d0
      table_hhh(2,4,6)= 1.d0/840.d0
      table_hhh(2,4,7)= -1.d0/105.d0
      table_hhh(2,4,8)= 1.d0/840.d0

      table_hhh(2,5,1)= -5.d0/84.d0
      table_hhh(2,5,2)= -1.d0/84.d0
      table_hhh(2,5,3)= -17.d0/420.d0
      table_hhh(2,5,4)= 1.d0/105.d0
      table_hhh(2,5,5)= 9.d0/70.d0
      table_hhh(2,5,6)= 1.d0/140.d0
      table_hhh(2,5,7)= -9.d0/70.d0
      table_hhh(2,5,8)= 3.d0/140.d0

      table_hhh(2,6,1)= 1.d0/168.d0
      table_hhh(2,6,2)= 0.d0
      table_hhh(2,6,3)= -1.d0/168.d0
      table_hhh(2,6,4)= 1.d0/840.d0
      table_hhh(2,6,5)= 1.d0/140.d0
      table_hhh(2,6,6)= 1.d0/120.d0
      table_hhh(2,6,7)= -1.d0/140.d0
      table_hhh(2,6,8)= -1.d0/840.d0

      table_hhh(2,7,1)= 5.d0/84.d0
      table_hhh(2,7,2)= 1.d0/84.d0
      table_hhh(2,7,3)= 17.d0/420.d0
      table_hhh(2,7,4)= -1.d0/105.d0
      table_hhh(2,7,5)= -9.d0/70.d0
      table_hhh(2,7,6)= -1.d0/140.d0
      table_hhh(2,7,7)= 9.d0/70.d0
      table_hhh(2,7,8)= -3.d0/140.d0

      table_hhh(2,8,1)= -11.d0/840.d0
      table_hhh(2,8,2)= -1.d0/420.d0
      table_hhh(2,8,3)= -1.d0/280.d0
      table_hhh(2,8,4)= 1.d0/840.d0
      table_hhh(2,8,5)= 3.d0/140.d0
      table_hhh(2,8,6)= -1.d0/840.d0
      table_hhh(2,8,7)= -3.d0/140.d0
      table_hhh(2,8,8)= 1.d0/168.d0


      table_hhh(3,1,1)= 9.d0/140.d0
      table_hhh(3,1,2)= 1.d0/72.d0
      table_hhh(3,1,3)= 9.d0/140.d0
      table_hhh(3,1,4)= -1.d0/72.d0
      table_hhh(3,1,5)= -1.d0/6.d0
      table_hhh(3,1,6)= -2.d0/105.d0
      table_hhh(3,1,7)= 1.d0/6.d0
      table_hhh(3,1,8)= -2.d0/105.d0

      table_hhh(3,2,1)= 1.d0/72.d0
      table_hhh(3,2,2)= 1.d0/315.d0
      table_hhh(3,2,3)= 43.d0/2520.d0
      table_hhh(3,2,4)= -1.d0/280.d0
      table_hhh(3,2,5)= -17.d0/420.d0
      table_hhh(3,2,6)= -1.d0/168.d0
      table_hhh(3,2,7)= 17.d0/420.d0
      table_hhh(3,2,8)= -1.d0/280.d0

      table_hhh(3,3,1)= 9.d0/140.d0
      table_hhh(3,3,2)= 43.d0/2520.d0
      table_hhh(3,3,3)= 43.d0/140.d0
      table_hhh(3,3,4)= -97.d0/2520.d0
      table_hhh(3,3,5)= -1.d0/3.d0
      table_hhh(3,3,6)= -17.d0/210.d0
      table_hhh(3,3,7)= 1.d0/3.d0
      table_hhh(3,3,8)= 5.d0/42.d0

      table_hhh(3,4,1)= -1.d0/72.d0
      table_hhh(3,4,2)= -1.d0/280.d0
      table_hhh(3,4,3)= -97.d0/2520.d0
      table_hhh(3,4,4)= 2.d0/315.d0
      table_hhh(3,4,5)= 5.d0/84.d0
      table_hhh(3,4,6)= 11.d0/840.d0
      table_hhh(3,4,7)= -5.d0/84.d0
      table_hhh(3,4,8)= -1.d0/168.d0

      table_hhh(3,5,1)= -1.d0/6.d0
      table_hhh(3,5,2)= -17.d0/420.d0
      table_hhh(3,5,3)= -1.d0/3.d0
      table_hhh(3,5,4)= 5.d0/84.d0
      table_hhh(3,5,5)= 3.d0/5.d0
      table_hhh(3,5,6)= 4.d0/35.d0
      table_hhh(3,5,7)= -3.d0/5.d0
      table_hhh(3,5,8)= -1.d0/70.d0

      table_hhh(3,6,1)= -2.d0/105.d0
      table_hhh(3,6,2)= -1.d0/168.d0
      table_hhh(3,6,3)= -17.d0/210.d0
      table_hhh(3,6,4)= 11.d0/840.d0
      table_hhh(3,6,5)= 4.d0/35.d0
      table_hhh(3,6,6)= 13.d0/420.d0
      table_hhh(3,6,7)= -4.d0/35.d0
      table_hhh(3,6,8)= -1.d0/60.d0

      table_hhh(3,7,1)= 1.d0/6.d0
      table_hhh(3,7,2)= 17.d0/420.d0
      table_hhh(3,7,3)= 1.d0/3.d0
      table_hhh(3,7,4)= -5.d0/84.d0
      table_hhh(3,7,5)= -3.d0/5.d0
      table_hhh(3,7,6)= -4.d0/35.d0
      table_hhh(3,7,7)= 3.d0/5.d0
      table_hhh(3,7,8)= 1.d0/70.d0

      table_hhh(3,8,1)= -2.d0/105.d0
      table_hhh(3,8,2)= -1.d0/280.d0
      table_hhh(3,8,3)= 5.d0/42.d0
      table_hhh(3,8,4)= -1.d0/168.d0
      table_hhh(3,8,5)= -1.d0/70.d0
      table_hhh(3,8,6)= -1.d0/60.d0
      table_hhh(3,8,7)= 1.d0/70.d0
      table_hhh(3,8,8)= 43.d0/420.d0


      table_hhh(4,1,1)= -43.d0/2520.d0
      table_hhh(4,1,2)= -1.d0/280.d0
      table_hhh(4,1,3)= -1.d0/72.d0
      table_hhh(4,1,4)= 1.d0/315.d0
      table_hhh(4,1,5)= 17.d0/420.d0
      table_hhh(4,1,6)= 1.d0/280.d0
      table_hhh(4,1,7)= -17.d0/420.d0
      table_hhh(4,1,8)= 1.d0/168.d0

      table_hhh(4,2,1)= -1.d0/280.d0
      table_hhh(4,2,2)= -1.d0/1260.d0
      table_hhh(4,2,3)= -1.d0/280.d0
      table_hhh(4,2,4)= 1.d0/1260.d0
      table_hhh(4,2,5)= 1.d0/105.d0
      table_hhh(4,2,6)= 1.d0/840.d0
      table_hhh(4,2,7)= -1.d0/105.d0
      table_hhh(4,2,8)= 1.d0/840.d0

      table_hhh(4,3,1)= -1.d0/72.d0
      table_hhh(4,3,2)= -1.d0/280.d0
      table_hhh(4,3,3)= -97.d0/2520.d0
      table_hhh(4,3,4)= 2.d0/315.d0
      table_hhh(4,3,5)= 5.d0/84.d0
      table_hhh(4,3,6)= 11.d0/840.d0
      table_hhh(4,3,7)= -5.d0/84.d0
      table_hhh(4,3,8)= -1.d0/168.d0

      table_hhh(4,4,1)= 1.d0/315.d0
      table_hhh(4,4,2)= 1.d0/1260.d0
      table_hhh(4,4,3)= 2.d0/315.d0
      table_hhh(4,4,4)= -1.d0/840.d0
      table_hhh(4,4,5)= -1.d0/84.d0
      table_hhh(4,4,6)= -1.d0/420.d0
      table_hhh(4,4,7)= 1.d0/84.d0
      table_hhh(4,4,8)= 0.d0

      table_hhh(4,5,1)= 17.d0/420.d0
      table_hhh(4,5,2)= 1.d0/105.d0
      table_hhh(4,5,3)= 5.d0/84.d0
      table_hhh(4,5,4)= -1.d0/84.d0
      table_hhh(4,5,5)= -9.d0/70.d0
      table_hhh(4,5,6)= -3.d0/140.d0
      table_hhh(4,5,7)= 9.d0/70.d0
      table_hhh(4,5,8)= -1.d0/140.d0

      table_hhh(4,6,1)= 1.d0/280.d0
      table_hhh(4,6,2)= 1.d0/840.d0
      table_hhh(4,6,3)= 11.d0/840.d0
      table_hhh(4,6,4)= -1.d0/420.d0
      table_hhh(4,6,5)= -3.d0/140.d0
      table_hhh(4,6,6)= -1.d0/168.d0
      table_hhh(4,6,7)= 3.d0/140.d0
      table_hhh(4,6,8)= 1.d0/840.d0

      table_hhh(4,7,1)= -17.d0/420.d0
      table_hhh(4,7,2)= -1.d0/105.d0
      table_hhh(4,7,3)= -5.d0/84.d0
      table_hhh(4,7,4)= 1.d0/84.d0
      table_hhh(4,7,5)= 9.d0/70.d0
      table_hhh(4,7,6)= 3.d0/140.d0
      table_hhh(4,7,7)= -9.d0/70.d0
      table_hhh(4,7,8)= 1.d0/140.d0

      table_hhh(4,8,1)= 1.d0/168.d0
      table_hhh(4,8,2)= 1.d0/840.d0
      table_hhh(4,8,3)= -1.d0/168.d0
      table_hhh(4,8,4)= 0.d0
      table_hhh(4,8,5)= -1.d0/140.d0
      table_hhh(4,8,6)= 1.d0/840.d0
      table_hhh(4,8,7)= 1.d0/140.d0
      table_hhh(4,8,8)= -1.d0/120.d0


      table_hhh(5,1,1)= -1.d0/3.d0
      table_hhh(5,1,2)= -5.d0/84.d0
      table_hhh(5,1,3)= -1.d0/6.d0
      table_hhh(5,1,4)= 17.d0/420.d0
      table_hhh(5,1,5)= 3.d0/5.d0
      table_hhh(5,1,6)= -1.d0/70.d0
      table_hhh(5,1,7)= -3.d0/5.d0
      table_hhh(5,1,8)= 4.d0/35.d0

      table_hhh(5,2,1)= -5.d0/84.d0
      table_hhh(5,2,2)= -1.d0/84.d0
      table_hhh(5,2,3)= -17.d0/420.d0
      table_hhh(5,2,4)= 1.d0/105.d0
      table_hhh(5,2,5)= 9.d0/70.d0
      table_hhh(5,2,6)= 1.d0/140.d0
      table_hhh(5,2,7)= -9.d0/70.d0
      table_hhh(5,2,8)= 3.d0/140.d0

      table_hhh(5,3,1)= -1.d0/6.d0
      table_hhh(5,3,2)= -17.d0/420.d0
      table_hhh(5,3,3)= -1.d0/3.d0
      table_hhh(5,3,4)= 5.d0/84.d0
      table_hhh(5,3,5)= 3.d0/5.d0
      table_hhh(5,3,6)= 4.d0/35.d0
      table_hhh(5,3,7)= -3.d0/5.d0
      table_hhh(5,3,8)= -1.d0/70.d0

      table_hhh(5,4,1)= 17.d0/420.d0
      table_hhh(5,4,2)= 1.d0/105.d0
      table_hhh(5,4,3)= 5.d0/84.d0
      table_hhh(5,4,4)= -1.d0/84.d0
      table_hhh(5,4,5)= -9.d0/70.d0
      table_hhh(5,4,6)= -3.d0/140.d0
      table_hhh(5,4,7)= 9.d0/70.d0
      table_hhh(5,4,8)= -1.d0/140.d0

      table_hhh(5,5,1)= 3.d0/5.d0
      table_hhh(5,5,2)= 9.d0/70.d0
      table_hhh(5,5,3)= 3.d0/5.d0
      table_hhh(5,5,4)= -9.d0/70.d0
      table_hhh(5,5,5)= -54.d0/35.d0
      table_hhh(5,5,6)= -6.d0/35.d0
      table_hhh(5,5,7)= 54.d0/35.d0
      table_hhh(5,5,8)= -6.d0/35.d0

      table_hhh(5,6,1)= -1.d0/70.d0
      table_hhh(5,6,2)= 1.d0/140.d0
      table_hhh(5,6,3)= 4.d0/35.d0
      table_hhh(5,6,4)= -3.d0/140.d0
      table_hhh(5,6,5)= -6.d0/35.d0
      table_hhh(5,6,6)= -3.d0/35.d0
      table_hhh(5,6,7)= 6.d0/35.d0
      table_hhh(5,6,8)= 1.d0/70.d0

      table_hhh(5,7,1)= -3.d0/5.d0
      table_hhh(5,7,2)= -9.d0/70.d0
      table_hhh(5,7,3)= -3.d0/5.d0
      table_hhh(5,7,4)= 9.d0/70.d0
      table_hhh(5,7,5)= 54.d0/35.d0
      table_hhh(5,7,6)= 6.d0/35.d0
      table_hhh(5,7,7)= -54.d0/35.d0
      table_hhh(5,7,8)= 6.d0/35.d0

      table_hhh(5,8,1)= 4.d0/35.d0
      table_hhh(5,8,2)= 3.d0/140.d0
      table_hhh(5,8,3)= -1.d0/70.d0
      table_hhh(5,8,4)= -1.d0/140.d0
      table_hhh(5,8,5)= -6.d0/35.d0
      table_hhh(5,8,6)= 1.d0/70.d0
      table_hhh(5,8,7)= 6.d0/35.d0
      table_hhh(5,8,8)= -3.d0/35.d0


      table_hhh(6,1,1)= 5.d0/42.d0
      table_hhh(6,1,2)= 1.d0/168.d0
      table_hhh(6,1,3)= -2.d0/105.d0
      table_hhh(6,1,4)= 1.d0/280.d0
      table_hhh(6,1,5)= -1.d0/70.d0
      table_hhh(6,1,6)= 43.d0/420.d0
      table_hhh(6,1,7)= 1.d0/70.d0
      table_hhh(6,1,8)= -1.d0/60.d0

      table_hhh(6,2,1)= 1.d0/168.d0
      table_hhh(6,2,2)= 0.d0
      table_hhh(6,2,3)= -1.d0/168.d0
      table_hhh(6,2,4)= 1.d0/840.d0
      table_hhh(6,2,5)= 1.d0/140.d0
      table_hhh(6,2,6)= 1.d0/120.d0
      table_hhh(6,2,7)= -1.d0/140.d0
      table_hhh(6,2,8)= -1.d0/840.d0

      table_hhh(6,3,1)= -2.d0/105.d0
      table_hhh(6,3,2)= -1.d0/168.d0
      table_hhh(6,3,3)= -17.d0/210.d0
      table_hhh(6,3,4)= 11.d0/840.d0
      table_hhh(6,3,5)= 4.d0/35.d0
      table_hhh(6,3,6)= 13.d0/420.d0
      table_hhh(6,3,7)= -4.d0/35.d0
      table_hhh(6,3,8)= -1.d0/60.d0

      table_hhh(6,4,1)= 1.d0/280.d0
      table_hhh(6,4,2)= 1.d0/840.d0
      table_hhh(6,4,3)= 11.d0/840.d0
      table_hhh(6,4,4)= -1.d0/420.d0
      table_hhh(6,4,5)= -3.d0/140.d0
      table_hhh(6,4,6)= -1.d0/168.d0
      table_hhh(6,4,7)= 3.d0/140.d0
      table_hhh(6,4,8)= 1.d0/840.d0

      table_hhh(6,5,1)= -1.d0/70.d0
      table_hhh(6,5,2)= 1.d0/140.d0
      table_hhh(6,5,3)= 4.d0/35.d0
      table_hhh(6,5,4)= -3.d0/140.d0
      table_hhh(6,5,5)= -6.d0/35.d0
      table_hhh(6,5,6)= -3.d0/35.d0
      table_hhh(6,5,7)= 6.d0/35.d0
      table_hhh(6,5,8)= 1.d0/70.d0

      table_hhh(6,6,1)= 43.d0/420.d0
      table_hhh(6,6,2)= 1.d0/120.d0
      table_hhh(6,6,3)= 13.d0/420.d0
      table_hhh(6,6,4)= -1.d0/168.d0
      table_hhh(6,6,5)= -3.d0/35.d0
      table_hhh(6,6,6)= 2.d0/35.d0
      table_hhh(6,6,7)= 3.d0/35.d0
      table_hhh(6,6,8)= -1.d0/105.d0

      table_hhh(6,7,1)= 1.d0/70.d0
      table_hhh(6,7,2)= -1.d0/140.d0
      table_hhh(6,7,3)= -4.d0/35.d0
      table_hhh(6,7,4)= 3.d0/140.d0
      table_hhh(6,7,5)= 6.d0/35.d0
      table_hhh(6,7,6)= 3.d0/35.d0
      table_hhh(6,7,7)= -6.d0/35.d0
      table_hhh(6,7,8)= -1.d0/70.d0

      table_hhh(6,8,1)= -1.d0/60.d0
      table_hhh(6,8,2)= -1.d0/840.d0
      table_hhh(6,8,3)= -1.d0/60.d0
      table_hhh(6,8,4)= 1.d0/840.d0
      table_hhh(6,8,5)= 1.d0/70.d0
      table_hhh(6,8,6)= -1.d0/105.d0
      table_hhh(6,8,7)= -1.d0/70.d0
      table_hhh(6,8,8)= -1.d0/105.d0


      table_hhh(7,1,1)= 1.d0/3.d0
      table_hhh(7,1,2)= 5.d0/84.d0
      table_hhh(7,1,3)= 1.d0/6.d0
      table_hhh(7,1,4)= -17.d0/420.d0
      table_hhh(7,1,5)= -3.d0/5.d0
      table_hhh(7,1,6)= 1.d0/70.d0
      table_hhh(7,1,7)= 3.d0/5.d0
      table_hhh(7,1,8)= -4.d0/35.d0

      table_hhh(7,2,1)= 5.d0/84.d0
      table_hhh(7,2,2)= 1.d0/84.d0
      table_hhh(7,2,3)= 17.d0/420.d0
      table_hhh(7,2,4)= -1.d0/105.d0
      table_hhh(7,2,5)= -9.d0/70.d0
      table_hhh(7,2,6)= -1.d0/140.d0
      table_hhh(7,2,7)= 9.d0/70.d0
      table_hhh(7,2,8)= -3.d0/140.d0

      table_hhh(7,3,1)= 1.d0/6.d0
      table_hhh(7,3,2)= 17.d0/420.d0
      table_hhh(7,3,3)= 1.d0/3.d0
      table_hhh(7,3,4)= -5.d0/84.d0
      table_hhh(7,3,5)= -3.d0/5.d0
      table_hhh(7,3,6)= -4.d0/35.d0
      table_hhh(7,3,7)= 3.d0/5.d0
      table_hhh(7,3,8)= 1.d0/70.d0

      table_hhh(7,4,1)= -17.d0/420.d0
      table_hhh(7,4,2)= -1.d0/105.d0
      table_hhh(7,4,3)= -5.d0/84.d0
      table_hhh(7,4,4)= 1.d0/84.d0
      table_hhh(7,4,5)= 9.d0/70.d0
      table_hhh(7,4,6)= 3.d0/140.d0
      table_hhh(7,4,7)= -9.d0/70.d0
      table_hhh(7,4,8)= 1.d0/140.d0

      table_hhh(7,5,1)= -3.d0/5.d0
      table_hhh(7,5,2)= -9.d0/70.d0
      table_hhh(7,5,3)= -3.d0/5.d0
      table_hhh(7,5,4)= 9.d0/70.d0
      table_hhh(7,5,5)= 54.d0/35.d0
      table_hhh(7,5,6)= 6.d0/35.d0
      table_hhh(7,5,7)= -54.d0/35.d0
      table_hhh(7,5,8)= 6.d0/35.d0

      table_hhh(7,6,1)= 1.d0/70.d0
      table_hhh(7,6,2)= -1.d0/140.d0
      table_hhh(7,6,3)= -4.d0/35.d0
      table_hhh(7,6,4)= 3.d0/140.d0
      table_hhh(7,6,5)= 6.d0/35.d0
      table_hhh(7,6,6)= 3.d0/35.d0
      table_hhh(7,6,7)= -6.d0/35.d0
      table_hhh(7,6,8)= -1.d0/70.d0

      table_hhh(7,7,1)= 3.d0/5.d0
      table_hhh(7,7,2)= 9.d0/70.d0
      table_hhh(7,7,3)= 3.d0/5.d0
      table_hhh(7,7,4)= -9.d0/70.d0
      table_hhh(7,7,5)= -54.d0/35.d0
      table_hhh(7,7,6)= -6.d0/35.d0
      table_hhh(7,7,7)= 54.d0/35.d0
      table_hhh(7,7,8)= -6.d0/35.d0

      table_hhh(7,8,1)= -4.d0/35.d0
      table_hhh(7,8,2)= -3.d0/140.d0
      table_hhh(7,8,3)= 1.d0/70.d0
      table_hhh(7,8,4)= 1.d0/140.d0
      table_hhh(7,8,5)= 6.d0/35.d0
      table_hhh(7,8,6)= -1.d0/70.d0
      table_hhh(7,8,7)= -6.d0/35.d0
      table_hhh(7,8,8)= 3.d0/35.d0


      table_hhh(8,1,1)= -17.d0/210.d0
      table_hhh(8,1,2)= -11.d0/840.d0
      table_hhh(8,1,3)= -2.d0/105.d0
      table_hhh(8,1,4)= 1.d0/168.d0
      table_hhh(8,1,5)= 4.d0/35.d0
      table_hhh(8,1,6)= -1.d0/60.d0
      table_hhh(8,1,7)= -4.d0/35.d0
      table_hhh(8,1,8)= 13.d0/420.d0

      table_hhh(8,2,1)= -11.d0/840.d0
      table_hhh(8,2,2)= -1.d0/420.d0
      table_hhh(8,2,3)= -1.d0/280.d0
      table_hhh(8,2,4)= 1.d0/840.d0
      table_hhh(8,2,5)= 3.d0/140.d0
      table_hhh(8,2,6)= -1.d0/840.d0
      table_hhh(8,2,7)= -3.d0/140.d0
      table_hhh(8,2,8)= 1.d0/168.d0

      table_hhh(8,3,1)= -2.d0/105.d0
      table_hhh(8,3,2)= -1.d0/280.d0
      table_hhh(8,3,3)= 5.d0/42.d0
      table_hhh(8,3,4)= -1.d0/168.d0
      table_hhh(8,3,5)= -1.d0/70.d0
      table_hhh(8,3,6)= -1.d0/60.d0
      table_hhh(8,3,7)= 1.d0/70.d0
      table_hhh(8,3,8)= 43.d0/420.d0

      table_hhh(8,4,1)= 1.d0/168.d0
      table_hhh(8,4,2)= 1.d0/840.d0
      table_hhh(8,4,3)= -1.d0/168.d0
      table_hhh(8,4,4)= 0.d0
      table_hhh(8,4,5)= -1.d0/140.d0
      table_hhh(8,4,6)= 1.d0/840.d0
      table_hhh(8,4,7)= 1.d0/140.d0
      table_hhh(8,4,8)= -1.d0/120.d0

      table_hhh(8,5,1)= 4.d0/35.d0
      table_hhh(8,5,2)= 3.d0/140.d0
      table_hhh(8,5,3)= -1.d0/70.d0
      table_hhh(8,5,4)= -1.d0/140.d0
      table_hhh(8,5,5)= -6.d0/35.d0
      table_hhh(8,5,6)= 1.d0/70.d0
      table_hhh(8,5,7)= 6.d0/35.d0
      table_hhh(8,5,8)= -3.d0/35.d0

      table_hhh(8,6,1)= -1.d0/60.d0
      table_hhh(8,6,2)= -1.d0/840.d0
      table_hhh(8,6,3)= -1.d0/60.d0
      table_hhh(8,6,4)= 1.d0/840.d0
      table_hhh(8,6,5)= 1.d0/70.d0
      table_hhh(8,6,6)= -1.d0/105.d0
      table_hhh(8,6,7)= -1.d0/70.d0
      table_hhh(8,6,8)= -1.d0/105.d0

      table_hhh(8,7,1)= -4.d0/35.d0
      table_hhh(8,7,2)= -3.d0/140.d0
      table_hhh(8,7,3)= 1.d0/70.d0
      table_hhh(8,7,4)= 1.d0/140.d0
      table_hhh(8,7,5)= 6.d0/35.d0
      table_hhh(8,7,6)= -1.d0/70.d0
      table_hhh(8,7,7)= -6.d0/35.d0
      table_hhh(8,7,8)= 3.d0/35.d0

      table_hhh(8,8,1)= 13.d0/420.d0
      table_hhh(8,8,2)= 1.d0/168.d0
      table_hhh(8,8,3)= 43.d0/420.d0
      table_hhh(8,8,4)= -1.d0/120.d0
      table_hhh(8,8,5)= -3.d0/35.d0
      table_hhh(8,8,6)= -1.d0/105.d0
      table_hhh(8,8,7)= 3.d0/35.d0
      table_hhh(8,8,8)= 2.d0/35.d0

    !  *** triple product hgh ***

      table_hgh(1,1,1)= 13.d0/35.d0
      table_hgh(1,1,2)= 11.d0/210.d0
      table_hgh(1,1,3)= 9.d0/70.d0
      table_hgh(1,1,4)= -13.d0/420.d0
      table_hgh(1,1,5)= -1.d0/2.d0
      table_hgh(1,1,6)= 1.d0/10.d0
      table_hgh(1,1,7)= 1.d0/2.d0
      table_hgh(1,1,8)= -1.d0/10.d0

      table_hgh(1,2,1)= 0.d0
      table_hgh(1,2,2)= 0.d0
      table_hgh(1,2,3)= 0.d0
      table_hgh(1,2,4)= 0.d0
      table_hgh(1,2,5)= 0.d0
      table_hgh(1,2,6)= 0.d0
      table_hgh(1,2,7)= 0.d0
      table_hgh(1,2,8)= 0.d0

      table_hgh(1,3,1)= -1.d0/10.d0
      table_hgh(1,3,2)= -1.d0/105.d0
      table_hgh(1,3,3)= 0.d0
      table_hgh(1,3,4)= 1.d0/840.d0
      table_hgh(1,3,5)= 9.d0/140.d0
      table_hgh(1,3,6)= -5.d0/84.d0
      table_hgh(1,3,7)= -9.d0/140.d0
      table_hgh(1,3,8)= 1.d0/42.d0

      table_hgh(1,4,1)= 13.d0/35.d0
      table_hgh(1,4,2)= 11.d0/210.d0
      table_hgh(1,4,3)= 9.d0/70.d0
      table_hgh(1,4,4)= -13.d0/420.d0
      table_hgh(1,4,5)= -1.d0/2.d0
      table_hgh(1,4,6)= 1.d0/10.d0
      table_hgh(1,4,7)= 1.d0/2.d0
      table_hgh(1,4,8)= -1.d0/10.d0

      table_hgh(1,5,1)= 47.d0/2520.d0
      table_hgh(1,5,2)= 1.d0/630.d0
      table_hgh(1,5,3)= 11.d0/5040.d0
      table_hgh(1,5,4)= -1.d0/2016.d0
      table_hgh(1,5,5)= -1.d0/80.d0
      table_hgh(1,5,6)= 19.d0/1680.d0
      table_hgh(1,5,7)= 1.d0/80.d0
      table_hgh(1,5,8)= -1.d0/336.d0

      table_hgh(1,6,1)= -1.d0/10.d0
      table_hgh(1,6,2)= -1.d0/105.d0
      table_hgh(1,6,3)= 0.d0
      table_hgh(1,6,4)= 1.d0/840.d0
      table_hgh(1,6,5)= 9.d0/140.d0
      table_hgh(1,6,6)= -5.d0/84.d0
      table_hgh(1,6,7)= -9.d0/140.d0
      table_hgh(1,6,8)= 1.d0/42.d0

      table_hgh(1,7,1)= -1.d0/420.d0
      table_hgh(1,7,2)= -1.d0/6048.d0
      table_hgh(1,7,3)= 0.d0
      table_hgh(1,7,4)= 1.d0/60480.d0
      table_hgh(1,7,5)= 11.d0/10080.d0
      table_hgh(1,7,6)= -17.d0/10080.d0
      table_hgh(1,7,7)= -11.d0/10080.d0
      table_hgh(1,7,8)= 1.d0/2520.d0

      table_hgh(1,8,1)= 47.d0/2520.d0
      table_hgh(1,8,2)= 1.d0/630.d0
      table_hgh(1,8,3)= 11.d0/5040.d0
      table_hgh(1,8,4)= -1.d0/2016.d0
      table_hgh(1,8,5)= -1.d0/80.d0
      table_hgh(1,8,6)= 19.d0/1680.d0
      table_hgh(1,8,7)= 1.d0/80.d0
      table_hgh(1,8,8)= -1.d0/336.d0


      table_hgh(2,1,1)= 11.d0/210.d0
      table_hgh(2,1,2)= 1.d0/105.d0
      table_hgh(2,1,3)= 13.d0/420.d0
      table_hgh(2,1,4)= -1.d0/140.d0
      table_hgh(2,1,5)= -1.d0/10.d0
      table_hgh(2,1,6)= 0.d0
      table_hgh(2,1,7)= 1.d0/10.d0
      table_hgh(2,1,8)= -1.d0/60.d0

      table_hgh(2,2,1)= 0.d0
      table_hgh(2,2,2)= 0.d0
      table_hgh(2,2,3)= 0.d0
      table_hgh(2,2,4)= 0.d0
      table_hgh(2,2,5)= 0.d0
      table_hgh(2,2,6)= 0.d0
      table_hgh(2,2,7)= 0.d0
      table_hgh(2,2,8)= 0.d0

      table_hgh(2,3,1)= -1.d0/105.d0
      table_hgh(2,3,2)= -1.d0/840.d0
      table_hgh(2,3,3)= 1.d0/840.d0
      table_hgh(2,3,4)= 0.d0
      table_hgh(2,3,5)= 1.d0/140.d0
      table_hgh(2,3,6)= -1.d0/210.d0
      table_hgh(2,3,7)= -1.d0/140.d0
      table_hgh(2,3,8)= 1.d0/280.d0

      table_hgh(2,4,1)= 11.d0/210.d0
      table_hgh(2,4,2)= 1.d0/105.d0
      table_hgh(2,4,3)= 13.d0/420.d0
      table_hgh(2,4,4)= -1.d0/140.d0
      table_hgh(2,4,5)= -1.d0/10.d0
      table_hgh(2,4,6)= 0.d0
      table_hgh(2,4,7)= 1.d0/10.d0
      table_hgh(2,4,8)= -1.d0/60.d0

      table_hgh(2,5,1)= 1.d0/630.d0
      table_hgh(2,5,2)= 1.d0/5040.d0
      table_hgh(2,5,3)= 1.d0/2016.d0
      table_hgh(2,5,4)= -1.d0/10080.d0
      table_hgh(2,5,5)= -1.d0/560.d0
      table_hgh(2,5,6)= 1.d0/1680.d0
      table_hgh(2,5,7)= 1.d0/560.d0
      table_hgh(2,5,8)= -1.d0/3360.d0

      table_hgh(2,6,1)= -1.d0/105.d0
      table_hgh(2,6,2)= -1.d0/840.d0
      table_hgh(2,6,3)= 1.d0/840.d0
      table_hgh(2,6,4)= 0.d0
      table_hgh(2,6,5)= 1.d0/140.d0
      table_hgh(2,6,6)= -1.d0/210.d0
      table_hgh(2,6,7)= -1.d0/140.d0
      table_hgh(2,6,8)= 1.d0/280.d0

      table_hgh(2,7,1)= -1.d0/6048.d0
      table_hgh(2,7,2)= -1.d0/60480.d0
      table_hgh(2,7,3)= 1.d0/60480.d0
      table_hgh(2,7,4)= 0.d0
      table_hgh(2,7,5)= 1.d0/10080.d0
      table_hgh(2,7,6)= -1.d0/10080.d0
      table_hgh(2,7,7)= -1.d0/10080.d0
      table_hgh(2,7,8)= 1.d0/20160.d0

      table_hgh(2,8,1)= 1.d0/630.d0
      table_hgh(2,8,2)= 1.d0/5040.d0
      table_hgh(2,8,3)= 1.d0/2016.d0
      table_hgh(2,8,4)= -1.d0/10080.d0
      table_hgh(2,8,5)= -1.d0/560.d0
      table_hgh(2,8,6)= 1.d0/1680.d0
      table_hgh(2,8,7)= 1.d0/560.d0
      table_hgh(2,8,8)= -1.d0/3360.d0


      table_hgh(3,1,1)= 9.d0/70.d0
      table_hgh(3,1,2)= 13.d0/420.d0
      table_hgh(3,1,3)= 13.d0/35.d0
      table_hgh(3,1,4)= -11.d0/210.d0
      table_hgh(3,1,5)= -1.d0/2.d0
      table_hgh(3,1,6)= -1.d0/10.d0
      table_hgh(3,1,7)= 1.d0/2.d0
      table_hgh(3,1,8)= 1.d0/10.d0

      table_hgh(3,2,1)= 0.d0
      table_hgh(3,2,2)= 0.d0
      table_hgh(3,2,3)= 0.d0
      table_hgh(3,2,4)= 0.d0
      table_hgh(3,2,5)= 0.d0
      table_hgh(3,2,6)= 0.d0
      table_hgh(3,2,7)= 0.d0
      table_hgh(3,2,8)= 0.d0

      table_hgh(3,3,1)= 0.d0
      table_hgh(3,3,2)= 1.d0/840.d0
      table_hgh(3,3,3)= 1.d0/10.d0
      table_hgh(3,3,4)= -1.d0/105.d0
      table_hgh(3,3,5)= -9.d0/140.d0
      table_hgh(3,3,6)= -1.d0/42.d0
      table_hgh(3,3,7)= 9.d0/140.d0
      table_hgh(3,3,8)= 5.d0/84.d0

      table_hgh(3,4,1)= 9.d0/70.d0
      table_hgh(3,4,2)= 13.d0/420.d0
      table_hgh(3,4,3)= 13.d0/35.d0
      table_hgh(3,4,4)= -11.d0/210.d0
      table_hgh(3,4,5)= -1.d0/2.d0
      table_hgh(3,4,6)= -1.d0/10.d0
      table_hgh(3,4,7)= 1.d0/2.d0
      table_hgh(3,4,8)= 1.d0/10.d0

      table_hgh(3,5,1)= 11.d0/5040.d0
      table_hgh(3,5,2)= 1.d0/2016.d0
      table_hgh(3,5,3)= 47.d0/2520.d0
      table_hgh(3,5,4)= -1.d0/630.d0
      table_hgh(3,5,5)= -1.d0/80.d0
      table_hgh(3,5,6)= -1.d0/336.d0
      table_hgh(3,5,7)= 1.d0/80.d0
      table_hgh(3,5,8)= 19.d0/1680.d0

      table_hgh(3,6,1)= 0.d0
      table_hgh(3,6,2)= 1.d0/840.d0
      table_hgh(3,6,3)= 1.d0/10.d0
      table_hgh(3,6,4)= -1.d0/105.d0
      table_hgh(3,6,5)= -9.d0/140.d0
      table_hgh(3,6,6)= -1.d0/42.d0
      table_hgh(3,6,7)= 9.d0/140.d0
      table_hgh(3,6,8)= 5.d0/84.d0

      table_hgh(3,7,1)= 0.d0
      table_hgh(3,7,2)= 1.d0/60480.d0
      table_hgh(3,7,3)= 1.d0/420.d0
      table_hgh(3,7,4)= -1.d0/6048.d0
      table_hgh(3,7,5)= -11.d0/10080.d0
      table_hgh(3,7,6)= -1.d0/2520.d0
      table_hgh(3,7,7)= 11.d0/10080.d0
      table_hgh(3,7,8)= 17.d0/10080.d0

      table_hgh(3,8,1)= 11.d0/5040.d0
      table_hgh(3,8,2)= 1.d0/2016.d0
      table_hgh(3,8,3)= 47.d0/2520.d0
      table_hgh(3,8,4)= -1.d0/630.d0
      table_hgh(3,8,5)= -1.d0/80.d0
      table_hgh(3,8,6)= -1.d0/336.d0
      table_hgh(3,8,7)= 1.d0/80.d0
      table_hgh(3,8,8)= 19.d0/1680.d0


      table_hgh(4,1,1)= -13.d0/420.d0
      table_hgh(4,1,2)= -1.d0/140.d0
      table_hgh(4,1,3)= -11.d0/210.d0
      table_hgh(4,1,4)= 1.d0/105.d0
      table_hgh(4,1,5)= 1.d0/10.d0
      table_hgh(4,1,6)= 1.d0/60.d0
      table_hgh(4,1,7)= -1.d0/10.d0
      table_hgh(4,1,8)= 0.d0

      table_hgh(4,2,1)= 0.d0
      table_hgh(4,2,2)= 0.d0
      table_hgh(4,2,3)= 0.d0
      table_hgh(4,2,4)= 0.d0
      table_hgh(4,2,5)= 0.d0
      table_hgh(4,2,6)= 0.d0
      table_hgh(4,2,7)= 0.d0
      table_hgh(4,2,8)= 0.d0

      table_hgh(4,3,1)= 1.d0/840.d0
      table_hgh(4,3,2)= 0.d0
      table_hgh(4,3,3)= -1.d0/105.d0
      table_hgh(4,3,4)= 1.d0/840.d0
      table_hgh(4,3,5)= 1.d0/140.d0
      table_hgh(4,3,6)= 1.d0/280.d0
      table_hgh(4,3,7)= -1.d0/140.d0
      table_hgh(4,3,8)= -1.d0/210.d0

      table_hgh(4,4,1)= -13.d0/420.d0
      table_hgh(4,4,2)= -1.d0/140.d0
      table_hgh(4,4,3)= -11.d0/210.d0
      table_hgh(4,4,4)= 1.d0/105.d0
      table_hgh(4,4,5)= 1.d0/10.d0
      table_hgh(4,4,6)= 1.d0/60.d0
      table_hgh(4,4,7)= -1.d0/10.d0
      table_hgh(4,4,8)= 0.d0

      table_hgh(4,5,1)= -1.d0/2016.d0
      table_hgh(4,5,2)= -1.d0/10080.d0
      table_hgh(4,5,3)= -1.d0/630.d0
      table_hgh(4,5,4)= 1.d0/5040.d0
      table_hgh(4,5,5)= 1.d0/560.d0
      table_hgh(4,5,6)= 1.d0/3360.d0
      table_hgh(4,5,7)= -1.d0/560.d0
      table_hgh(4,5,8)= -1.d0/1680.d0

      table_hgh(4,6,1)= 1.d0/840.d0
      table_hgh(4,6,2)= 0.d0
      table_hgh(4,6,3)= -1.d0/105.d0
      table_hgh(4,6,4)= 1.d0/840.d0
      table_hgh(4,6,5)= 1.d0/140.d0
      table_hgh(4,6,6)= 1.d0/280.d0
      table_hgh(4,6,7)= -1.d0/140.d0
      table_hgh(4,6,8)= -1.d0/210.d0

      table_hgh(4,7,1)= 1.d0/60480.d0
      table_hgh(4,7,2)= 0.d0
      table_hgh(4,7,3)= -1.d0/6048.d0
      table_hgh(4,7,4)= 1.d0/60480.d0
      table_hgh(4,7,5)= 1.d0/10080.d0
      table_hgh(4,7,6)= 1.d0/20160.d0
      table_hgh(4,7,7)= -1.d0/10080.d0
      table_hgh(4,7,8)= -1.d0/10080.d0

      table_hgh(4,8,1)= -1.d0/2016.d0
      table_hgh(4,8,2)= -1.d0/10080.d0
      table_hgh(4,8,3)= -1.d0/630.d0
      table_hgh(4,8,4)= 1.d0/5040.d0
      table_hgh(4,8,5)= 1.d0/560.d0
      table_hgh(4,8,6)= 1.d0/3360.d0
      table_hgh(4,8,7)= -1.d0/560.d0
      table_hgh(4,8,8)= -1.d0/1680.d0


      table_hgh(5,1,1)= -1.d0/2.d0
      table_hgh(5,1,2)= -1.d0/10.d0
      table_hgh(5,1,3)= -1.d0/2.d0
      table_hgh(5,1,4)= 1.d0/10.d0
      table_hgh(5,1,5)= 6.d0/5.d0
      table_hgh(5,1,6)= 1.d0/10.d0
      table_hgh(5,1,7)= -6.d0/5.d0
      table_hgh(5,1,8)= 1.d0/10.d0

      table_hgh(5,2,1)= 0.d0
      table_hgh(5,2,2)= 0.d0
      table_hgh(5,2,3)= 0.d0
      table_hgh(5,2,4)= 0.d0
      table_hgh(5,2,5)= 0.d0
      table_hgh(5,2,6)= 0.d0
      table_hgh(5,2,7)= 0.d0
      table_hgh(5,2,8)= 0.d0

      table_hgh(5,3,1)= 9.d0/140.d0
      table_hgh(5,3,2)= 1.d0/140.d0
      table_hgh(5,3,3)= -9.d0/140.d0
      table_hgh(5,3,4)= 1.d0/140.d0
      table_hgh(5,3,5)= 0.d0
      table_hgh(5,3,6)= 1.d0/20.d0
      table_hgh(5,3,7)= 0.d0
      table_hgh(5,3,8)= -1.d0/20.d0

      table_hgh(5,4,1)= -1.d0/2.d0
      table_hgh(5,4,2)= -1.d0/10.d0
      table_hgh(5,4,3)= -1.d0/2.d0
      table_hgh(5,4,4)= 1.d0/10.d0
      table_hgh(5,4,5)= 6.d0/5.d0
      table_hgh(5,4,6)= 1.d0/10.d0
      table_hgh(5,4,7)= -6.d0/5.d0
      table_hgh(5,4,8)= 1.d0/10.d0

      table_hgh(5,5,1)= -1.d0/80.d0
      table_hgh(5,5,2)= -1.d0/560.d0
      table_hgh(5,5,3)= -1.d0/80.d0
      table_hgh(5,5,4)= 1.d0/560.d0
      table_hgh(5,5,5)= 3.d0/140.d0
      table_hgh(5,5,6)= -1.d0/560.d0
      table_hgh(5,5,7)= -3.d0/140.d0
      table_hgh(5,5,8)= -1.d0/560.d0

      table_hgh(5,6,1)= 9.d0/140.d0
      table_hgh(5,6,2)= 1.d0/140.d0
      table_hgh(5,6,3)= -9.d0/140.d0
      table_hgh(5,6,4)= 1.d0/140.d0
      table_hgh(5,6,5)= 0.d0
      table_hgh(5,6,6)= 1.d0/20.d0
      table_hgh(5,6,7)= 0.d0
      table_hgh(5,6,8)= -1.d0/20.d0

      table_hgh(5,7,1)= 11.d0/10080.d0
      table_hgh(5,7,2)= 1.d0/10080.d0
      table_hgh(5,7,3)= -11.d0/10080.d0
      table_hgh(5,7,4)= 1.d0/10080.d0
      table_hgh(5,7,5)= 0.d0
      table_hgh(5,7,6)= 1.d0/1120.d0
      table_hgh(5,7,7)= 0.d0
      table_hgh(5,7,8)= -1.d0/1120.d0

      table_hgh(5,8,1)= -1.d0/80.d0
      table_hgh(5,8,2)= -1.d0/560.d0
      table_hgh(5,8,3)= -1.d0/80.d0
      table_hgh(5,8,4)= 1.d0/560.d0
      table_hgh(5,8,5)= 3.d0/140.d0
      table_hgh(5,8,6)= -1.d0/560.d0
      table_hgh(5,8,7)= -3.d0/140.d0
      table_hgh(5,8,8)= -1.d0/560.d0


      table_hgh(6,1,1)= 1.d0/10.d0
      table_hgh(6,1,2)= 0.d0
      table_hgh(6,1,3)= -1.d0/10.d0
      table_hgh(6,1,4)= 1.d0/60.d0
      table_hgh(6,1,5)= 1.d0/10.d0
      table_hgh(6,1,6)= 2.d0/15.d0
      table_hgh(6,1,7)= -1.d0/10.d0
      table_hgh(6,1,8)= -1.d0/30.d0

      table_hgh(6,2,1)= 0.d0
      table_hgh(6,2,2)= 0.d0
      table_hgh(6,2,3)= 0.d0
      table_hgh(6,2,4)= 0.d0
      table_hgh(6,2,5)= 0.d0
      table_hgh(6,2,6)= 0.d0
      table_hgh(6,2,7)= 0.d0
      table_hgh(6,2,8)= 0.d0

      table_hgh(6,3,1)= -5.d0/84.d0
      table_hgh(6,3,2)= -1.d0/210.d0
      table_hgh(6,3,3)= -1.d0/42.d0
      table_hgh(6,3,4)= 1.d0/280.d0
      table_hgh(6,3,5)= 1.d0/20.d0
      table_hgh(6,3,6)= -1.d0/30.d0
      table_hgh(6,3,7)= -1.d0/20.d0
      table_hgh(6,3,8)= 0.d0

      table_hgh(6,4,1)= 1.d0/10.d0
      table_hgh(6,4,2)= 0.d0
      table_hgh(6,4,3)= -1.d0/10.d0
      table_hgh(6,4,4)= 1.d0/60.d0
      table_hgh(6,4,5)= 1.d0/10.d0
      table_hgh(6,4,6)= 2.d0/15.d0
      table_hgh(6,4,7)= -1.d0/10.d0
      table_hgh(6,4,8)= -1.d0/30.d0

      table_hgh(6,5,1)= 19.d0/1680.d0
      table_hgh(6,5,2)= 1.d0/1680.d0
      table_hgh(6,5,3)= -1.d0/336.d0
      table_hgh(6,5,4)= 1.d0/3360.d0
      table_hgh(6,5,5)= -1.d0/560.d0
      table_hgh(6,5,6)= 1.d0/105.d0
      table_hgh(6,5,7)= 1.d0/560.d0
      table_hgh(6,5,8)= -1.d0/336.d0

      table_hgh(6,6,1)= -5.d0/84.d0
      table_hgh(6,6,2)= -1.d0/210.d0
      table_hgh(6,6,3)= -1.d0/42.d0
      table_hgh(6,6,4)= 1.d0/280.d0
      table_hgh(6,6,5)= 1.d0/20.d0
      table_hgh(6,6,6)= -1.d0/30.d0
      table_hgh(6,6,7)= -1.d0/20.d0
      table_hgh(6,6,8)= 0.d0

      table_hgh(6,7,1)= -17.d0/10080.d0
      table_hgh(6,7,2)= -1.d0/10080.d0
      table_hgh(6,7,3)= -1.d0/2520.d0
      table_hgh(6,7,4)= 1.d0/20160.d0
      table_hgh(6,7,5)= 1.d0/1120.d0
      table_hgh(6,7,6)= -1.d0/840.d0
      table_hgh(6,7,7)= -1.d0/1120.d0
      table_hgh(6,7,8)= 0.d0

      table_hgh(6,8,1)= 19.d0/1680.d0
      table_hgh(6,8,2)= 1.d0/1680.d0
      table_hgh(6,8,3)= -1.d0/336.d0
      table_hgh(6,8,4)= 1.d0/3360.d0
      table_hgh(6,8,5)= -1.d0/560.d0
      table_hgh(6,8,6)= 1.d0/105.d0
      table_hgh(6,8,7)= 1.d0/560.d0
      table_hgh(6,8,8)= -1.d0/336.d0


      table_hgh(7,1,1)= 1.d0/2.d0
      table_hgh(7,1,2)= 1.d0/10.d0
      table_hgh(7,1,3)= 1.d0/2.d0
      table_hgh(7,1,4)= -1.d0/10.d0
      table_hgh(7,1,5)= -6.d0/5.d0
      table_hgh(7,1,6)= -1.d0/10.d0
      table_hgh(7,1,7)= 6.d0/5.d0
      table_hgh(7,1,8)= -1.d0/10.d0

      table_hgh(7,2,1)= 0.d0
      table_hgh(7,2,2)= 0.d0
      table_hgh(7,2,3)= 0.d0
      table_hgh(7,2,4)= 0.d0
      table_hgh(7,2,5)= 0.d0
      table_hgh(7,2,6)= 0.d0
      table_hgh(7,2,7)= 0.d0
      table_hgh(7,2,8)= 0.d0

      table_hgh(7,3,1)= -9.d0/140.d0
      table_hgh(7,3,2)= -1.d0/140.d0
      table_hgh(7,3,3)= 9.d0/140.d0
      table_hgh(7,3,4)= -1.d0/140.d0
      table_hgh(7,3,5)= 0.d0
      table_hgh(7,3,6)= -1.d0/20.d0
      table_hgh(7,3,7)= 0.d0
      table_hgh(7,3,8)= 1.d0/20.d0

      table_hgh(7,4,1)= 1.d0/2.d0
      table_hgh(7,4,2)= 1.d0/10.d0
      table_hgh(7,4,3)= 1.d0/2.d0
      table_hgh(7,4,4)= -1.d0/10.d0
      table_hgh(7,4,5)= -6.d0/5.d0
      table_hgh(7,4,6)= -1.d0/10.d0
      table_hgh(7,4,7)= 6.d0/5.d0
      table_hgh(7,4,8)= -1.d0/10.d0

      table_hgh(7,5,1)= 1.d0/80.d0
      table_hgh(7,5,2)= 1.d0/560.d0
      table_hgh(7,5,3)= 1.d0/80.d0
      table_hgh(7,5,4)= -1.d0/560.d0
      table_hgh(7,5,5)= -3.d0/140.d0
      table_hgh(7,5,6)= 1.d0/560.d0
      table_hgh(7,5,7)= 3.d0/140.d0
      table_hgh(7,5,8)= 1.d0/560.d0

      table_hgh(7,6,1)= -9.d0/140.d0
      table_hgh(7,6,2)= -1.d0/140.d0
      table_hgh(7,6,3)= 9.d0/140.d0
      table_hgh(7,6,4)= -1.d0/140.d0
      table_hgh(7,6,5)= 0.d0
      table_hgh(7,6,6)= -1.d0/20.d0
      table_hgh(7,6,7)= 0.d0
      table_hgh(7,6,8)= 1.d0/20.d0

      table_hgh(7,7,1)= -11.d0/10080.d0
      table_hgh(7,7,2)= -1.d0/10080.d0
      table_hgh(7,7,3)= 11.d0/10080.d0
      table_hgh(7,7,4)= -1.d0/10080.d0
      table_hgh(7,7,5)= 0.d0
      table_hgh(7,7,6)= -1.d0/1120.d0
      table_hgh(7,7,7)= 0.d0
      table_hgh(7,7,8)= 1.d0/1120.d0

      table_hgh(7,8,1)= 1.d0/80.d0
      table_hgh(7,8,2)= 1.d0/560.d0
      table_hgh(7,8,3)= 1.d0/80.d0
      table_hgh(7,8,4)= -1.d0/560.d0
      table_hgh(7,8,5)= -3.d0/140.d0
      table_hgh(7,8,6)= 1.d0/560.d0
      table_hgh(7,8,7)= 3.d0/140.d0
      table_hgh(7,8,8)= 1.d0/560.d0


      table_hgh(8,1,1)= -1.d0/10.d0
      table_hgh(8,1,2)= -1.d0/60.d0
      table_hgh(8,1,3)= 1.d0/10.d0
      table_hgh(8,1,4)= 0.d0
      table_hgh(8,1,5)= 1.d0/10.d0
      table_hgh(8,1,6)= -1.d0/30.d0
      table_hgh(8,1,7)= -1.d0/10.d0
      table_hgh(8,1,8)= 2.d0/15.d0

      table_hgh(8,2,1)= 0.d0
      table_hgh(8,2,2)= 0.d0
      table_hgh(8,2,3)= 0.d0
      table_hgh(8,2,4)= 0.d0
      table_hgh(8,2,5)= 0.d0
      table_hgh(8,2,6)= 0.d0
      table_hgh(8,2,7)= 0.d0
      table_hgh(8,2,8)= 0.d0

      table_hgh(8,3,1)= 1.d0/42.d0
      table_hgh(8,3,2)= 1.d0/280.d0
      table_hgh(8,3,3)= 5.d0/84.d0
      table_hgh(8,3,4)= -1.d0/210.d0
      table_hgh(8,3,5)= -1.d0/20.d0
      table_hgh(8,3,6)= 0.d0
      table_hgh(8,3,7)= 1.d0/20.d0
      table_hgh(8,3,8)= 1.d0/30.d0

      table_hgh(8,4,1)= -1.d0/10.d0
      table_hgh(8,4,2)= -1.d0/60.d0
      table_hgh(8,4,3)= 1.d0/10.d0
      table_hgh(8,4,4)= 0.d0
      table_hgh(8,4,5)= 1.d0/10.d0
      table_hgh(8,4,6)= -1.d0/30.d0
      table_hgh(8,4,7)= -1.d0/10.d0
      table_hgh(8,4,8)= 2.d0/15.d0

      table_hgh(8,5,1)= -1.d0/336.d0
      table_hgh(8,5,2)= -1.d0/3360.d0
      table_hgh(8,5,3)= 19.d0/1680.d0
      table_hgh(8,5,4)= -1.d0/1680.d0
      table_hgh(8,5,5)= -1.d0/560.d0
      table_hgh(8,5,6)= -1.d0/336.d0
      table_hgh(8,5,7)= 1.d0/560.d0
      table_hgh(8,5,8)= 1.d0/105.d0

      table_hgh(8,6,1)= 1.d0/42.d0
      table_hgh(8,6,2)= 1.d0/280.d0
      table_hgh(8,6,3)= 5.d0/84.d0
      table_hgh(8,6,4)= -1.d0/210.d0
      table_hgh(8,6,5)= -1.d0/20.d0
      table_hgh(8,6,6)= 0.d0
      table_hgh(8,6,7)= 1.d0/20.d0
      table_hgh(8,6,8)= 1.d0/30.d0

      table_hgh(8,7,1)= 1.d0/2520.d0
      table_hgh(8,7,2)= 1.d0/20160.d0
      table_hgh(8,7,3)= 17.d0/10080.d0
      table_hgh(8,7,4)= -1.d0/10080.d0
      table_hgh(8,7,5)= -1.d0/1120.d0
      table_hgh(8,7,6)= 0.d0
      table_hgh(8,7,7)= 1.d0/1120.d0
      table_hgh(8,7,8)= 1.d0/840.d0

      table_hgh(8,8,1)= -1.d0/336.d0
      table_hgh(8,8,2)= -1.d0/3360.d0
      table_hgh(8,8,3)= 19.d0/1680.d0
      table_hgh(8,8,4)= -1.d0/1680.d0
      table_hgh(8,8,5)= -1.d0/560.d0
      table_hgh(8,8,6)= -1.d0/336.d0
      table_hgh(8,8,7)= 1.d0/560.d0
      table_hgh(8,8,8)= 1.d0/105.d0

    !  *** triple product hhg ***

      do k=1,8
         do j=1,8
            do i=1,8
               table_hhg(i,k,j)=table_hgh(i,j,k)
            enddo
         enddo
      enddo

    !  *** triple product hgg ***

      table_hgg(1,1,1)= 1.d0/2.d0
      table_hgg(1,1,2)= 0.d0
      table_hgg(1,1,3)= -1.d0/10.d0
      table_hgg(1,1,4)= 1.d0/2.d0
      table_hgg(1,1,5)= 1.d0/48.d0
      table_hgg(1,1,6)= -1.d0/10.d0
      table_hgg(1,1,7)= -1.d0/420.d0
      table_hgg(1,1,8)= 1.d0/48.d0

      table_hgg(1,2,1)= 0.d0
      table_hgg(1,2,2)= 0.d0
      table_hgg(1,2,3)= 0.d0
      table_hgg(1,2,4)= 0.d0
      table_hgg(1,2,5)= 0.d0
      table_hgg(1,2,6)= 0.d0
      table_hgg(1,2,7)= 0.d0
      table_hgg(1,2,8)= 0.d0

      table_hgg(1,3,1)= -1.d0/10.d0
      table_hgg(1,3,2)= 0.d0
      table_hgg(1,3,3)= 1.d0/24.d0
      table_hgg(1,3,4)= -1.d0/10.d0
      table_hgg(1,3,5)= -1.d0/140.d0
      table_hgg(1,3,6)= 1.d0/24.d0
      table_hgg(1,3,7)= 1.d0/960.d0
      table_hgg(1,3,8)= -1.d0/140.d0

      table_hgg(1,4,1)= 1.d0/2.d0
      table_hgg(1,4,2)= 0.d0
      table_hgg(1,4,3)= -1.d0/10.d0
      table_hgg(1,4,4)= 1.d0/2.d0
      table_hgg(1,4,5)= 1.d0/48.d0
      table_hgg(1,4,6)= -1.d0/10.d0
      table_hgg(1,4,7)= -1.d0/420.d0
      table_hgg(1,4,8)= 1.d0/48.d0

      table_hgg(1,5,1)= 1.d0/48.d0
      table_hgg(1,5,2)= 0.d0
      table_hgg(1,5,3)= -1.d0/140.d0
      table_hgg(1,5,4)= 1.d0/48.d0
      table_hgg(1,5,5)= 1.d0/640.d0
      table_hgg(1,5,6)= -1.d0/140.d0
      table_hgg(1,5,7)= -5.d0/24192.d0
      table_hgg(1,5,8)= 1.d0/640.d0

      table_hgg(1,6,1)= -1.d0/10.d0
      table_hgg(1,6,2)= 0.d0
      table_hgg(1,6,3)= 1.d0/24.d0
      table_hgg(1,6,4)= -1.d0/10.d0
      table_hgg(1,6,5)= -1.d0/140.d0
      table_hgg(1,6,6)= 1.d0/24.d0
      table_hgg(1,6,7)= 1.d0/960.d0
      table_hgg(1,6,8)= -1.d0/140.d0

      table_hgg(1,7,1)= -1.d0/420.d0
      table_hgg(1,7,2)= 0.d0
      table_hgg(1,7,3)= 1.d0/960.d0
      table_hgg(1,7,4)= -1.d0/420.d0
      table_hgg(1,7,5)= -5.d0/24192.d0
      table_hgg(1,7,6)= 1.d0/960.d0
      table_hgg(1,7,7)= 1.d0/32256.d0
      table_hgg(1,7,8)= -5.d0/24192.d0

      table_hgg(1,8,1)= 1.d0/48.d0
      table_hgg(1,8,2)= 0.d0
      table_hgg(1,8,3)= -1.d0/140.d0
      table_hgg(1,8,4)= 1.d0/48.d0
      table_hgg(1,8,5)= 1.d0/640.d0
      table_hgg(1,8,6)= -1.d0/140.d0
      table_hgg(1,8,7)= -5.d0/24192.d0
      table_hgg(1,8,8)= 1.d0/640.d0


      table_hgg(2,1,1)= 1.d0/12.d0
      table_hgg(2,1,2)= 0.d0
      table_hgg(2,1,3)= -1.d0/120.d0
      table_hgg(2,1,4)= 1.d0/12.d0
      table_hgg(2,1,5)= 1.d0/480.d0
      table_hgg(2,1,6)= -1.d0/120.d0
      table_hgg(2,1,7)= -1.d0/6720.d0
      table_hgg(2,1,8)= 1.d0/480.d0

      table_hgg(2,2,1)= 0.d0
      table_hgg(2,2,2)= 0.d0
      table_hgg(2,2,3)= 0.d0
      table_hgg(2,2,4)= 0.d0
      table_hgg(2,2,5)= 0.d0
      table_hgg(2,2,6)= 0.d0
      table_hgg(2,2,7)= 0.d0
      table_hgg(2,2,8)= 0.d0

      table_hgg(2,3,1)= -1.d0/120.d0
      table_hgg(2,3,2)= 0.d0
      table_hgg(2,3,3)= 1.d0/240.d0
      table_hgg(2,3,4)= -1.d0/120.d0
      table_hgg(2,3,5)= -1.d0/2240.d0
      table_hgg(2,3,6)= 1.d0/240.d0
      table_hgg(2,3,7)= 1.d0/13440.d0
      table_hgg(2,3,8)= -1.d0/2240.d0

      table_hgg(2,4,1)= 1.d0/12.d0
      table_hgg(2,4,2)= 0.d0
      table_hgg(2,4,3)= -1.d0/120.d0
      table_hgg(2,4,4)= 1.d0/12.d0
      table_hgg(2,4,5)= 1.d0/480.d0
      table_hgg(2,4,6)= -1.d0/120.d0
      table_hgg(2,4,7)= -1.d0/6720.d0
      table_hgg(2,4,8)= 1.d0/480.d0

      table_hgg(2,5,1)= 1.d0/480.d0
      table_hgg(2,5,2)= 0.d0
      table_hgg(2,5,3)= -1.d0/2240.d0
      table_hgg(2,5,4)= 1.d0/480.d0
      table_hgg(2,5,5)= 1.d0/8960.d0
      table_hgg(2,5,6)= -1.d0/2240.d0
      table_hgg(2,5,7)= -1.d0/96768.d0
      table_hgg(2,5,8)= 1.d0/8960.d0

      table_hgg(2,6,1)= -1.d0/120.d0
      table_hgg(2,6,2)= 0.d0
      table_hgg(2,6,3)= 1.d0/240.d0
      table_hgg(2,6,4)= -1.d0/120.d0
      table_hgg(2,6,5)= -1.d0/2240.d0
      table_hgg(2,6,6)= 1.d0/240.d0
      table_hgg(2,6,7)= 1.d0/13440.d0
      table_hgg(2,6,8)= -1.d0/2240.d0

      table_hgg(2,7,1)= -1.d0/6720.d0
      table_hgg(2,7,2)= 0.d0
      table_hgg(2,7,3)= 1.d0/13440.d0
      table_hgg(2,7,4)= -1.d0/6720.d0
      table_hgg(2,7,5)= -1.d0/96768.d0
      table_hgg(2,7,6)= 1.d0/13440.d0
      table_hgg(2,7,7)= 1.d0/580608.d0
      table_hgg(2,7,8)= -1.d0/96768.d0

      table_hgg(2,8,1)= 1.d0/480.d0
      table_hgg(2,8,2)= 0.d0
      table_hgg(2,8,3)= -1.d0/2240.d0
      table_hgg(2,8,4)= 1.d0/480.d0
      table_hgg(2,8,5)= 1.d0/8960.d0
      table_hgg(2,8,6)= -1.d0/2240.d0
      table_hgg(2,8,7)= -1.d0/96768.d0
      table_hgg(2,8,8)= 1.d0/8960.d0


      table_hgg(3,1,1)= 1.d0/2.d0
      table_hgg(3,1,2)= 0.d0
      table_hgg(3,1,3)= 1.d0/10.d0
      table_hgg(3,1,4)= 1.d0/2.d0
      table_hgg(3,1,5)= 1.d0/48.d0
      table_hgg(3,1,6)= 1.d0/10.d0
      table_hgg(3,1,7)= 1.d0/420.d0
      table_hgg(3,1,8)= 1.d0/48.d0

      table_hgg(3,2,1)= 0.d0
      table_hgg(3,2,2)= 0.d0
      table_hgg(3,2,3)= 0.d0
      table_hgg(3,2,4)= 0.d0
      table_hgg(3,2,5)= 0.d0
      table_hgg(3,2,6)= 0.d0
      table_hgg(3,2,7)= 0.d0
      table_hgg(3,2,8)= 0.d0

      table_hgg(3,3,1)= 1.d0/10.d0
      table_hgg(3,3,2)= 0.d0
      table_hgg(3,3,3)= 1.d0/24.d0
      table_hgg(3,3,4)= 1.d0/10.d0
      table_hgg(3,3,5)= 1.d0/140.d0
      table_hgg(3,3,6)= 1.d0/24.d0
      table_hgg(3,3,7)= 1.d0/960.d0
      table_hgg(3,3,8)= 1.d0/140.d0

      table_hgg(3,4,1)= 1.d0/2.d0
      table_hgg(3,4,2)= 0.d0
      table_hgg(3,4,3)= 1.d0/10.d0
      table_hgg(3,4,4)= 1.d0/2.d0
      table_hgg(3,4,5)= 1.d0/48.d0
      table_hgg(3,4,6)= 1.d0/10.d0
      table_hgg(3,4,7)= 1.d0/420.d0
      table_hgg(3,4,8)= 1.d0/48.d0

      table_hgg(3,5,1)= 1.d0/48.d0
      table_hgg(3,5,2)= 0.d0
      table_hgg(3,5,3)= 1.d0/140.d0
      table_hgg(3,5,4)= 1.d0/48.d0
      table_hgg(3,5,5)= 1.d0/640.d0
      table_hgg(3,5,6)= 1.d0/140.d0
      table_hgg(3,5,7)= 5.d0/24192.d0
      table_hgg(3,5,8)= 1.d0/640.d0

      table_hgg(3,6,1)= 1.d0/10.d0
      table_hgg(3,6,2)= 0.d0
      table_hgg(3,6,3)= 1.d0/24.d0
      table_hgg(3,6,4)= 1.d0/10.d0
      table_hgg(3,6,5)= 1.d0/140.d0
      table_hgg(3,6,6)= 1.d0/24.d0
      table_hgg(3,6,7)= 1.d0/960.d0
      table_hgg(3,6,8)= 1.d0/140.d0

      table_hgg(3,7,1)= 1.d0/420.d0
      table_hgg(3,7,2)= 0.d0
      table_hgg(3,7,3)= 1.d0/960.d0
      table_hgg(3,7,4)= 1.d0/420.d0
      table_hgg(3,7,5)= 5.d0/24192.d0
      table_hgg(3,7,6)= 1.d0/960.d0
      table_hgg(3,7,7)= 1.d0/32256.d0
      table_hgg(3,7,8)= 5.d0/24192.d0

      table_hgg(3,8,1)= 1.d0/48.d0
      table_hgg(3,8,2)= 0.d0
      table_hgg(3,8,3)= 1.d0/140.d0
      table_hgg(3,8,4)= 1.d0/48.d0
      table_hgg(3,8,5)= 1.d0/640.d0
      table_hgg(3,8,6)= 1.d0/140.d0
      table_hgg(3,8,7)= 5.d0/24192.d0
      table_hgg(3,8,8)= 1.d0/640.d0


      table_hgg(4,1,1)= -1.d0/12.d0
      table_hgg(4,1,2)= 0.d0
      table_hgg(4,1,3)= -1.d0/120.d0
      table_hgg(4,1,4)= -1.d0/12.d0
      table_hgg(4,1,5)= -1.d0/480.d0
      table_hgg(4,1,6)= -1.d0/120.d0
      table_hgg(4,1,7)= -1.d0/6720.d0
      table_hgg(4,1,8)= -1.d0/480.d0

      table_hgg(4,2,1)= 0.d0
      table_hgg(4,2,2)= 0.d0
      table_hgg(4,2,3)= 0.d0
      table_hgg(4,2,4)= 0.d0
      table_hgg(4,2,5)= 0.d0
      table_hgg(4,2,6)= 0.d0
      table_hgg(4,2,7)= 0.d0
      table_hgg(4,2,8)= 0.d0

      table_hgg(4,3,1)= -1.d0/120.d0
      table_hgg(4,3,2)= 0.d0
      table_hgg(4,3,3)= -1.d0/240.d0
      table_hgg(4,3,4)= -1.d0/120.d0
      table_hgg(4,3,5)= -1.d0/2240.d0
      table_hgg(4,3,6)= -1.d0/240.d0
      table_hgg(4,3,7)= -1.d0/13440.d0
      table_hgg(4,3,8)= -1.d0/2240.d0

      table_hgg(4,4,1)= -1.d0/12.d0
      table_hgg(4,4,2)= 0.d0
      table_hgg(4,4,3)= -1.d0/120.d0
      table_hgg(4,4,4)= -1.d0/12.d0
      table_hgg(4,4,5)= -1.d0/480.d0
      table_hgg(4,4,6)= -1.d0/120.d0
      table_hgg(4,4,7)= -1.d0/6720.d0
      table_hgg(4,4,8)= -1.d0/480.d0

      table_hgg(4,5,1)= -1.d0/480.d0
      table_hgg(4,5,2)= 0.d0
      table_hgg(4,5,3)= -1.d0/2240.d0
      table_hgg(4,5,4)= -1.d0/480.d0
      table_hgg(4,5,5)= -1.d0/8960.d0
      table_hgg(4,5,6)= -1.d0/2240.d0
      table_hgg(4,5,7)= -1.d0/96768.d0
      table_hgg(4,5,8)= -1.d0/8960.d0

      table_hgg(4,6,1)= -1.d0/120.d0
      table_hgg(4,6,2)= 0.d0
      table_hgg(4,6,3)= -1.d0/240.d0
      table_hgg(4,6,4)= -1.d0/120.d0
      table_hgg(4,6,5)= -1.d0/2240.d0
      table_hgg(4,6,6)= -1.d0/240.d0
      table_hgg(4,6,7)= -1.d0/13440.d0
      table_hgg(4,6,8)= -1.d0/2240.d0

      table_hgg(4,7,1)= -1.d0/6720.d0
      table_hgg(4,7,2)= 0.d0
      table_hgg(4,7,3)= -1.d0/13440.d0
      table_hgg(4,7,4)= -1.d0/6720.d0
      table_hgg(4,7,5)= -1.d0/96768.d0
      table_hgg(4,7,6)= -1.d0/13440.d0
      table_hgg(4,7,7)= -1.d0/580608.d0
      table_hgg(4,7,8)= -1.d0/96768.d0

      table_hgg(4,8,1)= -1.d0/480.d0
      table_hgg(4,8,2)= 0.d0
      table_hgg(4,8,3)= -1.d0/2240.d0
      table_hgg(4,8,4)= -1.d0/480.d0
      table_hgg(4,8,5)= -1.d0/8960.d0
      table_hgg(4,8,6)= -1.d0/2240.d0
      table_hgg(4,8,7)= -1.d0/96768.d0
      table_hgg(4,8,8)= -1.d0/8960.d0


      table_hgg(5,1,1)= -1.d0
      table_hgg(5,1,2)= 0.d0
      table_hgg(5,1,3)= 0.d0
      table_hgg(5,1,4)= -1.d0
      table_hgg(5,1,5)= -1.d0/40.d0
      table_hgg(5,1,6)= 0.d0
      table_hgg(5,1,7)= 0.d0
      table_hgg(5,1,8)= -1.d0/40.d0

      table_hgg(5,2,1)= 0.d0
      table_hgg(5,2,2)= 0.d0
      table_hgg(5,2,3)= 0.d0
      table_hgg(5,2,4)= 0.d0
      table_hgg(5,2,5)= 0.d0
      table_hgg(5,2,6)= 0.d0
      table_hgg(5,2,7)= 0.d0
      table_hgg(5,2,8)= 0.d0

      table_hgg(5,3,1)= 0.d0
      table_hgg(5,3,2)= 0.d0
      table_hgg(5,3,3)= -1.d0/20.d0
      table_hgg(5,3,4)= 0.d0
      table_hgg(5,3,5)= 0.d0
      table_hgg(5,3,6)= -1.d0/20.d0
      table_hgg(5,3,7)= -1.d0/1120.d0
      table_hgg(5,3,8)= 0.d0

      table_hgg(5,4,1)= -1.d0
      table_hgg(5,4,2)= 0.d0
      table_hgg(5,4,3)= 0.d0
      table_hgg(5,4,4)= -1.d0
      table_hgg(5,4,5)= -1.d0/40.d0
      table_hgg(5,4,6)= 0.d0
      table_hgg(5,4,7)= 0.d0
      table_hgg(5,4,8)= -1.d0/40.d0

      table_hgg(5,5,1)= -1.d0/40.d0
      table_hgg(5,5,2)= 0.d0
      table_hgg(5,5,3)= 0.d0
      table_hgg(5,5,4)= -1.d0/40.d0
      table_hgg(5,5,5)= -3.d0/2240.d0
      table_hgg(5,5,6)= 0.d0
      table_hgg(5,5,7)= 0.d0
      table_hgg(5,5,8)= -3.d0/2240.d0

      table_hgg(5,6,1)= 0.d0
      table_hgg(5,6,2)= 0.d0
      table_hgg(5,6,3)= -1.d0/20.d0
      table_hgg(5,6,4)= 0.d0
      table_hgg(5,6,5)= 0.d0
      table_hgg(5,6,6)= -1.d0/20.d0
      table_hgg(5,6,7)= -1.d0/1120.d0
      table_hgg(5,6,8)= 0.d0

      table_hgg(5,7,1)= 0.d0
      table_hgg(5,7,2)= 0.d0
      table_hgg(5,7,3)= -1.d0/1120.d0
      table_hgg(5,7,4)= 0.d0
      table_hgg(5,7,5)= 0.d0
      table_hgg(5,7,6)= -1.d0/1120.d0
      table_hgg(5,7,7)= -1.d0/48384.d0
      table_hgg(5,7,8)= 0.d0

      table_hgg(5,8,1)= -1.d0/40.d0
      table_hgg(5,8,2)= 0.d0
      table_hgg(5,8,3)= 0.d0
      table_hgg(5,8,4)= -1.d0/40.d0
      table_hgg(5,8,5)= -3.d0/2240.d0
      table_hgg(5,8,6)= 0.d0
      table_hgg(5,8,7)= 0.d0
      table_hgg(5,8,8)= -3.d0/2240.d0


      table_hgg(6,1,1)= 0.d0
      table_hgg(6,1,2)= 0.d0
      table_hgg(6,1,3)= -1.d0/12.d0
      table_hgg(6,1,4)= 0.d0
      table_hgg(6,1,5)= 1.d0/120.d0
      table_hgg(6,1,6)= -1.d0/12.d0
      table_hgg(6,1,7)= -1.d0/480.d0
      table_hgg(6,1,8)= 1.d0/120.d0

      table_hgg(6,2,1)= 0.d0
      table_hgg(6,2,2)= 0.d0
      table_hgg(6,2,3)= 0.d0
      table_hgg(6,2,4)= 0.d0
      table_hgg(6,2,5)= 0.d0
      table_hgg(6,2,6)= 0.d0
      table_hgg(6,2,7)= 0.d0
      table_hgg(6,2,8)= 0.d0

      table_hgg(6,3,1)= -1.d0/12.d0
      table_hgg(6,3,2)= 0.d0
      table_hgg(6,3,3)= 1.d0/60.d0
      table_hgg(6,3,4)= -1.d0/12.d0
      table_hgg(6,3,5)= -1.d0/160.d0
      table_hgg(6,3,6)= 1.d0/60.d0
      table_hgg(6,3,7)= 1.d0/1680.d0
      table_hgg(6,3,8)= -1.d0/160.d0

      table_hgg(6,4,1)= 0.d0
      table_hgg(6,4,2)= 0.d0
      table_hgg(6,4,3)= -1.d0/12.d0
      table_hgg(6,4,4)= 0.d0
      table_hgg(6,4,5)= 1.d0/120.d0
      table_hgg(6,4,6)= -1.d0/12.d0
      table_hgg(6,4,7)= -1.d0/480.d0
      table_hgg(6,4,8)= 1.d0/120.d0

      table_hgg(6,5,1)= 1.d0/120.d0
      table_hgg(6,5,2)= 0.d0
      table_hgg(6,5,3)= -1.d0/160.d0
      table_hgg(6,5,4)= 1.d0/120.d0
      table_hgg(6,5,5)= 1.d0/1120.d0
      table_hgg(6,5,6)= -1.d0/160.d0
      table_hgg(6,5,7)= -1.d0/5376.d0
      table_hgg(6,5,8)= 1.d0/1120.d0

      table_hgg(6,6,1)= -1.d0/12.d0
      table_hgg(6,6,2)= 0.d0
      table_hgg(6,6,3)= 1.d0/60.d0
      table_hgg(6,6,4)= -1.d0/12.d0
      table_hgg(6,6,5)= -1.d0/160.d0
      table_hgg(6,6,6)= 1.d0/60.d0
      table_hgg(6,6,7)= 1.d0/1680.d0
      table_hgg(6,6,8)= -1.d0/160.d0

      table_hgg(6,7,1)= -1.d0/480.d0
      table_hgg(6,7,2)= 0.d0
      table_hgg(6,7,3)= 1.d0/1680.d0
      table_hgg(6,7,4)= -1.d0/480.d0
      table_hgg(6,7,5)= -1.d0/5376.d0
      table_hgg(6,7,6)= 1.d0/1680.d0
      table_hgg(6,7,7)= 1.d0/48384.d0
      table_hgg(6,7,8)= -1.d0/5376.d0

      table_hgg(6,8,1)= 1.d0/120.d0
      table_hgg(6,8,2)= 0.d0
      table_hgg(6,8,3)= -1.d0/160.d0
      table_hgg(6,8,4)= 1.d0/120.d0
      table_hgg(6,8,5)= 1.d0/1120.d0
      table_hgg(6,8,6)= -1.d0/160.d0
      table_hgg(6,8,7)= -1.d0/5376.d0
      table_hgg(6,8,8)= 1.d0/1120.d0


      table_hgg(7,1,1)= 1.d0
      table_hgg(7,1,2)= 0.d0
      table_hgg(7,1,3)= 0.d0
      table_hgg(7,1,4)= 1.d0
      table_hgg(7,1,5)= 1.d0/40.d0
      table_hgg(7,1,6)= 0.d0
      table_hgg(7,1,7)= 0.d0
      table_hgg(7,1,8)= 1.d0/40.d0

      table_hgg(7,2,1)= 0.d0
      table_hgg(7,2,2)= 0.d0
      table_hgg(7,2,3)= 0.d0
      table_hgg(7,2,4)= 0.d0
      table_hgg(7,2,5)= 0.d0
      table_hgg(7,2,6)= 0.d0
      table_hgg(7,2,7)= 0.d0
      table_hgg(7,2,8)= 0.d0

      table_hgg(7,3,1)= 0.d0
      table_hgg(7,3,2)= 0.d0
      table_hgg(7,3,3)= 1.d0/20.d0
      table_hgg(7,3,4)= 0.d0
      table_hgg(7,3,5)= 0.d0
      table_hgg(7,3,6)= 1.d0/20.d0
      table_hgg(7,3,7)= 1.d0/1120.d0
      table_hgg(7,3,8)= 0.d0

      table_hgg(7,4,1)= 1.d0
      table_hgg(7,4,2)= 0.d0
      table_hgg(7,4,3)= 0.d0
      table_hgg(7,4,4)= 1.d0
      table_hgg(7,4,5)= 1.d0/40.d0
      table_hgg(7,4,6)= 0.d0
      table_hgg(7,4,7)= 0.d0
      table_hgg(7,4,8)= 1.d0/40.d0

      table_hgg(7,5,1)= 1.d0/40.d0
      table_hgg(7,5,2)= 0.d0
      table_hgg(7,5,3)= 0.d0
      table_hgg(7,5,4)= 1.d0/40.d0
      table_hgg(7,5,5)= 3.d0/2240.d0
      table_hgg(7,5,6)= 0.d0
      table_hgg(7,5,7)= 0.d0
      table_hgg(7,5,8)= 3.d0/2240.d0

      table_hgg(7,6,1)= 0.d0
      table_hgg(7,6,2)= 0.d0
      table_hgg(7,6,3)= 1.d0/20.d0
      table_hgg(7,6,4)= 0.d0
      table_hgg(7,6,5)= 0.d0
      table_hgg(7,6,6)= 1.d0/20.d0
      table_hgg(7,6,7)= 1.d0/1120.d0
      table_hgg(7,6,8)= 0.d0

      table_hgg(7,7,1)= 0.d0
      table_hgg(7,7,2)= 0.d0
      table_hgg(7,7,3)= 1.d0/1120.d0
      table_hgg(7,7,4)= 0.d0
      table_hgg(7,7,5)= 0.d0
      table_hgg(7,7,6)= 1.d0/1120.d0
      table_hgg(7,7,7)= 1.d0/48384.d0
      table_hgg(7,7,8)= 0.d0

      table_hgg(7,8,1)= 1.d0/40.d0
      table_hgg(7,8,2)= 0.d0
      table_hgg(7,8,3)= 0.d0
      table_hgg(7,8,4)= 1.d0/40.d0
      table_hgg(7,8,5)= 3.d0/2240.d0
      table_hgg(7,8,6)= 0.d0
      table_hgg(7,8,7)= 0.d0
      table_hgg(7,8,8)= 3.d0/2240.d0


      table_hgg(8,1,1)= 0.d0
      table_hgg(8,1,2)= 0.d0
      table_hgg(8,1,3)= 1.d0/12.d0
      table_hgg(8,1,4)= 0.d0
      table_hgg(8,1,5)= 1.d0/120.d0
      table_hgg(8,1,6)= 1.d0/12.d0
      table_hgg(8,1,7)= 1.d0/480.d0
      table_hgg(8,1,8)= 1.d0/120.d0

      table_hgg(8,2,1)= 0.d0
      table_hgg(8,2,2)= 0.d0
      table_hgg(8,2,3)= 0.d0
      table_hgg(8,2,4)= 0.d0
      table_hgg(8,2,5)= 0.d0
      table_hgg(8,2,6)= 0.d0
      table_hgg(8,2,7)= 0.d0
      table_hgg(8,2,8)= 0.d0

      table_hgg(8,3,1)= 1.d0/12.d0
      table_hgg(8,3,2)= 0.d0
      table_hgg(8,3,3)= 1.d0/60.d0
      table_hgg(8,3,4)= 1.d0/12.d0
      table_hgg(8,3,5)= 1.d0/160.d0
      table_hgg(8,3,6)= 1.d0/60.d0
      table_hgg(8,3,7)= 1.d0/1680.d0
      table_hgg(8,3,8)= 1.d0/160.d0

      table_hgg(8,4,1)= 0.d0
      table_hgg(8,4,2)= 0.d0
      table_hgg(8,4,3)= 1.d0/12.d0
      table_hgg(8,4,4)= 0.d0
      table_hgg(8,4,5)= 1.d0/120.d0
      table_hgg(8,4,6)= 1.d0/12.d0
      table_hgg(8,4,7)= 1.d0/480.d0
      table_hgg(8,4,8)= 1.d0/120.d0

      table_hgg(8,5,1)= 1.d0/120.d0
      table_hgg(8,5,2)= 0.d0
      table_hgg(8,5,3)= 1.d0/160.d0
      table_hgg(8,5,4)= 1.d0/120.d0
      table_hgg(8,5,5)= 1.d0/1120.d0
      table_hgg(8,5,6)= 1.d0/160.d0
      table_hgg(8,5,7)= 1.d0/5376.d0
      table_hgg(8,5,8)= 1.d0/1120.d0

      table_hgg(8,6,1)= 1.d0/12.d0
      table_hgg(8,6,2)= 0.d0
      table_hgg(8,6,3)= 1.d0/60.d0
      table_hgg(8,6,4)= 1.d0/12.d0
      table_hgg(8,6,5)= 1.d0/160.d0
      table_hgg(8,6,6)= 1.d0/60.d0
      table_hgg(8,6,7)= 1.d0/1680.d0
      table_hgg(8,6,8)= 1.d0/160.d0

      table_hgg(8,7,1)= 1.d0/480.d0
      table_hgg(8,7,2)= 0.d0
      table_hgg(8,7,3)= 1.d0/1680.d0
      table_hgg(8,7,4)= 1.d0/480.d0
      table_hgg(8,7,5)= 1.d0/5376.d0
      table_hgg(8,7,6)= 1.d0/1680.d0
      table_hgg(8,7,7)= 1.d0/48384.d0
      table_hgg(8,7,8)= 1.d0/5376.d0

      table_hgg(8,8,1)= 1.d0/120.d0
      table_hgg(8,8,2)= 0.d0
      table_hgg(8,8,3)= 1.d0/160.d0
      table_hgg(8,8,4)= 1.d0/120.d0
      table_hgg(8,8,5)= 1.d0/1120.d0
      table_hgg(8,8,6)= 1.d0/160.d0
      table_hgg(8,8,7)= 1.d0/5376.d0
      table_hgg(8,8,8)= 1.d0/1120.d0

!----------
      table_hq(1,1)= 11.d0/60.d0
      table_hq(1,2)= 1.d0/3.d0
      table_hq(1,3)= -1.d0/60.d0
      table_hq(1,4)= -9.d0/10.d0
      table_hq(1,5)= 4.d0/5.d0
      table_hq(1,6)= 1.d0/10.d0

      table_hq(2,1)= 1.d0/60.d0
      table_hq(2,2)= 1.d0/15.d0
      table_hq(2,3)= 0.d0
      table_hq(2,4)= -7.d0/60.d0
      table_hq(2,5)= 1.d0/15.d0
      table_hq(2,6)= 1.d0/20.d0

      table_hq(3,1)= -1.d0/60.d0
      table_hq(3,2)= 1.d0/3.d0
      table_hq(3,3)= 11.d0/60.d0
      table_hq(3,4)= -1.d0/10.d0
      table_hq(3,5)= -4.d0/5.d0
      table_hq(3,6)= 9.d0/10.d0

      table_hq(4,1)= 0.d0
      table_hq(4,2)= -1.d0/15.d0
      table_hq(4,3)= -1.d0/60.d0
      table_hq(4,4)= 1.d0/20.d0
      table_hq(4,5)= 1.d0/15.d0
      table_hq(4,6)= -7.d0/60.d0

      table_hq(5,1)= -1.d0/10.d0
      table_hq(5,2)= -4.d0/5.d0
      table_hq(5,3)= -1.d0/10.d0
      table_hq(5,4)= 1.d0
      table_hq(5,5)= 0.d0
      table_hq(5,6)= -1.d0

      table_hq(6,1)= 7.d0/60.d0
      table_hq(6,2)= -1.d0/15.d0
      table_hq(6,3)= -1.d0/20.d0
      table_hq(6,4)= -1.d0/3.d0
      table_hq(6,5)= 2.d0/3.d0
      table_hq(6,6)= -1.d0/3.d0

      table_hq(7,1)= 1.d0/10.d0
      table_hq(7,2)= 4.d0/5.d0
      table_hq(7,3)= 1.d0/10.d0
      table_hq(7,4)= -1.d0
      table_hq(7,5)= 0.d0
      table_hq(7,6)= 1.d0

      table_hq(8,1)= -1.d0/20.d0
      table_hq(8,2)= -1.d0/15.d0
      table_hq(8,3)= 7.d0/60.d0
      table_hq(8,4)= 1.d0/3.d0
      table_hq(8,5)= -2.d0/3.d0
      table_hq(8,6)= 1.d0/3.d0

!-------
      table_qq(1,1)= 2.d0/15.d0
      table_qq(1,2)= 1.d0/15.d0
      table_qq(1,3)= -1.d0/30.d0
      table_qq(1,4)= -1.d0/2.d0
      table_qq(1,5)= 2.d0/3.d0
      table_qq(1,6)= -1.d0/6.d0

      table_qq(2,1)= 1.d0/15.d0
      table_qq(2,2)= 8.d0/15.d0
      table_qq(2,3)= 1.d0/15.d0
      table_qq(2,4)= -2.d0/3.d0
      table_qq(2,5)= 0.d0
      table_qq(2,6)= 2.d0/3.d0

      table_qq(3,1)= -1.d0/30.d0
      table_qq(3,2)= 1.d0/15.d0
      table_qq(3,3)= 2.d0/15.d0
      table_qq(3,4)= 1.d0/6.d0
      table_qq(3,5)= -2.d0/3.d0
      table_qq(3,6)= 1.d0/2.d0

      table_qq(4,1)= -1.d0/2.d0
      table_qq(4,2)= -2.d0/3.d0
      table_qq(4,3)= 1.d0/6.d0
      table_qq(4,4)= 7.d0/3.d0
      table_qq(4,5)= -8.d0/3.d0
      table_qq(4,6)= 1.d0/3.d0

      table_qq(5,1)= 2.d0/3.d0
      table_qq(5,2)= 0.d0
      table_qq(5,3)= -2.d0/3.d0
      table_qq(5,4)= -8.d0/3.d0
      table_qq(5,5)= 16.d0/3.d0
      table_qq(5,6)= -8.d0/3.d0

      table_qq(6,1)= -1.d0/6.d0
      table_qq(6,2)= 2.d0/3.d0
      table_qq(6,3)= 1.d0/2.d0
      table_qq(6,4)= 1.d0/3.d0
      table_qq(6,5)= -8.d0/3.d0
      table_qq(6,6)= 7.d0/3.d0

!------

      table_hqh(1,1,1)= 11.d0/63.d0
      table_hqh(1,1,2)= 1.d0/63.d0
      table_hqh(1,1,3)= 11.d0/1260.d0
      table_hqh(1,1,4)= -1.d0/315.d0
      table_hqh(1,1,5)= -4.d0/35.d0
      table_hqh(1,1,6)= 11.d0/105.d0
      table_hqh(1,1,7)= 4.d0/35.d0
      table_hqh(1,1,8)= -1.d0/28.d0

      table_hqh(1,2,1)= 2.d0/9.d0
      table_hqh(1,2,2)= 5.d0/126.d0
      table_hqh(1,2,3)= 1.d0/9.d0
      table_hqh(1,2,4)= -17.d0/630.d0
      table_hqh(1,2,5)= -2.d0/5.d0
      table_hqh(1,2,6)= 1.d0/105.d0
      table_hqh(1,2,7)= 2.d0/5.d0
      table_hqh(1,2,8)= -8.d0/105.d0

      table_hqh(1,3,1)= -8.d0/315.d0
      table_hqh(1,3,2)= -1.d0/315.d0
      table_hqh(1,3,3)= 11.d0/1260.d0
      table_hqh(1,3,4)= -1.d0/1260.d0
      table_hqh(1,3,5)= 1.d0/70.d0
      table_hqh(1,3,6)= -1.d0/70.d0
      table_hqh(1,3,7)= -1.d0/70.d0
      table_hqh(1,3,8)= 1.d0/84.d0

      table_hqh(1,4,1)= -27.d0/35.d0
      table_hqh(1,4,2)= -19.d0/210.d0
      table_hqh(1,4,3)= -9.d0/70.d0
      table_hqh(1,4,4)= 1.d0/28.d0
      table_hqh(1,4,5)= 53.d0/70.d0
      table_hqh(1,4,6)= -71.d0/210.d0
      table_hqh(1,4,7)= -53.d0/70.d0
      table_hqh(1,4,8)= 41.d0/210.d0

      table_hqh(1,5,1)= 4.d0/5.d0
      table_hqh(1,5,2)= 8.d0/105.d0
      table_hqh(1,5,3)= 0.d0
      table_hqh(1,5,4)= -1.d0/105.d0
      table_hqh(1,5,5)= -18.d0/35.d0
      table_hqh(1,5,6)= 10.d0/21.d0
      table_hqh(1,5,7)= 18.d0/35.d0
      table_hqh(1,5,8)= -4.d0/21.d0

      table_hqh(1,6,1)= -1.d0/35.d0
      table_hqh(1,6,2)= 1.d0/70.d0
      table_hqh(1,6,3)= 9.d0/70.d0
      table_hqh(1,6,4)= -11.d0/420.d0
      table_hqh(1,6,5)= -17.d0/70.d0
      table_hqh(1,6,6)= -29.d0/210.d0
      table_hqh(1,6,7)= 17.d0/70.d0
      table_hqh(1,6,8)= -1.d0/210.d0


      table_hqh(2,1,1)= 1.d0/63.d0
      table_hqh(2,1,2)= 1.d0/504.d0
      table_hqh(2,1,3)= 1.d0/1260.d0
      table_hqh(2,1,4)= -1.d0/2520.d0
      table_hqh(2,1,5)= -1.d0/70.d0
      table_hqh(2,1,6)= 1.d0/140.d0
      table_hqh(2,1,7)= 1.d0/70.d0
      table_hqh(2,1,8)= -1.d0/210.d0

      table_hqh(2,2,1)= 5.d0/126.d0
      table_hqh(2,2,2)= 1.d0/126.d0
      table_hqh(2,2,3)= 17.d0/630.d0
      table_hqh(2,2,4)= -2.d0/315.d0
      table_hqh(2,2,5)= -3.d0/35.d0
      table_hqh(2,2,6)= -1.d0/210.d0
      table_hqh(2,2,7)= 3.d0/35.d0
      table_hqh(2,2,8)= -1.d0/70.d0

      table_hqh(2,3,1)= -1.d0/315.d0
      table_hqh(2,3,2)= -1.d0/2520.d0
      table_hqh(2,3,3)= 1.d0/315.d0
      table_hqh(2,3,4)= -1.d0/2520.d0
      table_hqh(2,3,5)= 0.d0
      table_hqh(2,3,6)= -1.d0/420.d0
      table_hqh(2,3,7)= 0.d0
      table_hqh(2,3,8)= 1.d0/420.d0

      table_hqh(2,4,1)= -19.d0/210.d0
      table_hqh(2,4,2)= -1.d0/70.d0
      table_hqh(2,4,3)= -11.d0/420.d0
      table_hqh(2,4,4)= 1.d0/140.d0
      table_hqh(2,4,5)= 9.d0/70.d0
      table_hqh(2,4,6)= -2.d0/105.d0
      table_hqh(2,4,7)= -9.d0/70.d0
      table_hqh(2,4,8)= 13.d0/420.d0

      table_hqh(2,5,1)= 8.d0/105.d0
      table_hqh(2,5,2)= 1.d0/105.d0
      table_hqh(2,5,3)= -1.d0/105.d0
      table_hqh(2,5,4)= 0.d0
      table_hqh(2,5,5)= -2.d0/35.d0
      table_hqh(2,5,6)= 4.d0/105.d0
      table_hqh(2,5,7)= 2.d0/35.d0
      table_hqh(2,5,8)= -1.d0/35.d0

      table_hqh(2,6,1)= 1.d0/70.d0
      table_hqh(2,6,2)= 1.d0/210.d0
      table_hqh(2,6,3)= 1.d0/28.d0
      table_hqh(2,6,4)= -1.d0/140.d0
      table_hqh(2,6,5)= -1.d0/14.d0
      table_hqh(2,6,6)= -2.d0/105.d0
      table_hqh(2,6,7)= 1.d0/14.d0
      table_hqh(2,6,8)= -1.d0/420.d0


      table_hqh(3,1,1)= 11.d0/1260.d0
      table_hqh(3,1,2)= 1.d0/1260.d0
      table_hqh(3,1,3)= -8.d0/315.d0
      table_hqh(3,1,4)= 1.d0/315.d0
      table_hqh(3,1,5)= 1.d0/70.d0
      table_hqh(3,1,6)= 1.d0/84.d0
      table_hqh(3,1,7)= -1.d0/70.d0
      table_hqh(3,1,8)= -1.d0/70.d0

      table_hqh(3,2,1)= 1.d0/9.d0
      table_hqh(3,2,2)= 17.d0/630.d0
      table_hqh(3,2,3)= 2.d0/9.d0
      table_hqh(3,2,4)= -5.d0/126.d0
      table_hqh(3,2,5)= -2.d0/5.d0
      table_hqh(3,2,6)= -8.d0/105.d0
      table_hqh(3,2,7)= 2.d0/5.d0
      table_hqh(3,2,8)= 1.d0/105.d0

      table_hqh(3,3,1)= 11.d0/1260.d0
      table_hqh(3,3,2)= 1.d0/315.d0
      table_hqh(3,3,3)= 11.d0/63.d0
      table_hqh(3,3,4)= -1.d0/63.d0
      table_hqh(3,3,5)= -4.d0/35.d0
      table_hqh(3,3,6)= -1.d0/28.d0
      table_hqh(3,3,7)= 4.d0/35.d0
      table_hqh(3,3,8)= 11.d0/105.d0

      table_hqh(3,4,1)= -9.d0/70.d0
      table_hqh(3,4,2)= -11.d0/420.d0
      table_hqh(3,4,3)= 1.d0/35.d0
      table_hqh(3,4,4)= 1.d0/70.d0
      table_hqh(3,4,5)= 17.d0/70.d0
      table_hqh(3,4,6)= 1.d0/210.d0
      table_hqh(3,4,7)= -17.d0/70.d0
      table_hqh(3,4,8)= 29.d0/210.d0

      table_hqh(3,5,1)= 0.d0
      table_hqh(3,5,2)= -1.d0/105.d0
      table_hqh(3,5,3)= -4.d0/5.d0
      table_hqh(3,5,4)= 8.d0/105.d0
      table_hqh(3,5,5)= 18.d0/35.d0
      table_hqh(3,5,6)= 4.d0/21.d0
      table_hqh(3,5,7)= -18.d0/35.d0
      table_hqh(3,5,8)= -10.d0/21.d0

      table_hqh(3,6,1)= 9.d0/70.d0
      table_hqh(3,6,2)= 1.d0/28.d0
      table_hqh(3,6,3)= 27.d0/35.d0
      table_hqh(3,6,4)= -19.d0/210.d0
      table_hqh(3,6,5)= -53.d0/70.d0
      table_hqh(3,6,6)= -41.d0/210.d0
      table_hqh(3,6,7)= 53.d0/70.d0
      table_hqh(3,6,8)= 71.d0/210.d0


      table_hqh(4,1,1)= -1.d0/315.d0
      table_hqh(4,1,2)= -1.d0/2520.d0
      table_hqh(4,1,3)= 1.d0/315.d0
      table_hqh(4,1,4)= -1.d0/2520.d0
      table_hqh(4,1,5)= 0.d0
      table_hqh(4,1,6)= -1.d0/420.d0
      table_hqh(4,1,7)= 0.d0
      table_hqh(4,1,8)= 1.d0/420.d0

      table_hqh(4,2,1)= -17.d0/630.d0
      table_hqh(4,2,2)= -2.d0/315.d0
      table_hqh(4,2,3)= -5.d0/126.d0
      table_hqh(4,2,4)= 1.d0/126.d0
      table_hqh(4,2,5)= 3.d0/35.d0
      table_hqh(4,2,6)= 1.d0/70.d0
      table_hqh(4,2,7)= -3.d0/35.d0
      table_hqh(4,2,8)= 1.d0/210.d0

      table_hqh(4,3,1)= -1.d0/1260.d0
      table_hqh(4,3,2)= -1.d0/2520.d0
      table_hqh(4,3,3)= -1.d0/63.d0
      table_hqh(4,3,4)= 1.d0/504.d0
      table_hqh(4,3,5)= 1.d0/70.d0
      table_hqh(4,3,6)= 1.d0/210.d0
      table_hqh(4,3,7)= -1.d0/70.d0
      table_hqh(4,3,8)= -1.d0/140.d0

      table_hqh(4,4,1)= 1.d0/28.d0
      table_hqh(4,4,2)= 1.d0/140.d0
      table_hqh(4,4,3)= 1.d0/70.d0
      table_hqh(4,4,4)= -1.d0/210.d0
      table_hqh(4,4,5)= -1.d0/14.d0
      table_hqh(4,4,6)= -1.d0/420.d0
      table_hqh(4,4,7)= 1.d0/14.d0
      table_hqh(4,4,8)= -2.d0/105.d0

      table_hqh(4,5,1)= -1.d0/105.d0
      table_hqh(4,5,2)= 0.d0
      table_hqh(4,5,3)= 8.d0/105.d0
      table_hqh(4,5,4)= -1.d0/105.d0
      table_hqh(4,5,5)= -2.d0/35.d0
      table_hqh(4,5,6)= -1.d0/35.d0
      table_hqh(4,5,7)= 2.d0/35.d0
      table_hqh(4,5,8)= 4.d0/105.d0

      table_hqh(4,6,1)= -11.d0/420.d0
      table_hqh(4,6,2)= -1.d0/140.d0
      table_hqh(4,6,3)= -19.d0/210.d0
      table_hqh(4,6,4)= 1.d0/70.d0
      table_hqh(4,6,5)= 9.d0/70.d0
      table_hqh(4,6,6)= 13.d0/420.d0
      table_hqh(4,6,7)= -9.d0/70.d0
      table_hqh(4,6,8)= -2.d0/105.d0


      table_hqh(5,1,1)= -4.d0/35.d0
      table_hqh(5,1,2)= -1.d0/70.d0
      table_hqh(5,1,3)= 1.d0/70.d0
      table_hqh(5,1,4)= 0.d0
      table_hqh(5,1,5)= 3.d0/35.d0
      table_hqh(5,1,6)= -2.d0/35.d0
      table_hqh(5,1,7)= -3.d0/35.d0
      table_hqh(5,1,8)= 3.d0/70.d0

      table_hqh(5,2,1)= -2.d0/5.d0
      table_hqh(5,2,2)= -3.d0/35.d0
      table_hqh(5,2,3)= -2.d0/5.d0
      table_hqh(5,2,4)= 3.d0/35.d0
      table_hqh(5,2,5)= 36.d0/35.d0
      table_hqh(5,2,6)= 4.d0/35.d0
      table_hqh(5,2,7)= -36.d0/35.d0
      table_hqh(5,2,8)= 4.d0/35.d0

      table_hqh(5,3,1)= 1.d0/70.d0
      table_hqh(5,3,2)= 0.d0
      table_hqh(5,3,3)= -4.d0/35.d0
      table_hqh(5,3,4)= 1.d0/70.d0
      table_hqh(5,3,5)= 3.d0/35.d0
      table_hqh(5,3,6)= 3.d0/70.d0
      table_hqh(5,3,7)= -3.d0/35.d0
      table_hqh(5,3,8)= -2.d0/35.d0

      table_hqh(5,4,1)= 53.d0/70.d0
      table_hqh(5,4,2)= 9.d0/70.d0
      table_hqh(5,4,3)= 17.d0/70.d0
      table_hqh(5,4,4)= -1.d0/14.d0
      table_hqh(5,4,5)= -6.d0/5.d0
      table_hqh(5,4,6)= 1.d0/10.d0
      table_hqh(5,4,7)= 6.d0/5.d0
      table_hqh(5,4,8)= -3.d0/10.d0

      table_hqh(5,5,1)= -18.d0/35.d0
      table_hqh(5,5,2)= -2.d0/35.d0
      table_hqh(5,5,3)= 18.d0/35.d0
      table_hqh(5,5,4)= -2.d0/35.d0
      table_hqh(5,5,5)= 0.d0
      table_hqh(5,5,6)= -2.d0/5.d0
      table_hqh(5,5,7)= 0.d0
      table_hqh(5,5,8)= 2.d0/5.d0

      table_hqh(5,6,1)= -17.d0/70.d0
      table_hqh(5,6,2)= -1.d0/14.d0
      table_hqh(5,6,3)= -53.d0/70.d0
      table_hqh(5,6,4)= 9.d0/70.d0
      table_hqh(5,6,5)= 6.d0/5.d0
      table_hqh(5,6,6)= 3.d0/10.d0
      table_hqh(5,6,7)= -6.d0/5.d0
      table_hqh(5,6,8)= -1.d0/10.d0


      table_hqh(6,1,1)= 11.d0/105.d0
      table_hqh(6,1,2)= 1.d0/140.d0
      table_hqh(6,1,3)= 1.d0/84.d0
      table_hqh(6,1,4)= -1.d0/420.d0
      table_hqh(6,1,5)= -2.d0/35.d0
      table_hqh(6,1,6)= 1.d0/14.d0
      table_hqh(6,1,7)= 2.d0/35.d0
      table_hqh(6,1,8)= -1.d0/84.d0

      table_hqh(6,2,1)= 1.d0/105.d0
      table_hqh(6,2,2)= -1.d0/210.d0
      table_hqh(6,2,3)= -8.d0/105.d0
      table_hqh(6,2,4)= 1.d0/70.d0
      table_hqh(6,2,5)= 4.d0/35.d0
      table_hqh(6,2,6)= 2.d0/35.d0
      table_hqh(6,2,7)= -4.d0/35.d0
      table_hqh(6,2,8)= -1.d0/105.d0

      table_hqh(6,3,1)= -1.d0/70.d0
      table_hqh(6,3,2)= -1.d0/420.d0
      table_hqh(6,3,3)= -1.d0/28.d0
      table_hqh(6,3,4)= 1.d0/210.d0
      table_hqh(6,3,5)= 3.d0/70.d0
      table_hqh(6,3,6)= 1.d0/210.d0
      table_hqh(6,3,7)= -3.d0/70.d0
      table_hqh(6,3,8)= -1.d0/84.d0

      table_hqh(6,4,1)= -71.d0/210.d0
      table_hqh(6,4,2)= -2.d0/105.d0
      table_hqh(6,4,3)= 1.d0/210.d0
      table_hqh(6,4,4)= -1.d0/420.d0
      table_hqh(6,4,5)= 1.d0/10.d0
      table_hqh(6,4,6)= -4.d0/15.d0
      table_hqh(6,4,7)= -1.d0/10.d0
      table_hqh(6,4,8)= 1.d0/30.d0

      table_hqh(6,5,1)= 10.d0/21.d0
      table_hqh(6,5,2)= 4.d0/105.d0
      table_hqh(6,5,3)= 4.d0/21.d0
      table_hqh(6,5,4)= -1.d0/35.d0
      table_hqh(6,5,5)= -2.d0/5.d0
      table_hqh(6,5,6)= 4.d0/15.d0
      table_hqh(6,5,7)= 2.d0/5.d0
      table_hqh(6,5,8)= 0.d0

      table_hqh(6,6,1)= -29.d0/210.d0
      table_hqh(6,6,2)= -2.d0/105.d0
      table_hqh(6,6,3)= -41.d0/210.d0
      table_hqh(6,6,4)= 13.d0/420.d0
      table_hqh(6,6,5)= 3.d0/10.d0
      table_hqh(6,6,6)= 0.d0
      table_hqh(6,6,7)= -3.d0/10.d0
      table_hqh(6,6,8)= -1.d0/30.d0


      table_hqh(7,1,1)= 4.d0/35.d0
      table_hqh(7,1,2)= 1.d0/70.d0
      table_hqh(7,1,3)= -1.d0/70.d0
      table_hqh(7,1,4)= 0.d0
      table_hqh(7,1,5)= -3.d0/35.d0
      table_hqh(7,1,6)= 2.d0/35.d0
      table_hqh(7,1,7)= 3.d0/35.d0
      table_hqh(7,1,8)= -3.d0/70.d0

      table_hqh(7,2,1)= 2.d0/5.d0
      table_hqh(7,2,2)= 3.d0/35.d0
      table_hqh(7,2,3)= 2.d0/5.d0
      table_hqh(7,2,4)= -3.d0/35.d0
      table_hqh(7,2,5)= -36.d0/35.d0
      table_hqh(7,2,6)= -4.d0/35.d0
      table_hqh(7,2,7)= 36.d0/35.d0
      table_hqh(7,2,8)= -4.d0/35.d0

      table_hqh(7,3,1)= -1.d0/70.d0
      table_hqh(7,3,2)= 0.d0
      table_hqh(7,3,3)= 4.d0/35.d0
      table_hqh(7,3,4)= -1.d0/70.d0
      table_hqh(7,3,5)= -3.d0/35.d0
      table_hqh(7,3,6)= -3.d0/70.d0
      table_hqh(7,3,7)= 3.d0/35.d0
      table_hqh(7,3,8)= 2.d0/35.d0

      table_hqh(7,4,1)= -53.d0/70.d0
      table_hqh(7,4,2)= -9.d0/70.d0
      table_hqh(7,4,3)= -17.d0/70.d0
      table_hqh(7,4,4)= 1.d0/14.d0
      table_hqh(7,4,5)= 6.d0/5.d0
      table_hqh(7,4,6)= -1.d0/10.d0
      table_hqh(7,4,7)= -6.d0/5.d0
      table_hqh(7,4,8)= 3.d0/10.d0

      table_hqh(7,5,1)= 18.d0/35.d0
      table_hqh(7,5,2)= 2.d0/35.d0
      table_hqh(7,5,3)= -18.d0/35.d0
      table_hqh(7,5,4)= 2.d0/35.d0
      table_hqh(7,5,5)= 0.d0
      table_hqh(7,5,6)= 2.d0/5.d0
      table_hqh(7,5,7)= 0.d0
      table_hqh(7,5,8)= -2.d0/5.d0

      table_hqh(7,6,1)= 17.d0/70.d0
      table_hqh(7,6,2)= 1.d0/14.d0
      table_hqh(7,6,3)= 53.d0/70.d0
      table_hqh(7,6,4)= -9.d0/70.d0
      table_hqh(7,6,5)= -6.d0/5.d0
      table_hqh(7,6,6)= -3.d0/10.d0
      table_hqh(7,6,7)= 6.d0/5.d0
      table_hqh(7,6,8)= 1.d0/10.d0


      table_hqh(8,1,1)= -1.d0/28.d0
      table_hqh(8,1,2)= -1.d0/210.d0
      table_hqh(8,1,3)= -1.d0/70.d0
      table_hqh(8,1,4)= 1.d0/420.d0
      table_hqh(8,1,5)= 3.d0/70.d0
      table_hqh(8,1,6)= -1.d0/84.d0
      table_hqh(8,1,7)= -3.d0/70.d0
      table_hqh(8,1,8)= 1.d0/210.d0

      table_hqh(8,2,1)= -8.d0/105.d0
      table_hqh(8,2,2)= -1.d0/70.d0
      table_hqh(8,2,3)= 1.d0/105.d0
      table_hqh(8,2,4)= 1.d0/210.d0
      table_hqh(8,2,5)= 4.d0/35.d0
      table_hqh(8,2,6)= -1.d0/105.d0
      table_hqh(8,2,7)= -4.d0/35.d0
      table_hqh(8,2,8)= 2.d0/35.d0

      table_hqh(8,3,1)= 1.d0/84.d0
      table_hqh(8,3,2)= 1.d0/420.d0
      table_hqh(8,3,3)= 11.d0/105.d0
      table_hqh(8,3,4)= -1.d0/140.d0
      table_hqh(8,3,5)= -2.d0/35.d0
      table_hqh(8,3,6)= -1.d0/84.d0
      table_hqh(8,3,7)= 2.d0/35.d0
      table_hqh(8,3,8)= 1.d0/14.d0

      table_hqh(8,4,1)= 41.d0/210.d0
      table_hqh(8,4,2)= 13.d0/420.d0
      table_hqh(8,4,3)= 29.d0/210.d0
      table_hqh(8,4,4)= -2.d0/105.d0
      table_hqh(8,4,5)= -3.d0/10.d0
      table_hqh(8,4,6)= 1.d0/30.d0
      table_hqh(8,4,7)= 3.d0/10.d0
      table_hqh(8,4,8)= 0.d0

      table_hqh(8,5,1)= -4.d0/21.d0
      table_hqh(8,5,2)= -1.d0/35.d0
      table_hqh(8,5,3)= -10.d0/21.d0
      table_hqh(8,5,4)= 4.d0/105.d0
      table_hqh(8,5,5)= 2.d0/5.d0
      table_hqh(8,5,6)= 0.d0
      table_hqh(8,5,7)= -2.d0/5.d0
      table_hqh(8,5,8)= -4.d0/15.d0

      table_hqh(8,6,1)= -1.d0/210.d0
      table_hqh(8,6,2)= -1.d0/420.d0
      table_hqh(8,6,3)= 71.d0/210.d0
      table_hqh(8,6,4)= -2.d0/105.d0
      table_hqh(8,6,5)= -1.d0/10.d0
      table_hqh(8,6,6)= -1.d0/30.d0
      table_hqh(8,6,7)= 1.d0/10.d0
      table_hqh(8,6,8)= 4.d0/15.d0

!------
      do k=1,8
         do j=1,6
            do i=1,8
               table_hhq(i,k,j)=table_hqh(i,j,k)
            enddo
         enddo
      enddo

!------
      table_hqq(1,1,1)= 13.d0/105.d0
      table_hqq(1,1,2)= 8.d0/105.d0
      table_hqq(1,1,3)= -1.d0/60.d0
      table_hqq(1,1,4)= -13.d0/28.d0
      table_hqq(1,1,5)= 59.d0/105.d0
      table_hqq(1,1,6)= -41.d0/420.d0

      table_hqq(1,2,1)= 8.d0/105.d0
      table_hqq(1,2,2)= 4.d0/15.d0
      table_hqq(1,2,3)= -1.d0/105.d0
      table_hqq(1,2,4)= -53.d0/105.d0
      table_hqq(1,2,5)= 12.d0/35.d0
      table_hqq(1,2,6)= 17.d0/105.d0

      table_hqq(1,3,1)= -1.d0/60.d0
      table_hqq(1,3,2)= -1.d0/105.d0
      table_hqq(1,3,3)= 1.d0/105.d0
      table_hqq(1,3,4)= 29.d0/420.d0
      table_hqq(1,3,5)= -11.d0/105.d0
      table_hqq(1,3,6)= 1.d0/28.d0

      table_hqq(1,4,1)= -13.d0/28.d0
      table_hqq(1,4,2)= -53.d0/105.d0
      table_hqq(1,4,3)= 29.d0/420.d0
      table_hqq(1,4,4)= 59.d0/30.d0
      table_hqq(1,4,5)= -32.d0/15.d0
      table_hqq(1,4,6)= 1.d0/6.d0

      table_hqq(1,5,1)= 59.d0/105.d0
      table_hqq(1,5,2)= 12.d0/35.d0
      table_hqq(1,5,3)= -11.d0/105.d0
      table_hqq(1,5,4)= -32.d0/15.d0
      table_hqq(1,5,5)= 8.d0/3.d0
      table_hqq(1,5,6)= -8.d0/15.d0

      table_hqq(1,6,1)= -41.d0/420.d0
      table_hqq(1,6,2)= 17.d0/105.d0
      table_hqq(1,6,3)= 1.d0/28.d0
      table_hqq(1,6,4)= 1.d0/6.d0
      table_hqq(1,6,5)= -8.d0/15.d0
      table_hqq(1,6,6)= 11.d0/30.d0


      table_hqq(2,1,1)= 1.d0/105.d0
      table_hqq(2,1,2)= 1.d0/105.d0
      table_hqq(2,1,3)= -1.d0/420.d0
      table_hqq(2,1,4)= -17.d0/420.d0
      table_hqq(2,1,5)= 1.d0/21.d0
      table_hqq(2,1,6)= -1.d0/140.d0

      table_hqq(2,2,1)= 1.d0/105.d0
      table_hqq(2,2,2)= 2.d0/35.d0
      table_hqq(2,2,3)= 0.d0
      table_hqq(2,2,4)= -3.d0/35.d0
      table_hqq(2,2,5)= 4.d0/105.d0
      table_hqq(2,2,6)= 1.d0/21.d0

      table_hqq(2,3,1)= -1.d0/420.d0
      table_hqq(2,3,2)= 0.d0
      table_hqq(2,3,3)= 1.d0/420.d0
      table_hqq(2,3,4)= 1.d0/105.d0
      table_hqq(2,3,5)= -2.d0/105.d0
      table_hqq(2,3,6)= 1.d0/105.d0

      table_hqq(2,4,1)= -17.d0/420.d0
      table_hqq(2,4,2)= -3.d0/35.d0
      table_hqq(2,4,3)= 1.d0/105.d0
      table_hqq(2,4,4)= 13.d0/60.d0
      table_hqq(2,4,5)= -1.d0/5.d0
      table_hqq(2,4,6)= -1.d0/60.d0

      table_hqq(2,5,1)= 1.d0/21.d0
      table_hqq(2,5,2)= 4.d0/105.d0
      table_hqq(2,5,3)= -2.d0/105.d0
      table_hqq(2,5,4)= -1.d0/5.d0
      table_hqq(2,5,5)= 4.d0/15.d0
      table_hqq(2,5,6)= -1.d0/15.d0

      table_hqq(2,6,1)= -1.d0/140.d0
      table_hqq(2,6,2)= 1.d0/21.d0
      table_hqq(2,6,3)= 1.d0/105.d0
      table_hqq(2,6,4)= -1.d0/60.d0
      table_hqq(2,6,5)= -1.d0/15.d0
      table_hqq(2,6,6)= 1.d0/12.d0


      table_hqq(3,1,1)= 1.d0/105.d0
      table_hqq(3,1,2)= -1.d0/105.d0
      table_hqq(3,1,3)= -1.d0/60.d0
      table_hqq(3,1,4)= -1.d0/28.d0
      table_hqq(3,1,5)= 11.d0/105.d0
      table_hqq(3,1,6)= -29.d0/420.d0

      table_hqq(3,2,1)= -1.d0/105.d0
      table_hqq(3,2,2)= 4.d0/15.d0
      table_hqq(3,2,3)= 8.d0/105.d0
      table_hqq(3,2,4)= -17.d0/105.d0
      table_hqq(3,2,5)= -12.d0/35.d0
      table_hqq(3,2,6)= 53.d0/105.d0

      table_hqq(3,3,1)= -1.d0/60.d0
      table_hqq(3,3,2)= 8.d0/105.d0
      table_hqq(3,3,3)= 13.d0/105.d0
      table_hqq(3,3,4)= 41.d0/420.d0
      table_hqq(3,3,5)= -59.d0/105.d0
      table_hqq(3,3,6)= 13.d0/28.d0

      table_hqq(3,4,1)= -1.d0/28.d0
      table_hqq(3,4,2)= -17.d0/105.d0
      table_hqq(3,4,3)= 41.d0/420.d0
      table_hqq(3,4,4)= 11.d0/30.d0
      table_hqq(3,4,5)= -8.d0/15.d0
      table_hqq(3,4,6)= 1.d0/6.d0

      table_hqq(3,5,1)= 11.d0/105.d0
      table_hqq(3,5,2)= -12.d0/35.d0
      table_hqq(3,5,3)= -59.d0/105.d0
      table_hqq(3,5,4)= -8.d0/15.d0
      table_hqq(3,5,5)= 8.d0/3.d0
      table_hqq(3,5,6)= -32.d0/15.d0

      table_hqq(3,6,1)= -29.d0/420.d0
      table_hqq(3,6,2)= 53.d0/105.d0
      table_hqq(3,6,3)= 13.d0/28.d0
      table_hqq(3,6,4)= 1.d0/6.d0
      table_hqq(3,6,5)= -32.d0/15.d0
      table_hqq(3,6,6)= 59.d0/30.d0


      table_hqq(4,1,1)= -1.d0/420.d0
      table_hqq(4,1,2)= 0.d0
      table_hqq(4,1,3)= 1.d0/420.d0
      table_hqq(4,1,4)= 1.d0/105.d0
      table_hqq(4,1,5)= -2.d0/105.d0
      table_hqq(4,1,6)= 1.d0/105.d0

      table_hqq(4,2,1)= 0.d0
      table_hqq(4,2,2)= -2.d0/35.d0
      table_hqq(4,2,3)= -1.d0/105.d0
      table_hqq(4,2,4)= 1.d0/21.d0
      table_hqq(4,2,5)= 4.d0/105.d0
      table_hqq(4,2,6)= -3.d0/35.d0

      table_hqq(4,3,1)= 1.d0/420.d0
      table_hqq(4,3,2)= -1.d0/105.d0
      table_hqq(4,3,3)= -1.d0/105.d0
      table_hqq(4,3,4)= -1.d0/140.d0
      table_hqq(4,3,5)= 1.d0/21.d0
      table_hqq(4,3,6)= -17.d0/420.d0

      table_hqq(4,4,1)= 1.d0/105.d0
      table_hqq(4,4,2)= 1.d0/21.d0
      table_hqq(4,4,3)= -1.d0/140.d0
      table_hqq(4,4,4)= -1.d0/12.d0
      table_hqq(4,4,5)= 1.d0/15.d0
      table_hqq(4,4,6)= 1.d0/60.d0

      table_hqq(4,5,1)= -2.d0/105.d0
      table_hqq(4,5,2)= 4.d0/105.d0
      table_hqq(4,5,3)= 1.d0/21.d0
      table_hqq(4,5,4)= 1.d0/15.d0
      table_hqq(4,5,5)= -4.d0/15.d0
      table_hqq(4,5,6)= 1.d0/5.d0

      table_hqq(4,6,1)= 1.d0/105.d0
      table_hqq(4,6,2)= -3.d0/35.d0
      table_hqq(4,6,3)= -17.d0/420.d0
      table_hqq(4,6,4)= 1.d0/60.d0
      table_hqq(4,6,5)= 1.d0/5.d0
      table_hqq(4,6,6)= -13.d0/60.d0


      table_hqq(5,1,1)= -1.d0/14.d0
      table_hqq(5,1,2)= -2.d0/35.d0
      table_hqq(5,1,3)= 1.d0/35.d0
      table_hqq(5,1,4)= 3.d0/10.d0
      table_hqq(5,1,5)= -2.d0/5.d0
      table_hqq(5,1,6)= 1.d0/10.d0

      table_hqq(5,2,1)= -2.d0/35.d0
      table_hqq(5,2,2)= -24.d0/35.d0
      table_hqq(5,2,3)= -2.d0/35.d0
      table_hqq(5,2,4)= 4.d0/5.d0
      table_hqq(5,2,5)= 0.d0
      table_hqq(5,2,6)= -4.d0/5.d0

      table_hqq(5,3,1)= 1.d0/35.d0
      table_hqq(5,3,2)= -2.d0/35.d0
      table_hqq(5,3,3)= -1.d0/14.d0
      table_hqq(5,3,4)= -1.d0/10.d0
      table_hqq(5,3,5)= 2.d0/5.d0
      table_hqq(5,3,6)= -3.d0/10.d0

      table_hqq(5,4,1)= 3.d0/10.d0
      table_hqq(5,4,2)= 4.d0/5.d0
      table_hqq(5,4,3)= -1.d0/10.d0
      table_hqq(5,4,4)= -9.d0/5.d0
      table_hqq(5,4,5)= 8.d0/5.d0
      table_hqq(5,4,6)= 1.d0/5.d0

      table_hqq(5,5,1)= -2.d0/5.d0
      table_hqq(5,5,2)= 0.d0
      table_hqq(5,5,3)= 2.d0/5.d0
      table_hqq(5,5,4)= 8.d0/5.d0
      table_hqq(5,5,5)= -16.d0/5.d0
      table_hqq(5,5,6)= 8.d0/5.d0

      table_hqq(5,6,1)= 1.d0/10.d0
      table_hqq(5,6,2)= -4.d0/5.d0
      table_hqq(5,6,3)= -3.d0/10.d0
      table_hqq(5,6,4)= 1.d0/5.d0
      table_hqq(5,6,5)= 8.d0/5.d0
      table_hqq(5,6,6)= -9.d0/5.d0


      table_hqq(6,1,1)= 17.d0/210.d0
      table_hqq(6,1,2)= 4.d0/105.d0
      table_hqq(6,1,3)= -1.d0/420.d0
      table_hqq(6,1,4)= -17.d0/60.d0
      table_hqq(6,1,5)= 1.d0/3.d0
      table_hqq(6,1,6)= -1.d0/20.d0

      table_hqq(6,2,1)= 4.d0/105.d0
      table_hqq(6,2,2)= -8.d0/105.d0
      table_hqq(6,2,3)= -1.d0/35.d0
      table_hqq(6,2,4)= -1.d0/15.d0
      table_hqq(6,2,5)= 4.d0/15.d0
      table_hqq(6,2,6)= -1.d0/5.d0

      table_hqq(6,3,1)= -1.d0/420.d0
      table_hqq(6,3,2)= -1.d0/35.d0
      table_hqq(6,3,3)= -2.d0/105.d0
      table_hqq(6,3,4)= 1.d0/60.d0
      table_hqq(6,3,5)= 1.d0/15.d0
      table_hqq(6,3,6)= -1.d0/12.d0

      table_hqq(6,4,1)= -17.d0/60.d0
      table_hqq(6,4,2)= -1.d0/15.d0
      table_hqq(6,4,3)= 1.d0/60.d0
      table_hqq(6,4,4)= 14.d0/15.d0
      table_hqq(6,4,5)= -6.d0/5.d0
      table_hqq(6,4,6)= 4.d0/15.d0

      table_hqq(6,5,1)= 1.d0/3.d0
      table_hqq(6,5,2)= 4.d0/15.d0
      table_hqq(6,5,3)= 1.d0/15.d0
      table_hqq(6,5,4)= -6.d0/5.d0
      table_hqq(6,5,5)= 16.d0/15.d0
      table_hqq(6,5,6)= 2.d0/15.d0

      table_hqq(6,6,1)= -1.d0/20.d0
      table_hqq(6,6,2)= -1.d0/5.d0
      table_hqq(6,6,3)= -1.d0/12.d0
      table_hqq(6,6,4)= 4.d0/15.d0
      table_hqq(6,6,5)= 2.d0/15.d0
      table_hqq(6,6,6)= -2.d0/5.d0


      table_hqq(7,1,1)= 1.d0/14.d0
      table_hqq(7,1,2)= 2.d0/35.d0
      table_hqq(7,1,3)= -1.d0/35.d0
      table_hqq(7,1,4)= -3.d0/10.d0
      table_hqq(7,1,5)= 2.d0/5.d0
      table_hqq(7,1,6)= -1.d0/10.d0

      table_hqq(7,2,1)= 2.d0/35.d0
      table_hqq(7,2,2)= 24.d0/35.d0
      table_hqq(7,2,3)= 2.d0/35.d0
      table_hqq(7,2,4)= -4.d0/5.d0
      table_hqq(7,2,5)= 0.d0
      table_hqq(7,2,6)= 4.d0/5.d0

      table_hqq(7,3,1)= -1.d0/35.d0
      table_hqq(7,3,2)= 2.d0/35.d0
      table_hqq(7,3,3)= 1.d0/14.d0
      table_hqq(7,3,4)= 1.d0/10.d0
      table_hqq(7,3,5)= -2.d0/5.d0
      table_hqq(7,3,6)= 3.d0/10.d0

      table_hqq(7,4,1)= -3.d0/10.d0
      table_hqq(7,4,2)= -4.d0/5.d0
      table_hqq(7,4,3)= 1.d0/10.d0
      table_hqq(7,4,4)= 9.d0/5.d0
      table_hqq(7,4,5)= -8.d0/5.d0
      table_hqq(7,4,6)= -1.d0/5.d0

      table_hqq(7,5,1)= 2.d0/5.d0
      table_hqq(7,5,2)= 0.d0
      table_hqq(7,5,3)= -2.d0/5.d0
      table_hqq(7,5,4)= -8.d0/5.d0
      table_hqq(7,5,5)= 16.d0/5.d0
      table_hqq(7,5,6)= -8.d0/5.d0

      table_hqq(7,6,1)= -1.d0/10.d0
      table_hqq(7,6,2)= 4.d0/5.d0
      table_hqq(7,6,3)= 3.d0/10.d0
      table_hqq(7,6,4)= -1.d0/5.d0
      table_hqq(7,6,5)= -8.d0/5.d0
      table_hqq(7,6,6)= 9.d0/5.d0


      table_hqq(8,1,1)= -2.d0/105.d0
      table_hqq(8,1,2)= -1.d0/35.d0
      table_hqq(8,1,3)= -1.d0/420.d0
      table_hqq(8,1,4)= 1.d0/12.d0
      table_hqq(8,1,5)= -1.d0/15.d0
      table_hqq(8,1,6)= -1.d0/60.d0

      table_hqq(8,2,1)= -1.d0/35.d0
      table_hqq(8,2,2)= -8.d0/105.d0
      table_hqq(8,2,3)= 4.d0/105.d0
      table_hqq(8,2,4)= 1.d0/5.d0
      table_hqq(8,2,5)= -4.d0/15.d0
      table_hqq(8,2,6)= 1.d0/15.d0

      table_hqq(8,3,1)= -1.d0/420.d0
      table_hqq(8,3,2)= 4.d0/105.d0
      table_hqq(8,3,3)= 17.d0/210.d0
      table_hqq(8,3,4)= 1.d0/20.d0
      table_hqq(8,3,5)= -1.d0/3.d0
      table_hqq(8,3,6)= 17.d0/60.d0

      table_hqq(8,4,1)= 1.d0/12.d0
      table_hqq(8,4,2)= 1.d0/5.d0
      table_hqq(8,4,3)= 1.d0/20.d0
      table_hqq(8,4,4)= -2.d0/5.d0
      table_hqq(8,4,5)= 2.d0/15.d0
      table_hqq(8,4,6)= 4.d0/15.d0

      table_hqq(8,5,1)= -1.d0/15.d0
      table_hqq(8,5,2)= -4.d0/15.d0
      table_hqq(8,5,3)= -1.d0/3.d0
      table_hqq(8,5,4)= 2.d0/15.d0
      table_hqq(8,5,5)= 16.d0/15.d0
      table_hqq(8,5,6)= -6.d0/5.d0

      table_hqq(8,6,1)= -1.d0/60.d0
      table_hqq(8,6,2)= 1.d0/15.d0
      table_hqq(8,6,3)= 17.d0/60.d0
      table_hqq(8,6,4)= 4.d0/15.d0
      table_hqq(8,6,5)= -6.d0/5.d0
      table_hqq(8,6,6)= 14.d0/15.d0


      table_lq(1,1)= 1.d0/6.d0
      table_lq(1,2)= 1.d0/3.d0
      table_lq(1,3)= 0.d0
      table_lq(1,4)=-5.d0/6.d0
      table_lq(1,5)= 2.d0/3.d0
      table_lq(1,6)= 1.d0/6.d0

      table_lq(2,1)= 0.d0
      table_lq(2,2)= 1.d0/3.d0
      table_lq(2,3)= 1.d0/6.d0
      table_lq(2,4)=-1.d0/6.d0
      table_lq(2,5)=-2.d0/3.d0
      table_lq(2,6)= 5.d0/6.d0

      table_lq(3,1)=-1.d0/6.d0
      table_lq(3,2)=-2.d0/3.d0
      table_lq(3,3)=-1.d0/6.d0
      table_lq(3,4)= 1.d0
      table_lq(3,5)= 0.d0
      table_lq(3,6)=-1.d0

      table_lq(4,1)= 1.d0/6.d0
      table_lq(4,2)= 2.d0/3.d0
      table_lq(4,3)= 1.d0/6.d0
      table_lq(4,4)=-1.d0
      table_lq(4,5)= 0.d0
      table_lq(4,6)= 1.d0

      do j=1,4
         do i=1,6
            table_ql(i,j)=table_lq(j,i)
         enddo
      enddo

      table_lql(1,1,1)= 3.d0/20.d0
      table_lql(1,1,2)= 1.d0/60.d0
      table_lql(1,1,3)=-1.d0/6.d0
      table_lql(1,1,4)= 1.d0/6.d0

      table_lql(1,2,1)= 1.d0/5.d0
      table_lql(1,2,2)= 2.d0/15.d0
      table_lql(1,2,3)=-1.d0/3.d0
      table_lql(1,2,4)= 1.d0/3.d0

      table_lql(1,3,1)=-1.d0/60.d0
      table_lql(1,3,2)= 1.d0/60.d0
      table_lql(1,3,3)= 0.d0
      table_lql(1,3,4)= 0.d0

      table_lql(1,4,1)=-2.d0/3.d0
      table_lql(1,4,2)=-1.d0/6.d0
      table_lql(1,4,3)= 5.d0/6.d0
      table_lql(1,4,4)=-5.d0/6.d0

      table_lql(1,5,1)= 2.d0/3.d0
      table_lql(1,5,2)= 0.d0
      table_lql(1,5,3)=-2.d0/3.d0
      table_lql(1,5,4)= 2.d0/3.d0

      table_lql(1,6,1)= 0.d0
      table_lql(1,6,2)= 1.d0/6.d0
      table_lql(1,6,3)=-1.d0/6.d0
      table_lql(1,6,4)= 1.d0/6.d0

      table_lql(2,1,1)= 1.d0/60.d0
      table_lql(2,1,2)=-1.d0/60.d0
      table_lql(2,1,3)= 0.d0
      table_lql(2,1,4)= 0.d0

      table_lql(2,2,1)= 2.d0/15.d0
      table_lql(2,2,2)= 1.d0/5.d0
      table_lql(2,2,3)=-1.d0/3.d0
      table_lql(2,2,4)= 1.d0/3.d0

      table_lql(2,3,1)= 1.d0/60.d0
      table_lql(2,3,2)= 3.d0/20.d0
      table_lql(2,3,3)=-1.d0/6.d0
      table_lql(2,3,4)= 1.d0/6.d0

      table_lql(2,4,1)=-1.d0/6.d0
      table_lql(2,4,2)= 0.d0
      table_lql(2,4,3)= 1.d0/6.d0
      table_lql(2,4,4)=-1.d0/6.d0

      table_lql(2,5,1)= 0.d0
      table_lql(2,5,2)=-2.d0/3.d0
      table_lql(2,5,3)= 2.d0/3.d0
      table_lql(2,5,4)=-2.d0/3.d0

      table_lql(2,6,1)= 1.d0/6.d0
      table_lql(2,6,2)= 2.d0/3.d0
      table_lql(2,6,3)=-5.d0/6.d0
      table_lql(2,6,4)= 5.d0/6.d0

      table_lql(3,1,1)=-1.d0/6.d0
      table_lql(3,1,2)= 0.d0
      table_lql(3,1,3)= 1.d0/6.d0
      table_lql(3,1,4)=-1.d0/6.d0

      table_lql(3,2,1)=-1.d0/3.d0
      table_lql(3,2,2)=-1.d0/3.d0
      table_lql(3,2,3)= 2.d0/3.d0
      table_lql(3,2,4)=-2.d0/3.d0

      table_lql(3,3,1)= 0.d0
      table_lql(3,3,2)=-1.d0/6.d0
      table_lql(3,3,3)= 1.d0/6.d0
      table_lql(3,3,4)=-1.d0/6.d0

      table_lql(3,4,1)= 5.d0/6.d0
      table_lql(3,4,2)= 1.d0/6.d0
      table_lql(3,4,3)=-1.d0
      table_lql(3,4,4)= 1.d0

      table_lql(3,5,1)=-2.d0/3.d0
      table_lql(3,5,2)= 2.d0/3.d0
      table_lql(3,5,3)= 0.d0
      table_lql(3,5,4)= 0.d0

      table_lql(3,6,1)=-1.d0/6.d0
      table_lql(3,6,2)=-5.d0/6.d0
      table_lql(3,6,3)= 1.d0
      table_lql(3,6,4)=-1.d0

      table_lql(4,1,1)= 1.d0/6.d0
      table_lql(4,1,2)= 0.d0
      table_lql(4,1,3)=-1.d0/6.d0
      table_lql(4,1,4)= 1.d0/6.d0

      table_lql(4,2,1)= 1.d0/3.d0
      table_lql(4,2,2)= 1.d0/3.d0
      table_lql(4,2,3)=-2.d0/3.d0
      table_lql(4,2,4)= 2.d0/3.d0

      table_lql(4,3,1)= 0.d0
      table_lql(4,3,2)= 1.d0/6.d0
      table_lql(4,3,3)=-1.d0/6.d0
      table_lql(4,3,4)= 1.d0/6.d0

      table_lql(4,4,1)=-5.d0/6.d0
      table_lql(4,4,2)=-1.d0/6.d0
      table_lql(4,4,3)= 1.d0
      table_lql(4,4,4)=-1.d0

      table_lql(4,5,1)= 2.d0/3.d0
      table_lql(4,5,2)=-2.d0/3.d0
      table_lql(4,5,3)= 0.d0
      table_lql(4,5,4)= 0.d0

      table_lql(4,6,1)= 1.d0/6.d0
      table_lql(4,6,2)= 5.d0/6.d0
      table_lql(4,6,3)=-1.d0
      table_lql(4,6,4)= 1.d0

      do k=1,6
         do j=1,4
            do i=1,4
               table_llq(i,j,k)=table_lql(i,k,j)
            enddo
         enddo
      enddo

      table_lqq(1,1,1)= 7.d0/60.d0
      table_lqq(1,1,2)= 1.d0/15.d0
      table_lqq(1,1,3)=-1.d0/60.d0
      table_lqq(1,1,4)=-13.d0/30.d0
      table_lqq(1,1,5)= 8.d0/15.d0
      table_lqq(1,1,6)=-1.d0/10.d0

      table_lqq(1,2,1)= 1.d0/15.d0
      table_lqq(1,2,2)= 4.d0/15.d0
      table_lqq(1,2,3)= 0.d0
      table_lqq(1,2,4)=-7.d0/15.d0
      table_lqq(1,2,5)= 4.d0/15.d0
      table_lqq(1,2,6)= 1.d0/5.d0

      table_lqq(1,3,1)=-1.d0/60.d0
      table_lqq(1,3,2)= 0.d0
      table_lqq(1,3,3)= 1.d0/60.d0
      table_lqq(1,3,4)= 1.d0/15.d0
      table_lqq(1,3,5)=-2.d0/15.d0
      table_lqq(1,3,6)= 1.d0/15.d0

      table_lqq(1,4,1)=-13.d0/30.d0
      table_lqq(1,4,2)=-7.d0/15.d0
      table_lqq(1,4,3)= 1.d0/15.d0
      table_lqq(1,4,4)= 11.d0/6.d0
      table_lqq(1,4,5)=-2.d0
      table_lqq(1,4,6)= 1.d0/6.d0

      table_lqq(1,5,1)= 8.d0/15.d0
      table_lqq(1,5,2)= 4.d0/15.d0
      table_lqq(1,5,3)=-2.d0/15.d0
      table_lqq(1,5,4)=-2.d0
      table_lqq(1,5,5)= 8.d0/3.d0
      table_lqq(1,5,6)=-2.d0/3.d0

      table_lqq(1,6,1)=-1.d0/10.d0
      table_lqq(1,6,2)= 1.d0/5.d0
      table_lqq(1,6,3)= 1.d0/15.d0
      table_lqq(1,6,4)= 1.d0/6.d0
      table_lqq(1,6,5)=-2.d0/3.d0
      table_lqq(1,6,6)= 1.d0/2.d0

      table_lqq(2,1,1)= 1.d0/60.d0
      table_lqq(2,1,2)= 0.d0
      table_lqq(2,1,3)=-1.d0/60.d0
      table_lqq(2,1,4)=-1.d0/15.d0
      table_lqq(2,1,5)= 2.d0/15.d0
      table_lqq(2,1,6)=-1.d0/15.d0

      table_lqq(2,2,1)= 0.d0
      table_lqq(2,2,2)= 4.d0/15.d0
      table_lqq(2,2,3)= 1.d0/15.d0
      table_lqq(2,2,4)=-1.d0/5.d0
      table_lqq(2,2,5)=-4.d0/15.d0
      table_lqq(2,2,6)= 7.d0/15.d0

      table_lqq(2,3,1)=-1.d0/60.d0
      table_lqq(2,3,2)= 1.d0/15.d0
      table_lqq(2,3,3)= 7.d0/60.d0
      table_lqq(2,3,4)= 1.d0/10.d0
      table_lqq(2,3,5)=-8.d0/15.d0
      table_lqq(2,3,6)= 13.d0/30.d0

      table_lqq(2,4,1)=-1.d0/15.d0
      table_lqq(2,4,2)=-1.d0/5.d0
      table_lqq(2,4,3)= 1.d0/10.d0
      table_lqq(2,4,4)= 1.d0/2.d0
      table_lqq(2,4,5)=-2.d0/3.d0
      table_lqq(2,4,6)= 1.d0/6.d0

      table_lqq(2,5,1)= 2.d0/15.d0
      table_lqq(2,5,2)=-4.d0/15.d0
      table_lqq(2,5,3)=-8.d0/15.d0
      table_lqq(2,5,4)=-2.d0/3.d0
      table_lqq(2,5,5)= 8.d0/3.d0
      table_lqq(2,5,6)=-2.d0

      table_lqq(2,6,1)=-1.d0/15.d0
      table_lqq(2,6,2)= 7.d0/15.d0
      table_lqq(2,6,3)= 13.d0/30.d0
      table_lqq(2,6,4)= 1.d0/6.d0
      table_lqq(2,6,5)=-2.d0
      table_lqq(2,6,6)= 11.d0/6.d0

      table_lqq(3,1,1)=-2.d0/15.d0
      table_lqq(3,1,2)=-1.d0/15.d0
      table_lqq(3,1,3)= 1.d0/30.d0
      table_lqq(3,1,4)= 1.d0/2.d0
      table_lqq(3,1,5)=-2.d0/3.d0
      table_lqq(3,1,6)= 1.d0/6.d0

      table_lqq(3,2,1)=-1.d0/15.d0
      table_lqq(3,2,2)=-8.d0/15.d0
      table_lqq(3,2,3)=-1.d0/15.d0
      table_lqq(3,2,4)= 2.d0/3.d0
      table_lqq(3,2,5)= 0.d0
      table_lqq(3,2,6)=-2.d0/3.d0

      table_lqq(3,3,1)= 1.d0/30.d0
      table_lqq(3,3,2)=-1.d0/15.d0
      table_lqq(3,3,3)=-2.d0/15.d0
      table_lqq(3,3,4)=-1.d0/6.d0
      table_lqq(3,3,5)= 2.d0/3.d0
      table_lqq(3,3,6)=-1.d0/2.d0

      table_lqq(3,4,1)= 1.d0/2.d0
      table_lqq(3,4,2)= 2.d0/3.d0
      table_lqq(3,4,3)=-1.d0/6.d0
      table_lqq(3,4,4)=-7.d0/3.d0
      table_lqq(3,4,5)= 8.d0/3.d0
      table_lqq(3,4,6)=-1.d0/3.d0

      table_lqq(3,5,1)=-2.d0/3.d0
      table_lqq(3,5,2)= 0.d0
      table_lqq(3,5,3)= 2.d0/3.d0
      table_lqq(3,5,4)= 8.d0/3.d0
      table_lqq(3,5,5)=-16.d0/3.d0
      table_lqq(3,5,6)= 8.d0/3.d0

      table_lqq(3,6,1)= 1.d0/6.d0
      table_lqq(3,6,2)=-2.d0/3.d0
      table_lqq(3,6,3)=-1.d0/2.d0
      table_lqq(3,6,4)=-1.d0/3.d0
      table_lqq(3,6,5)= 8.d0/3.d0
      table_lqq(3,6,6)=-7.d0/3.d0

      table_lqq(4,1,1)= 2.d0/15.d0
      table_lqq(4,1,2)= 1.d0/15.d0
      table_lqq(4,1,3)=-1.d0/30.d0
      table_lqq(4,1,4)=-1.d0/2.d0
      table_lqq(4,1,5)= 2.d0/3.d0
      table_lqq(4,1,6)=-1.d0/6.d0

      table_lqq(4,2,1)= 1.d0/15.d0
      table_lqq(4,2,2)= 8.d0/15.d0
      table_lqq(4,2,3)= 1.d0/15.d0
      table_lqq(4,2,4)=-2.d0/3.d0
      table_lqq(4,2,5)= 0.d0
      table_lqq(4,2,6)= 2.d0/3.d0

      table_lqq(4,3,1)=-1.d0/30.d0
      table_lqq(4,3,2)= 1.d0/15.d0
      table_lqq(4,3,3)= 2.d0/15.d0
      table_lqq(4,3,4)= 1.d0/6.d0
      table_lqq(4,3,5)=-2.d0/3.d0
      table_lqq(4,3,6)= 1.d0/2.d0

      table_lqq(4,4,1)=-1.d0/2.d0
      table_lqq(4,4,2)=-2.d0/3.d0
      table_lqq(4,4,3)= 1.d0/6.d0
      table_lqq(4,4,4)= 7.d0/3.d0
      table_lqq(4,4,5)=-8.d0/3.d0
      table_lqq(4,4,6)= 1.d0/3.d0

      table_lqq(4,5,1)= 2.d0/3.d0
      table_lqq(4,5,2)= 0.d0
      table_lqq(4,5,3)=-2.d0/3.d0
      table_lqq(4,5,4)=-8.d0/3.d0
      table_lqq(4,5,5)= 16.d0/3.d0
      table_lqq(4,5,6)=-8.d0/3.d0

      table_lqq(4,6,1)=-1.d0/6.d0
      table_lqq(4,6,2)= 2.d0/3.d0
      table_lqq(4,6,3)= 1.d0/2.d0
      table_lqq(4,6,4)= 1.d0/3.d0
      table_lqq(4,6,5)=-8.d0/3.d0
      table_lqq(4,6,6)= 7.d0/3.d0
      
      return
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

END MODULE libfem
