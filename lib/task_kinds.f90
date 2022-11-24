! task_kinds.f90

module task_kinds
  USE ISO_FORTRAN_ENV
  integer, parameter :: qkind=REAL128
  integer, parameter :: rkind=REAL64
  integer, parameter :: ikind=INT32
  integer, parameter :: dp=REAL64
  integer, parameter :: sp=REAL32
  integer, parameter :: long=INT64
!  integer, parameter :: rkind=selected_real_kind(12,100)
!  integer, parameter :: ikind=selected_int_kind(8)
end module task_kinds
