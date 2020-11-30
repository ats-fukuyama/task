! task_kinds.f90

module task_kinds
  implicit none
  integer, parameter :: rkind=selected_real_kind(12,100)
  integer, parameter :: ikind=selected_int_kind(8)
  integer, parameter :: dp=rkind
  integer, parameter :: sp=selected_real_kind(6)
end module task_kinds
