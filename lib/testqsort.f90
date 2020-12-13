program test_qsort
  USE task_kinds,ONLY: dp
  use libqsort
  implicit none
  integer, parameter :: r = 10
  integer, dimension(1:r) :: array = &        ! (1:r)
     (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
  complex(dp), dimension(1:r) :: arrayx = &        ! (1:r)
     (/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),(4.D0,0.D0), &
       (5.d0,0.d0),(6.d0,0.d0),(7.d0,0.d0),(8.d0,0.d0),(9.D0,0.D0)/)
  print *, "array is ", array
  call qsort_ic(array,arrayx)
  print *, "sorted array is ", array
  print *, "sorted associated array is ", arrayx
end program test_qsort
