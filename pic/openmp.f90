program openmp
  implicit none
  integer :: i
!$omp parallel do
do i = 1,5
 write(*,*) "Hello OPENMPI world"
end do
!$omp end parallel do

end program
