c
      program main
      use equunit_mod
      use trunit_mod
      use eqinp_mod

      integer ii

c
                         write(6,*)'==step:0'
      call equ_init
      call tr_init
      call eqinp
c
                         write(6,*)'==step:1'
      call equ_prof
      call tr_prof
                         write(6,*)'==step:2'
      call equ_calc
c
      do ii=1,10
                         write(6,*)'==step:3'
      call tr_exec
c
      call equ_calc
      enddo
c
      stop
      end program main

