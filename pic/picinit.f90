!  ***** TASK/PIC INIT *****

MODULE picinit

PUBLIC pic_init

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE pic_init

    USE piccomm_parm
    IMPLICIT NONE

!.......................................................................
!............. npx   : number of particles in x                  .......
!............. npy   : number of particles in y                  .......
!............. np = npx * npy : total number of each particles   .......
!............. nx    : number of grids in x                      .......
!............. ny    : number of grids in y                      .......
!............. iend  : total time steps                          .......
!............. nhmod : output energies in each nhmod steps       .......
!.......................................................................

      npx = 100
      npy = 100
      npz = 100
      nx = 128
      ny = 128
      nz = 10
      iend = 1000
      nhmod = 1

!----- set some parameters -------------------------------
      dt     =    0.2d0     !: time step size
      me     =    1.0d0     !: electron mass
      mi     = 1836.0d0     !: ion mass
      chrge  =   -1.0d0     !: electron charge
      chrgi  =    1.0d0     !: ion charge
      te     =    1.0d0     !: electron temperature
      ti     =    1.0d0     !: ion temperature
!      ez     =    0.0d0     !: z electric field
      bx     =    0.0d0     !: x magnetic flux density
      by     =    0.0d0     !: y magnetic flux density
      bz     =    0.0d0     !: z magnetic flux density
      c      =    1.0d0
      omega  =    1.0d0
      eps = 0.0000000000000001d0    !: constants to define boundary condition
!---------------------------------------------------------

    RETURN
  END SUBROUTINE pic_init
END MODULE picinit
