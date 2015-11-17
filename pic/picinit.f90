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
      nx = 128
      ny = 128
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
      bxmin  =    0.0d0     !: max value of background x magnetic flux density
      bxmax  =    0.0d0     !: min value of background x magnetic flux density
      bymin  =    0.0d0     !: min value of background x magnetic flux density
      bymax  =    0.0d0     !: min value of background x magnetic flux density
      bzmin  =    0.0d0     !: min value of background x magnetic flux density
      bzmax  =    0.0d0     !: z component of background magnetic flux density
      c      =    1.0d0
      omega  =    1.0d0
      eps = 0.0000000000000001d0    !: constants to define boundary condition
      jxbg   =    0.0d0    ! x component of background current density 
      jybg   =    0.0d0    ! y component of background current density
      jzbg   =    0.0d0    ! z component of background current density
      f      =    0.0d0    ! frequency of background current density
      thetax =    0.0d0    ! first phase of jxbg
      thetay =    0.0d0    ! first phase of jybg
      thetaz =    0.0d0    ! first phase of jzbg
!---------------------------------------------------------

    RETURN
  END SUBROUTINE pic_init
END MODULE picinit
