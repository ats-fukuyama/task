!  ***** TASK/PIC INIT *****

MODULE picinit

PUBLIC pic_init

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE pic_init

    USE piccomm_parm
    IMPLICIT NONE

!.......................................................................
!............. npxmax: number of particles in x                  .......
!............. npymax: number of particles in y                  .......
!............. np = npx * npy : total number of each particles   .......
!............. nxmax : number of grids in x                      .......
!............. nymax : number of grids in y                      .......
!............. ntmax : total time steps                          .......
!............. ntstep: short summary step intercal               .......
!.......................................................................

      npxmax = 100
      npymax = 100
      nxmax = 128
      nymax = 128
      ntmax = 1000
      ntstep= 1

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
      bymin  =    0.0d0     !: min value of background y magnetic flux density
      bymax  =    0.0d0     !: max value of background y magnetic flux density
      bzmin  =    0.0d0     !: min value of background z magnetic flux density
      bzmax  =    0.0d0     !: max value of background z magnetic flux density
      vcfact =   10.0d0     !: c^2/omegape^2 debye^2
      eps    =   1.D-16     !: constants to define boundary condition
      omega  =    3.0d0     !: antena frequency
      jxant  =    0.0d0     ! x component of antenna current density 
      jyant  =    0.0d0     ! y component of antenna current density
      jzant  =    0.0d0     ! z component of antenn current density
      phxant =    0.0d0     ! initial phase of jxant
      phyant =    0.0d0     ! initial phase of jyant
      phzant =    0.0d0     ! initial phase of jzant
!---------------------------------------------------------

    RETURN
  END SUBROUTINE pic_init
END MODULE picinit
