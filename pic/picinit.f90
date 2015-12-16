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
!............. ntstep: line output interval                      .......
!............. ntgstep: global date save interval                .......
!............. ntpstep: profile date save interval               .......
!............. npomax : number of particles to follow orbits     .......
!............. npostep: particle data save internal (1,1+npostep,...)...
!............. ntostep: orbit data save interval                 .......
!.......................................................................

      npxmax = 100
      npymax = 100
      nxmax = 64
      nymax = 64
      ntmax = 1000
      ntstep= 10
      ntgstep= 1
      ntpstep= 100
      npomax = 0
      npostep= 1
      ntostep= 1

!----- set some parameters -------------------------------

      dt     =    0.02d0    !: time step size
      me     =    1.0d0     !: electron mass
      mi     = 1836.0d0     !: ion mass
      chrge  =   -1.0d0     !: electron charge
      chrgi  =    1.0d0     !: ion charge
      te     =    1.0d0     !: electron temperature
      ti     =    1.0d0     !: ion temperature
      densx  =    0.0d0     !: x density gradient
      bxmin  =    0.0d0     !: max value of background x magnetic flux density
      bxmax  =    0.0d0     !: min value of background x magnetic flux density
      bymin  =    0.0d0     !: min value of background y magnetic flux density
      bymax  =    0.0d0     !: max value of background y magnetic flux density
      bzmin  =    0.0d0     !: min value of background z magnetic flux density
      bzmax  =    0.0d0     !: max value of background z magnetic flux density
      vcfact =   10.0d0     !: c^2/omegape^2 debye^2
      eps    =   1.D-16     !: constants to define boundary condition
      omega  =    0.0d0     !: antena frequency
      jxant  =    0.0d0     ! x component of antenna current density 
      jyant  =    0.0d0     ! y component of antenna current density
      jzant  =    0.0d0     ! z component of antenn current density
      phxant =    0.0d0     ! initial phase of jxant
      phyant =    0.0d0     ! initial phase of jyant
      phzant =    0.0d0     ! initial phase of jzant

!---------------------------------------------------------

      model_push = 15       ! force to push particles
                            !   0: no force
                            !   1: electrostatic E field
                            !   2: electromagnetic E field
                            !   4: static magnetic B field
                            !   8: electromagnetic B field
                            !  15: all field
      model_boundary = 0    ! boundary condition 
                            !   0: periodic boundary
                            !   1: conducting wall and particle reflection
                            !   2: absorbing boundary and particle reflection
      model_antenna = 0     ! antenna location
                            !   0: near xmin uniform
                            !   1: near xmax uniform
                            !   2: near ymin uniform
                            !   3: near ymax uniform
      model_wg = 0          ! waveguide location
                            !   0: near xmin uniform
                            !   1: near xmax uniform
                            !   2: near ymin uniform
                            !   3: near ymax uniform
      xmin_wg = 0.D0
      xmax_wg = 0.D0
      ymin_wg = 0.D0
      ymax_wg = 0.D0
      amp_wg  = 0.D0        ! amplitude of wave vector potential
      ph_wg   = 0.D0        ! phase difference between x/ymin and  x/ymax 2*pi
      rot_wg  = 0.D0        ! angle of wave electric fields: 0 for y-direction
                            !                               90 for z-direction
      eli_wg  = 0.D0        ! angle of elipticity of wave electric field

      model_matrix0 = 0     ! ksp default
      model_matrix1 = 4     ! ksp default
      model_matrix2 = 5     ! ksp default
      tolerance_matrix=1.D-7! tolerance for ksp

    RETURN
  END SUBROUTINE pic_init
END MODULE picinit
