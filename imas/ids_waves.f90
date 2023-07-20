! ids_waves.f90

MODULE ids_waves

  USE task_kinds
  
  PRIVATE
  PUBLIC get_ids_waves
  PUBLIC clear_ids_waves
  PUBLIC put_ids_waves

CONTAINS

  SUBROUTINE get_ids_waves

    USE ids_schemas,ONLY: ids_waves
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_waves):: waves_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: nwave_max,nbeam_max,nmode_max
    INTEGER:: nwave,nbeam,nmode
    REAL(rkind):: freq,ang_t,rk0,rn_tor

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx,'public',trim(machine),'3')
    call ids_get(idx,'waves',waves_ids)
    call imas_close(idx)
    
    ! --- set size of coherent_wave ---

    nwave_max = size(waves_ids%coherent_wave)

    nbeam_max=0
    nmode_max=0
    DO nwave=1,nwave_max
       IF(waves_ids%coherent_wave(nwave)%wave_solver_type%index.EQ.1) &
            nbeam_max=nbeam_max+1
       IF(waves_ids%coherent_wave(nwave)%wave_solver_type%index.EQ.2) &
            nmode_max=nmode_max+1
    END DO

  END SUBROUTINE get_ids_waves

  SUBROUTINE clean_ids_waves
    IMPLICIT NONE
  END SUBROUTINE clean_ids_waves

  ! ----------------------------
  ! SAVE RESULTS INTO WAVES IDS
  ! ----------------------------

  SUBROUTINE put_ids_waves

    USE ids_schemas,ONLY: ids_waves,ids_ec_launchers,ids_equilibrium
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_waves):: waves_ids
    type(ids_equilibrium):: equilibrium_ids
    type(ids_ec_launchers):: ec_launchers_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: nwave_max,nbeam_max,nmode_max
    INTEGER:: nwave,nbeam,nmode

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx,'public',trim(machine),'3')
    call ids_get(idx,'equilibrium',equilibrium_ids)
    call ids_get(idx,'ec_launchers',ec_launchers_ids)
    call ids_get(idx,'waves',waves_ids)
    call imas_close(idx)

    ! GLOBAL QUANTITIES
    waves_ids%ids_properties%homogeneous_time = 1
    allocate(waves_ids%time(1))
    waves_ids%time = equilibrium_ids%time
    allocate(waves_ids%code%name(1))
    waves_ids%code%name = 'TASK-WR'

    ! --- set number of beam and allocate ---

    nbeam_max = nray_max
    allocate(waves_ids%coherent_wave(nbeam_max))

    ! --- loop over nbeam: coherent_wave structure ---

    do nbeam=1,nbeam_max
       ! --- identifier ---
       ALLOCATE(waves_ids%coherent_wave(nbeam)%identifier%type%name(1))
                waves_ids%coherent_wave(nbeam)%identifier%type%name='EC'
                waves_ids%coherent_wave(nbeam)%identifier%type%index=1
       ALLOCATE(waves_ids%coherent_wave(nbeam)%identifier%antenna_name(1))
                waves_ids%coherent_wave(nbeam)%identifier%antenna_name &
                     = ec_launchers_ids%beam(nbeam)%name
                waves_ids%coherent_wave(nbeam)%identifier%index_in_antenna=1
       ! --- wave_solver_type ---
       ALLOCATE(waves_ids%coherent_wave(nbeam)%wave_solver_type%name(1))
                waves_ids%coherent_wave(nbeam)%wave_solver_type%name='WR'
                waves_ids%coherent_wave(nbeam)%wave_solver_type%index=1
                
       ! --- global quantities ---
       ntime_max=1
       DO ntime=1,ntime_max         
          waves_ids%coherent_wave(nbeam)%global_quantities(ntime)%frequency &
               = ec_launchers_ids%beam(nbeam)%frequency%data(ntime)
          ALLOCATE(waves_ids%coherent_wave(nbeam)%global_quantities(ntime) &
               %n_tor(ntim))
          freq=ec_launchers_ids%beam(nbeam)%frequency%data(ntime)
          ang_t=ec_launchers_ids%beam(nbeam)%steering_angle_tor(ntime)
          rk0=2.D0*PI*freq/VC*SIN(ang_t)
          rn_tor=rk0*ec_launchers_ids%beam(nbeam)%launching_position%r(ntime)
          waves_ids%coherent_wave(nbeam)%global_quantities(ntime)%n_tor(1) &
               =NINT(rn_tor)
          
       DO nstep=1,nstep_max         
       allocate(waves_ids%coherent_wave(nbeam)%global_quantities(1))
       allocate(waves_ids%coherent_wave(nbeam)%identifier%antenna_name(1))
       waves_ids%coherent_wave(nbeam)%identifier%antenna_name &
            = ec_launchers_ids%beam(nbeam)%name
       waves_ids%coherent_wave(nbeam)%global_quantities(1)%frequency &
            = ec_launchers_ids%beam(nbeam)%frequency%data(1)
       waves_ids%coherent_wave(nbeam)%global_quantities(1)%electrons &
            %power_thermal = extrascal(nbeam,2)
       waves_ids%coherent_wave(nbeam)%global_quantities(1)%power &
            = extrascal(nbeam,2)
       waves_ids%coherent_wave(nbeam)%global_quantities(1)%current_tor &
            = extrascal(nbeam,3)
       waves_ids%coherent_wave(nbeam)%identifier%type%description &
            = 'TASK-WR'
       ! if statement needed to filter empty sources for hcd2core_sources
       if(ec_launchers_ids%beam(nbeam)%power_launched%data(1)>0) then
          waves_ids%coherent_wave(nbeam)%identifier%type%name = 'EC'
          waves_ids%coherent_wave(nbeam)%identifier%type%index = 1
          waves_ids%coherent_wave(nbeam)%wave_solver_type%index = 1
       endif
       
       nstep_max=nstpmax_nray(nbeam)
       ALLOCATE(waves_ids%coherent_wave(nbeam)%
       
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %grid%psi(npnt))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %grid%rho_tor(npnt))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %grid%rho_tor_norm(npnt))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %power_density(npnt))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %electrons%power_density_thermal(npnt))
       allocate(waves_ids%coherent_wave(nbeam)%profiles_1d(1) &
            %current_parallel_density(npnt))
       waves_ids%coherent_wave(nbeam)%profiles_1d(1)%grid%psi  &
            = profout(nbeam,1,1:npnt)**2*(psiedge-psiax)+psiax
       waves_ids%coherent_wave(nbeam)%profiles_1d(1)%power_density &
            = 1.e6*profout(nbeam,2,1:npnt)
       waves_ids%coherent_wave(nbeam)%profiles_1d(1)%electrons% &
            power_density_thermal = 1.e6*profout(nbeam,2,1:npnt)
       waves_ids%coherent_wave(nbeam)%profiles_1d(1)%current_parallel_density &
            = -1.e6*profout(nbeam,3,1:npnt)*iplasma/abs(iplasma)
    enddo ! LOOP OVER BEAMS (LAUNCHERS)

    waves_ids%code%output_flag = 0 ! NO ERROR

222 format(6(1x,1pe13.5))
225 format(3(1x,1pe13.5))
999 format(1p,3e13.5,i5,e13.5)

    RETURN
  END SUBROUTINE put_ids_waves
END MODULE ids_waves
