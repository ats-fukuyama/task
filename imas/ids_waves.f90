! ids_waves.f90

MODULE ids_waves

  USE task_kinds
  INTEGER:: nrgmax,nzgmax
  REAL(rkind):: b0,rip,psip0,psipa
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: rg,zg
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: br,bt,bz,psip2d
  
  PRIVATE
  PUBLIC get_ids_waves
  PUBLIC clear_ids_waves

CONTAINS

  SUBROUTINE get_ids_waves

    USE ids_schemas,ONLY: ids_waves
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_waves):: equilibrium_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: nrgmax,nzgmax

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx, &
         'public',trim(machine),'3')
    call ids_get(idx,'equilibrium',equilibrium_ids)
    call imas_close(idx)
    
    ! --- set size of R and Z grid ---

    nrgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1)
    nzgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2)

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(zg)) DEALLOCATE(zg)
    IF(ALLOCATED(br)) DEALLOCATE(br)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    IF(ALLOCATED(psip2d)) DEALLOCATE(psip2d)
    ALLOCATE(RG(nrgmax))
    ALLOCATE(ZG(nzgmax))
    allocate(br(nrgmax,nzgmax))
    allocate(bt(nrgmax,nzgmax))
    allocate(bz(nrgmax,nzgmax))
    allocate(psip2d(nrgmax,nzgmax))
    rg    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1
    zg    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2
    br    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_r
    bt    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_tor
    bz    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_z
    psip2d= equilibrium_ids%time_slice(1)%profiles_2d(1)%psi
    psipa = equilibrium_ids%time_slice(1)%global_quantities%psi_boundary
    psip0 = equilibrium_ids%time_slice(1)%global_quantities%psi_axis
    b0    = -equilibrium_ids%time_slice(1)%global_quantities%magnetic_axis%b_field_tor
    rip   = equilibrium_ids%time_slice(1)%global_quantities%ip
  END SUBROUTINE get_ids_waves

  SUBROUTINE clean_ids_waves
    IMPLICIT NONE
    deallocate(psip2d)
    deallocate(rg,zg)
    deallocate(br)
    deallocate(bt)
    deallocate(bz)
  END SUBROUTINE clean_ids_waves
END MODULE ids_waves








  ! ----------------------------
  ! SAVE RESULTS INTO WAVES IDS
  ! ----------------------------

  ! GLOBAL QUANTITIES
  waves_ids%ids_properties%homogeneous_time = 1
  allocate(waves_ids%time(1))
  waves_ids%time = equilibrium_ids%time
  allocate(waves_ids%code%name(1))
  waves_ids%code%name = 'TORBEAM'

  ! LOOP OVER BEAMS (LAUNCHERS)
  !if(nbeam.gt.10) nbeam = 10 ! MSR waiting for IMAS-3271
  allocate(waves_ids%coherent_wave(nbeam))
  do ibeam=1,nbeam
     allocate(waves_ids%coherent_wave(ibeam)%global_quantities(1))
     allocate(waves_ids%coherent_wave(ibeam)%identifier%antenna_name(1))
     waves_ids%coherent_wave(ibeam)%identifier%antenna_name                      = ec_launchers_ids%launcher(ibeam)%name
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%frequency               = ec_launchers_ids%launcher(ibeam)%frequency%data(1)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%electrons%power_thermal = extrascal(ibeam,2)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%power                   = extrascal(ibeam,2)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%current_tor             = extrascal(ibeam,3)
     waves_ids%coherent_wave(ibeam)%identifier%type%description                  = 'TORBEAM'
     ! if statement needed to filter empty sources for hcd2core_sources
     if(ec_launchers_ids%launcher(ibeam)%power_launched%data(1)>0) then
        waves_ids%coherent_wave(ibeam)%identifier%type%name                         = 'EC'
        waves_ids%coherent_wave(ibeam)%identifier%type%index                        = 1
        waves_ids%coherent_wave(ibeam)%wave_solver_type%index                       = 1 ! BEAM/RAY TRACING
     endif
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor_norm(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%power_density(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%electrons%power_density_thermal(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%current_parallel_density(npnt))
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi                        = profout(ibeam,1,1:npnt)**2*(psiedge-psiax)+psiax
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%power_density                   = 1.e6*profout(ibeam,2,1:npnt)
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%electrons%power_density_thermal = 1.e6*profout(ibeam,2,1:npnt)
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%current_parallel_density        = -1.e6*profout(ibeam,3,1:npnt)*iplasma/abs(iplasma)

     ! rho_tor and rho_tor_norm to be filled for hcd2core_sources not to complain
     call interpos(-equilibrium_ids%time_slice(1)%profiles_1d%psi, &
          equilibrium_ids%time_slice(1)%profiles_1d%rho_tor, &
          size(equilibrium_ids%time_slice(1)%profiles_1d%psi), &
          size(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi), &
          xout=-waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi, &
          yout=waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor)
     call interpos(-equilibrium_ids%time_slice(1)%profiles_1d%psi, &
          equilibrium_ids%time_slice(1)%profiles_1d%rho_tor_norm, &
          size(equilibrium_ids%time_slice(1)%profiles_1d%psi), &
          size(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi), &
          xout=-waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi, &
          yout=waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor_norm)

     ! LOOP OVER RAYS
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(5))
     iend = min(npointsout(ibeam),ntraj)
     if(ec_launchers_ids%launcher(ibeam)%power_launched%data(1)>0) then
        do iray = 1,5 ! 5 rays (to not mix with input beams)
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%r(iend))
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%z(iend))
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%phi(iend))
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%r   = 1.e-2*trajout(ibeam,1+3*(iray-1),1:iend)
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%z   = 1.e-2*trajout(ibeam,2+3*(iray-1),1:iend)
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%phi = trajout(ibeam,3+3*(iray-1),1:iend)
        enddo
     endif

     ! SOME QUANTITIES ARE KNOWN ONLY ON THE CENTRAL RAY
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%electrons%power(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%parallel(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%perpendicular(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing%beam(1)%length(iend))
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_initial                 = extrascal(ibeam,1)*1.e6
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%electrons%power               = extradata(ibeam,1,1:iend)*1.e6
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%parallel      = extradata(ibeam,2,1:iend)
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%perpendicular = extradata(ibeam,3,1:iend)
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%length                        = extradata(ibeam,4,1:iend)

  enddo ! LOOP OVER BEAMS (LAUNCHERS)

  waves_ids%code%output_flag = 0 ! NO ERROR

222 format(6(1x,1pe13.5))
225 format(3(1x,1pe13.5))
999 format(1p,3e13.5,i5,e13.5)

