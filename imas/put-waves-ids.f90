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

