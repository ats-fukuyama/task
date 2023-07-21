  ! READ INPUT XML FILE
  call parse_torbeam_codeparam(codeparam_torbeam%parameters_value,codeparam_torbeam_data)

  ! Initialize antenna data as TORBEAM input:

  !Int parameters
  npow         = codeparam_torbeam_data%npow
  ncd          = codeparam_torbeam_data%ncd
  ncdroutine   = codeparam_torbeam_data%ncdroutine
  nprofv       = codeparam_torbeam_data%nprofv
  noout        = codeparam_torbeam_data%noout
  nrela        = codeparam_torbeam_data%nrela
  nmaxh        = codeparam_torbeam_data%nmaxh
  nabsroutine  = codeparam_torbeam_data%nabsroutine
  nastra       = codeparam_torbeam_data%nastra
  nprofcalc    = codeparam_torbeam_data%nprofcalc
  ncdharm      = codeparam_torbeam_data%ncdharm
  npnts_extrap = codeparam_torbeam_data%npnts_extrap
  nfreq_extrap = codeparam_torbeam_data%nfreq_extrap
  nrel         = codeparam_torbeam_data%nrel

  !Float parameters
  xrtol   = codeparam_torbeam_data%xrtol
  xatol   = codeparam_torbeam_data%xatol
  xstep   = codeparam_torbeam_data%xstep
  rhostop = codeparam_torbeam_data%rhostop
  xzsrch  = codeparam_torbeam_data%xzsrch

