  ! COUND NUMBER OF INPUT BEAMS (EVEN THOSE WITH NO POWER)
  if (associated(ec_launchers_ids%launcher)) then
     nbeam = size(ec_launchers_ids%launcher)
  else
     nbeam = 0
  endif

  ! ALLOCATIONS
  allocate(profout(nbeam,3,npnt))
  allocate(trajout(nbeam,15,ntraj))
  allocate(npointsout(nbeam))
  allocate(extrascal(nbeam,3))      ! increase the second dimension if there are more data to be saved
  allocate(extradata(nbeam,4,ndat)) ! increase the second dimension if there are more data to be saved
  profout    = 0.
  trajout    = 0.
  npointsout = 0
  extrascal  = 0.
  extradata  = 0.

