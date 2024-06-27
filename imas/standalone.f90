program standalone

  use ids_schemas, only: ids_equilibrium
  use ids_routines, only: imas_open_env,imas_create_env,imas_close,ids_get,ids_put, ids_deallocate
  use physics_module_level1
  
  implicit none

  integer:: idx,error_flag
  type(ids_equilibrium):: equilibrium_in,equilibrium_out
  character(len=200):: user,machine
  character(len=:), pointer:: error_message

  ! DEFINE LOCAL DATABASE
  call getenv('USER',user)
  machine = 'ITER'

  ! OPEN INPUT DATAFILE FROM OFFICIAL IMAS SCENARIO DATABASE
  write(*,*) '=> Read input IDSs'
  call imas_open_env('ids',131024,41,idx,'public',trim(machine),'3')
  call ids_get(idx,'equilibrium',equilibrium_in)
  call imas_close(idx)
  write(*,*) 'Finished reading input IDSs'

  ! EXECUTE PHYSICS CODE
  call physics_model_level1(equilibrium_in,equilibrium_out,error_flag,error_message)

  if(error_flag.eq.0) then

     ! EXPORT RESULTS TO LOCAL DATABASE
     write(*,*) '=> Export output IDSs to local database'
     call imas_create_env('ids',131024,2,0,0,idx,trim(user),trim(machine),'3')
     call ids_put(idx,'equilibrium',equilibrium_out)
     call imas_close(idx)
     write(*,*) 'Done exporting.'
     write(*,*) ' '
     write(*,*) 'End of standalone'

  else

     write(*,*) error_message
     write(*,*) '=> Program stopped.'

  endif

  ! RELEASING THE MEMORY ALLOCATED FOR EACH IDS
  call ids_deallocate(equilibrium_in)
  call ids_deallocate(equilibrium_out)
  
  
end program standalone
