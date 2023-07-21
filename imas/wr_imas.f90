MODULE task_wr_sub

CONTAINS

  ! IMAS module of task/wr

  SUBROUTINE task_wr_ids( &
       equilibrium_ids, &
       core_profiles_ids, &
       ec_launchers_ids, &
       waves_ids, &
       codeparam_task_wr, &
       output_flag,output_message)
      
    USE ids_schemas, ONLY: &
         equilibrium_ids, &
         core_profiles_ids, &
         ec_launchers_ids, &
         waves_ids
    USE ids_routines
    USE codeparam_input_task_wr
    USE wrcomm
    USE wrsetup
    USE wrexec
    IMPLICIT NONE
    TYPE(ids_equilibrium),INTENT(IN):: equilibrium_in
    TYPE(ids_core_profiles),INTENT(IN):: core_profile_in
    TYPE(ids_ec_launchers),INTENT(IN):: ec_launchers_in
    TYPE(ids_waves),INTENT(IN):: waves_out
    TYPE(ids_parameters_input):: codeparam_task_wr
    TYPE(type_codeparam_data):: codeparam_data_task_wr
    LOGICAL,INTENT(OUT):: output_flag
    CHARACTER{LEN=128),INTENT(OUT):: output_messave

    CALL wr_init
      
    CALL assign_codeparm_task_wr( &
         codeparam_task_wr%parameters_value, &
         codeparam_data_task_wr)

    CALL wr_allocate

    CALL get_ids_equilibrium(equilibrium_in)
    CALL get_ids_core_profile(core_profile_in)
    CALL get_ec_launchers(ec_launchers_in)

    CALL wr_setup
    CALL wr_exec

    CALL put_waves(waves_out)

    CALL wr_deallocate

    RETURN
  END SUBROUTINE task_wr_ids
END MODULE task_wr_sub
