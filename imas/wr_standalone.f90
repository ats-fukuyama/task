! task_wr_standalone.f90

PROGRAM task_wr_standalone

  USE ids_schemas, ONLY: &
       equilibrium_ids, &
       core_profiles_ids, &
       ec_launchers_ids, &
       waves_ids
  USE f90_file_reader, ONLY: file2buffer
  USE codeparam_task_wr
  IMPLICIT NONE
  TYPE(ids_parameter_input):: codeparam_task_wr
  INTEGER:: iounit=1

  ! READ XML INPUT FILE
  CALL file2buffer('input/task_wr.xml', &
       iounit,codeparam_task_wr%parameters_value)

  ! EXECUTE TASK_WR
  CALL task_wr_sub(
