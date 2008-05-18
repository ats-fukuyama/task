PROGRAM main

  USE hfreya_all
  USE mcnbi_all
  IMPLICIT NONE

  CHARACTER(LEN=4):: argc

  CALL GETARG(0,argc)
  IF(argc =="")THEN
    CALL GETARG(2,argc)
  ELSE 
    CALL GETARG(1,argc)
  END IF

  OPEN(5,FILE='hfreya.in', ACTION='READ')
  OPEN(6,FILE='hfreya.out',ACTION='WRITE')

  CALL hfreya

  CLOSE(5)
  CLOSE(6)

  OPEN(5,FILE='mcnbi.in', ACTION='READ')
  OPEN(6,FILE='mcnbi.out',ACTION='WRITE')
  REWIND(9)
  REWIND(10)
  REWIND(15)

  CALL mcnbi

  CLOSE(5)
  CLOSE(6)

  OPEN(5,FILE='fit.in', ACTION='READ')
  OPEN(6,FILE='fit.out',ACTION='WRITE')
  REWIND(20)
  CALL fit(20)

  IF(argc == "pb4") THEN
    REWIND(5)
    REWIND(30)
    CALL fit(30)

    REWIND(5)
    REWIND(40)
    CALL fit(40)

    REWIND(60)
    REWIND(70)
    REWIND(80)
    CALL depsum
  END IF

  CLOSE(5)
  CLOSE(6)

  STOP
END PROGRAM
