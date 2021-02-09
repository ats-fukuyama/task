! sakmain.f90

PROGRAM sak

  USE sakcomm
  USE sak1
  USE sak2
  USE sak3
  IMPLICIT NONE
  INTEGER:: id
  EXTERNAL GSOPEN,GSCLOS

  WRITE(6,*) '## TASK/SAK 2020/11/21'
  CALL GSOPEN

  id=3
  SELECT CASE(id)
  CASE(1)
     CALL sak_1
  CASE(2)
     CALL sak_2
  CASE(3)
     CALL sak_3
  END SELECT
  
  CALL GSCLOS

  STOP
END PROGRAM sak
