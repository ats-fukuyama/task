! sakmain.f90

PROGRAM sak

  USE sakcomm
  USE sak1
  USE sak2
  USE sak3
  USE sak4
  USE sak5
  IMPLICIT NONE
  INTEGER:: id
  EXTERNAL GSOPEN,GSCLOS

  WRITE(6,*) '## TASK/SAK 2021/02/25'
  CALL GSOPEN

  id=5
1 CONTINUE
  WRITE(6,'(A)') '## INPUT id [1..5]:'
  READ(5,*,ERR=1,END=9000) id
  SELECT CASE(id)
  CASE(1)
     CALL sak_1
  CASE(2)
     CALL sak_2
  CASE(3)
     CALL sak_3
  CASE(4)
     CALL sak_4
  CASE(5)
     CALL sak_5
  END SELECT

9000 CONTINUE
  CALL GSCLOS

  STOP
END PROGRAM sak
