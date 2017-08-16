!
! program to convert OPEN-ADAS ADF11 files (Unresolved only) to TASK data files
!
PROGRAM conv_adf11

  USE ADF11
  IMPLICIT NONE
  INTEGER:: IERR,IZ0,ICLASS,IS,ND
  REAL(dp):: PN,PT,DR

  CALL READ_ADF11(IERR)
  IF(IERR.NE.0) THEN
     WRITE(6,*) 'XX CCONV_ADF11: READ_ADF11'
     STOP
  END IF

  CALL SAVE_ADF11_bin(IERR)
  IF(IERR.NE.0) THEN
     WRITE(6,*) 'XX CCONV_ADF11: SAVE_ADF11_bin'
     STOP
  END IF

  STOP
END PROGRAM conv_adf11

