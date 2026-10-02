program xtr
  USE bpsd
  USE equnit_mod
  USE equunit_mod
  USE trunit
  USE trcomm, ONLY: NT
  USE tr_bpsd,ONLY: tr_bpsd_set, tr_bpsd_get
  IMPLICIT NONE
  INTEGER(ikind):: ierr
  
  CALL eq_init
  CALL equ_init
  CALL tr_init

  CALL eq_parm(1,'eqparm',ierr)
  CALL equ_parm(1,'equparm',ierr)
  CALL tr_parm(1,'trparm',ierr)

  CALL tr_view

  CALL bpsd_load(ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'XX bpsd_load error: ierr=',ierr
     STOP
  ENDIF

  CALL tr_bpsd_get(ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'XX tr_bpsd_get error: ierr=',ierr
     STOP
  ENDIF

  NT=0
  CALL TRLOOP

  CALL bpsd_save(ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'XX bpsd_save error: ierr=',ierr
     STOP
  ENDIF
  STOP
END PROGRAM xtr

