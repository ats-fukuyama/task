!
!               ############# TASK/PIC #############
!
!                Particle-in-cell component of TASK
!                Based on pic2des_p by Prof. Hiroshi Naitou
!
!-----------------------------------------------------------------------

PROGRAM pic_main
  USE piccomm,ONLY: nrank,nsize
  USE picinit,ONLY: pic_init
  USE picparm,ONLY: pic_parm
  USE picmenu,ONLY: pic_menu
  USE libmtx

  IMPLICIT none
  INCLUDE 'mpif.h'
  INTEGER :: ierr

  CALL mtx_initialize
  write(6,*) 'nsize,nrank,ncomm=',nsize,nrank,ncomm
  IF(nrank.EQ.9) THEN
     WRITE(6,*) '## TASK/PIC 2015/11/22'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  END IF

  CALL pic_init
  IF(nrank.EQ.0) THEN
     CALL pic_parm(1,'picparm.nl',ierr)
  END IF
  CALL mtx_broadcast1_integer(ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: xxparm.nl'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL pic_menu

  CLOSE(7)
  IF(nrank.EQ.0) CALL GSCLOS
  CALL mtx_finalize
  CALL pic_deallocate
  STOP
END PROGRAM pic_main
