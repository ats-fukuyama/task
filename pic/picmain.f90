!
!               ############# TASK/PIC #############
!
!                Particle-in-cell component of TASK
!                Based on pic2des_p by Prof. Hiroshi Naitou
!
!-----------------------------------------------------------------------

PROGRAM pic_main
  USE piccomm,ONLY: nrank,nsize,ncomm,pic_deallocate
  USE picinit,ONLY: pic_init
  USE picparm,ONLY: pic_parm
  USE picmenu,ONLY: pic_menu
  USE libmtx
  USE libmpi

  IMPLICIT none
  INCLUDE 'mpif.h'
  INTEGER :: ierr

  CALL mtx_initialize
  write(6,*) 'rank,nsize,ncomm=',nrank,nsize,ncomm
  IF(nrank.EQ.0) THEN
     WRITE(6,*) '## TASK/PIC 2015/12/08'
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

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
     CLOSE(7)
  END IF
  CALL mtx_finalize
  STOP
END PROGRAM pic_main
