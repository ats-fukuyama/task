!
!               ############# TASK/PIC #############
!
!                Particle-in-cell component of TASK
!                Based on pic2des_p by Prof. Hiroshi Naitou
!
!-----------------------------------------------------------------------

PROGRAM pic_main
  USE piccomm,ONLY: myid,nodes
  USE picinit,ONLY: pic_init
  USE picparm,ONLY: pic_parm
  USE picmenu,ONLY: pic_menu
  USE libmpi

  IMPLICIT none
  INCLUDE 'mpif.h'
  INTEGER :: ierr

  CALL GSOPEN
  WRITE(6,*) '## TASK/PIC 2015/05/16'
  OPEN(7,STATUS='SCRATCH')
  
!----- start parallel calculation
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nodes,ierr)
     if( myid .eq. 0 ) write(*,*) '*** number of nodes = ***', nodes
  CALL mtx_mpi(mpi_comm_world,myid,nodes)

  CALL pic_init
  CALL pic_parm(1,'picparm.nl',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: xxparm.nl'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL pic_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM pic_main
