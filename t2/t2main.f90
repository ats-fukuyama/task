!C    *************************** TASK.T2 ***************************
!C    *                                                              * 
!C    *  Diffusive-Advective 2D Transport in a Tokamak              *
!C    *                                                             *
!C    *  PROGRAMMED BY                                              *
!C    *                                                             *
!C    *        H. SETO                                              *
!C    *        DEPARTMENT OF NUCLEAR ENGINEERING                    *
!C    *        GRADUATE SCHOOL OF ENGINEERING                       *
!C    *        KYOTO UNIVERSITY                                     *
!C    *        KYOTO 615-8540, JAPAN                                *
!C    *                                                             *
!C    *  PROGRAM VERSION                                            *
!C    *                                                             *
!C    *                                                             *
!C    ***************************************************************

PROGRAM T2

  USE T2CNST, ONLY: i0ikind
  USE T2INIT, ONLY: T2_INIT
  USE T2PARM, ONLY: T2_PARM
  USE T2MENU, ONLY: T2_MENU
    USE LIBMTX, ONLY: &
         MTX_INITIALIZE,MTX_FINALIZE  
  IMPLICIT NONE
  INTEGER(i0ikind)  :: i0err
  
  CALL GSOPEN
  CALL MTX_INITIALIZE

  WRITE(6,*)'***** TASK/T2 V0.01 2013-12-10 *****'
  OPEN(7,STATUS='SCRATCH')

  CALL T2_INIT                       !C Initialize program
  CALL T2_PARM(1,'t2parm.nl',i0err)  !C Read parameter file
  IF(i0err /= 0 .AND. i0err /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: t2parm.nl'
     WRITE(6,*) '     i0err = ',i0err
     STOP
  END IF

  CALL T2_MENU

  CLOSE(7)

  CALL MTX_FINALIZE
  CALL GSCLOS
  

  STOP
END PROGRAM T2
