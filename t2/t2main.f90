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
!C    *        KYOTO 615-8530, JAPAN                                *
!C    *                                                             *
!C    *  PROGRAM VERSION                                            *
!C    *                                                             *
!C    *                                                             *
!C    ***************************************************************
PROGRAM T2


  USE T2COMM, ONLY: T2COMM_DEALLOCATE
  USE T2INIT, ONLY: T2_INIT
  USE T2NGRA, ONLY: T2_NGRA
  
  USE T2VGRA, ONLY: T2_VGRA
  USE T2LOOP, ONLY: T2_LOOP
  USE T2INTG, ONLY: T2_INTG
  IMPLICIT NONE
  
  WRITE(6,*)'***** TASK/T2 V0.00 2013-MM-DD            *****'

  WRITE(6,*)'INITIALIZE PROGRAM'
  CALL T2_INIT       !C INITIALIZE PROGRAM
  WRITE(6,*)'GENERATE   NODE     GRAPH'
  CALL T2_NGRA       !C GENERATE   NODE     GRAPH
  WRITE(6,*)'GENERATE   VARIABLE GRAPH'
  CALL T2_VGRA       !C GENERATE   VARIABLE GRAPH
  WRITE(6,*)'CALCULATE INTEGRAND ARRAYS'
  CALL T2_INTG
  WRITE(6,*)'CALCULATE TEMPORAL EVOLUTION'
  CALL T2_LOOP
  WRITE(6,*)'TERMINATE  PROGRAM'
  CALL T2COMM_DEALLOCATE !C TERMINATE  PROGRAM

  WRITE(6,*)'***** TASK/T2 WILL NOW TERMINATE NORMALLY *****'
  STOP

END PROGRAM T2
  
