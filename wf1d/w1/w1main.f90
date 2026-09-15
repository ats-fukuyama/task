!
!     ************************* TASK.W1 *************************
!
!     ICRF WAVE PROPAGATION AND ABSORPTION
!
!          DIFFERENTIAL OR INTEGRO-DIFFERENTIAL ANALYSIS
!          SLAB PLASMA MODEL
!          MULTIPLE LAYER OR FINITE ELEMENT METHOD
!
!     PROGRAMED BY
!
!          A. FUKUYAMA, S. NISHIYAMA, T. MATSUISHI AND k. HAMAMATSU(*)
!
!          DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
!          FACULTY OF ENGINEERING
!          OKAYAMA UNIVERSITY
!          OKAYA 700, JAPAN
!
!      (*) DEPARTMENT OF LARGE TOKAMAK RESEARCH
!          NAKA FUSION RESEARCH ESTABLISHMENT
!          JAPAN ATOMIC ENERGY RESEARCH INSTITUTE
!          NAKA, IBARAKI 311-01, JAPAN
!
!     PROGRAM VERSION
!
!          V1.68 : 90/02/17
!          V1.70 : 90/03/14
!          V1.71 : 90/06/08 : ICDTYP ADDED
!          V1.72 : 90/06/16 : WVYSIZ ADDED
!          V1.73 : 90/06/18 : W1CLCD CORRECTED, WVYSIZ CORRECTED
!          V2.00 : 90/09/08 : B, A, HELICITY, ABSORBED BOUNDARY ADDED
!          V2.01 : 90/09/12 : HELICITY CURRENT ADDED
!          V2.02 : 90/09/14 : CORRECTION
!          V2.03 : 90/09/20 : TALPHA,TAU CORRECTED
!          V2.04 : 90/09/22 : NEW GRAPH, QH RECORRECTED
!          V2.05 : 90/10/16 : TRAPPED PARTICLE EFFECTS
!          V2.06 : 91/03/29 : NONRESONANT FORCE MODIFIED
!          V2.10 : 91/06/18 : UNDERFLOW ADJUSTMENT FOR SX
!          V2.11 : 92/03/18 : WORKING ON HP WS
!          V2.12 : 94/01/12 : SPLIT SOURCE FILES, INCLUDE FILE, ERF
!
!     *************************************************************

PROGRAM w1_main
  USE w1comm
  USE w1init
  USE w1parm
  USE w1menu

  IMPLICIT NONE
  INTEGER(ikind)  :: ierr

  CALL GSOPEN
  WRITE(6,*) '## TASK/W1 2018/07/01'
  OPEN(7,STATUS='SCRATCH')

  CALL w1_init
  CALL w1_parm(1,'w1parm',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: w1parm'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL w1_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM w1_main

