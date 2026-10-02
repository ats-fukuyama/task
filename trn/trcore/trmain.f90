!     *************************** TASK.TR ***************************
!     *                                                             *
!     *  Diffusive 1D Transport in a Tokamak                        *
!     *                                                             *
!     *  PROGRAMMED BY                                              *
!     *                                                             *
!     *        A. FUKUYAMA                                          *
!     *        DEPARTMENT OF NUCLEAR ENGINEERING                    *
!     *        GRADUATE SCHOOL OF ENGINEERING                       *
!     *        KYOTO UNIVERSITY                                     *
!     *        KYOTO 606-8501, JAPAN                                *
!     *                                                             *
!     *  ACKNOWLEDGEMENTS FOR                                       *
!     *                                                             *
!     *        T. KASAI                                             *
!     *        M. ISHIZU                                            *
!     *        S. TAKATSUKA                                         *
!     *        S. TAKATSUKA                                         *
!     *        M. HONDA                                             *
!     *        T. IKARI                                             *
!     *                                                             *
!     *  PROGRAM VERSION                                            *
!     *                                                             *
!     *        V1.00 : 89/04/26                                     *
!     *        V1.40 : 91/07/24                                     *
!     *        V2.00 : 93/06/18                                     *
!     *        V2.21 : 00/04/19                                     *
!     *        V3.00 : 11/11/20                                     *
!     *                                                             *
!     ***************************************************************

PROGRAM tr  

  USE trcomm, ONLY : rkind,ikind
  USE equnit_mod
!  USE equunit_mod
  USE trinit
  USE trmenu
  IMPLICIT NONE
  INTEGER(ikind):: ierr

!     ------ INITIALIZATION ------

  WRITE(6,*) '***** TASK/TR  3.00  2011-11-20 *****'

  CALL GSOPEN

  CALL eq_init
!  CALL equ_init
  CALL tr_init

  CALL tr_parm(1,'trparm',ierr)
  CALL eq_parm(1,'eqparm',ierr)
!  CALL equ_parm(1,'equparm',ierr)
  
  CALL tr_menu

  CALL tr_term

  CALL GSCLOS
  STOP
END PROGRAM tr
