! wrview.f90

MODULE wrview

  PRIVATE
  PUBLIC wr_view

CONTAINS

!     ****** SHOW PARAMETERS ******

  SUBROUTINE WR_VIEW

    USE wrcomm_parm
    IMPLICIT NONE

    WRITE(6,603) 'NRAYMAX     ',NRAYMAX, &
                 'NSTPMAX     ',NSTPMAX
    WRITE(6,603) 'NRSMAX      ',NRSMAX, &
                 'NRLMAX      ',NRLMAX
    WRITE(6,601) 'SMAX  ',SMAX  ,'DELS  ',DELS  , &
                 'UUMIN ',UUMIN
    WRITE(6,601) 'EPSRAY',EPSRAY,'DELRAY',DELRAY, &
                 'DELDER',DELDER
    WRITE(6,601) 'DELKR ',DELKR, 'EPSNW ',EPSNW
    WRITE(6,602) 'LMAXNW',LMAXNW
    WRITE(6,602) 'MDLWRI',MDLWRI,'MDLWRG',MDLWRG, &
                 'MDLWRP',MDLWRP
    WRITE(6,602) 'MDLWRQ',MDLWRQ,'MDLWRW',MDLWRW
    WRITE(6,603) 'nres_max    ',nres_max, &
                 'nres_type   ',nres_type,&
                 'mode_beam   ',mode_beam
    WRITE(6,603) 'NPMAX_DP    ',NPMAX_DP, &
                 'NTHMAX_DP   ',NTHMAX_DP,&
                 'NRMAX_DP    ',NRMAX_DP
    RETURN

601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
603 FORMAT(1H ,A12,'=',I7,4X  :2X,A12,'=',I7,4X  : &
            2X,A12,'=',I7)
  END SUBROUTINE WR_VIEW
END MODULE wrview
