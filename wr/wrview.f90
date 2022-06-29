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
    WRITE(6,602) 'SMAX  ',SMAX  ,'DELS  ',DELS  , &
                 'UUMIN ',UUMIN
    WRITE(6,602) 'EPSRAY',EPSRAY,'DELRAY',DELRAY, &
                 'DELDER',DELDER
    WRITE(6,602) 'DELKR ',DELKR, 'EPSNW ',EPSNW
    WRITE(6,602) 'DELKR ',DELKR, 'EPSNW ',EPSNW
    WRITE(6,601) 'LMAXNW',LMAXNW
    WRITE(6,601) 'MDLWRI',MDLWRI,'MDLWRG',MDLWRG, &
                 'MDLWRP',MDLWRP
    WRITE(6,601) 'MDLWRQ',MDLWRQ,'MDLWRW',MDLWRW
    WRITE(6,603) 'nres_max    ',nres_max, &
                 'nres_type   ',nres_type,&
                 'mode_beam   ',mode_beam
    WRITE(6,603) 'NPMAX_DP    ',NPMAX_DP, &
                 'NTHMAX_DP   ',NTHMAX_DP,&
                 'NRMAX_DP    ',NRMAX_DP
    WRITE(6,604) 'pne_threshold   ',pne_threshold
    WRITE(6,603) 'Rmax_wr ',Rmax_wr,'Rmin_wr ',Rmin_wr, &
                 'Zmax_wr ',Zmax_wr,'Zmin_wr ',Zmin_wr
    RETURN

601 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
602 FORMAT(1H ,A6,'=',ES12.4:2X,A6,'=',ES12.4: &
            2X,A6,'=',ES12.4:2X,A6,'=',ES12.4)
603 FORMAT(1H ,A12,'=',I7,4X  :2X,A12,'=',I7,4X  : &
            2X,A12,'=',I7)
604 FORMAT(1H ,A16'=',ES12.4:2X,A16,'=',ES12.4)
  END SUBROUTINE WR_VIEW
END MODULE wrview
