

!     ***** TASK/TR MENU *****

      SUBROUTINE TRMENU
      use T3D_FILE_IO
   
      
      use trunit
      use t3d_er_param
    
      IMPLICIT NONE
      INTEGER(4)       :: IERR
      EXTERNAL TRPARM
      integer(4):: i

!     ------ SELECTION OF TASK TYPE ------

      IERR=0
      
if(flg_SHOT_LOOP_CONTROL == 1) then

    call PARAMETER_SURVEY_LOOP

else
    flg_ERF_Difference_Coptimize=0
    flg_t3d_er_calc_timing =0
    CNT_TREXE_ITERATION =-1

    call T3D_FILENAME_MAKER

    CALL tr_prof(ierr)
    CALL TRLOOP
endif

    
    RETURN
    END SUBROUTINE TRMENU
