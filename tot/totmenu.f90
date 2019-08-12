! totmenu.f90

MODULE totmenu

  PRIVATE
  PUBLIC tot_menu

CONTAINS

!     ***** TOT MENU *****

  SUBROUTINE tot_menu

    USE plmenu,ONLY:pl_menu
    USE dpmenu,ONLY:dp_menu
    USE fpmenu,ONLY:fp_menu
    USE timenu,ONLY:ti_menu
    USE wrmenu,ONLY:wr_menu
    USE commpi
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=2) :: KID
    CHARACTER(LEN=1) :: KID1,KID2

1   CONTINUE
    IF(nrank.EQ.0) THEN
       WRITE(6,601)
601    FORMAT('## TASK MENU: EQ/TR/WR/WM/FP/DP/PL/TI  Q/QUIT')
       READ(5,'(A2)') KID
       KID1=KID(1:1)
       KID2=KID(2:2)
       CALL GUCPTL(KID1)
       CALL GUCPTL(KID2)
    ENDIF
    CALL mtx_broadcast1_character(KID1)
    CALL mtx_broadcast1_character(KID2)
    KID=KID1//KID2

    IF (KID1.EQ.'Q') GOTO 9000

    SELECT CASE(KID)
    CASE('EQ')
       CALL eqmenu
    CASE('TR')
       CALL trmenu
    CASE('WR')
       CALL wr_menu
    CASE('WM')
       CALL wmmenu
    CASE('FP')
       CALL fp_menu
    CASE('DP')
       CALL dp_menu
    CASE('PL')
       CALL pl_menu
    CASE('TI')
       CALL ti_menu
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE tot_menu
END MODULE totmenu
