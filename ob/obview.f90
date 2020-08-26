! obview.f90

MODULE obview

CONTAINS

  !     ****** SHOW PARAMETERS ******

  SUBROUTINE ob_view

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER:: nobt

    WRITE(6,603) &
         'nobt_max    ',nobt_max, &
         'nstp_max    ',nstp_max, &
         'ns_ob       ',ns_ob
    WRITE(6,603) &
         'lmax_nw     ',lmax_nw, &
         'mdlobp      ',mdlobp, &
         'mdlobi      ',mdlobi
    WRITE(6,603) &
         'mdlobq      ',mdlobq, &
         'mdlobt      ',mdlobt, &
         'mdlobc      ',mdlobc
    WRITE(6,603) &
         'mdlobw      ',mdlobw, &
         'mdlobg      ',mdlobg
    WRITE(6,604) &
         'tmax        ',tmax, &
         'delt        ',delt
    WRITE(6,604) &
         'eps_obt     ',eps_obt, &
         'del_obt     ',del_obt, &
         'eps_nw      ',eps_nw
    SELECT CASE(mdlobi)
    CASE(0)
       WRITE(6,'(A)') &
            'nobt  penergy_in  pcangle_in  zeta_in     psipn_in    theta_in'
       WRITE(6,'(I4,1P5E12.4)') &
            (nobt,penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
             psipn_in(nobt),theta_in(nobt),nobt=1,nobt_max)
    CASE(1)
       WRITE(6,'(A)') &
            'nobt  penergy_in  pcangle_in  zeta_in     rr_in       zz_in'
       WRITE(6,'(I4,1P5E12.4)') &
            (nobt,penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
             rr_in(nobt),zz_in(nobt),nobt=1,nobt_max)
    CASE DEFAULT
       WRITE(6,*) 'XX obview: undefined mdlobi: mdlobi=',mdlobi
    END SELECT
    WRITE(6,603) &
         'nrmax_ob    ',nrmax_ob, &
         'nthmax_ob   ',nthmax_ob, &
         'nsumax_ob   ',nsumax_ob
    RETURN

601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
603 FORMAT(1H ,A12,'=',I7,4X  :2X,A12,'=',I7,4X  : &
            2X,A12,'=',I7)
604 FORMAT(1H ,A12,'=',1PE11.3:2X,A12,'=',1PE11.3: &
            2X,A12,'=',1PE11.3)
  END SUBROUTINE ob_view
END MODULE obview
