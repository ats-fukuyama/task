MODULE trgout
! --------------------------------------------------------------------------
!   The variables declared below should be refered 
!    from ONLY 'trn/trgraphics' directory
! ---------------------------------------------------------------------------

  USE trcomm, ONLY : ikind,rkind

  PRIVATE
  PUBLIC tr_gout,tr_gr_time

CONTAINS

  SUBROUTINE tr_gout
! ***********************************************************
!     TASK/TR graphic outputs control routine
! ***********************************************************

!    USE trcomm, ONLY : NGR,NGT
    USE trgrad, ONLY: tr_gr_radial
    USE trgtmp, ONLY: tr_gr_temporal
    USE trgcom, ONLY: tr_gr_comp
    USE trgdgn, ONLY: tr_gr_diagnostic
    USE trgexp, ONLY: tr_gr_exp
    IMPLICIT NONE

    INTEGER(ikind), SAVE :: init = 0, inqg
    CHARACTER(LEN=5) :: kig
    CHARACTER(LEN=3) :: kk
    CHARACTER(LEN=1) :: k1, k2, k3, k4, k5
    INTEGER(ikind)   :: iosts

    INTEGER(ikind) :: NGR,NGT

    IF(init.EQ.0) THEN
       inqg=0
       init=1
    ENDIF

    DO
       WRITE(6,*) '# Graph select : R1-13, T1-2, N1-12, D1, U1-2'
       WRITE(6,*) '#  Menu select : S/save  L/load  H/help  ',&
                                   'C/clear  I/inq  X/exit'
       READ(5,'(A5)',iostat=iosts) KIG
       if(iosts > 0) then
          cycle
       elseif(iosts < 0) then
          exit
       end if
       k1=kig(1:1)
       k2=kig(2:2)
       k3=kig(3:3)
       k4=kig(4:4)
       k5=kig(5:5)
       ! Capital letter to small letter (in GSAF lib)
       CALL GUCPTL(k1)
       CALL GUCPTL(k2)
       CALL GUCPTL(k3)
       CALL GUCPTL(k4)
       CALL GUCPTL(k5)
       kk=k3//k4//k5

       SELECT CASE(k1)
!!$         CASE('H') ! help
!!$            CALL TRHELP('G')
!!$            CYCLE
       CASE('C') ! clear
          NGR=0
          NGT=0
       CASE('I') ! inquire
          IF(inqg.EQ.0) THEN
             inqg=4
             WRITE(6,*) '## Graphic scale inquire mode : ON'
          ELSE
             inqg=0
             WRITE(6,*) '## Graphic scale inquire mode : OFF'
          ENDIF
       CASE('S') ! save
          CALL tr_gr_save
          CYCLE
       CASE('L') ! load
          CALL tr_gr_load
          CYCLE
       CASE('R') ! snap shot of radial profile
          CALL tr_gr_radial(k2,k3)
          CYCLE
       CASE('T')
          CALL tr_gr_temporal(k2)
          CYCLE
       CASE('N')
          CALL tr_gr_comp(k2,k3)
          CYCLE
       CASE('D')
          CALL tr_gr_diagnostic(k2)
          CYCLE
       CASE('U')
          CALL tr_gr_exp(k2,k3)
          CYCLE
!!$         CASE('Y')
!!$            CALL TRGRY0(k2,inqg)
!!$         CASE('Z')
!!$            CALL TRGRX0(k2,inqg)
!!$         CASE('G')
!!$            CALL TRGRG0(k2,inqg)
!!$         CASE('A')
!!$            CALL TRGRA0(k2,inqg)
!!$         CASE('E')
!!$            CALL TRGRE0(k2,inqg)
!!$         CASE('M')
!!$            CALL TRCOMP(k2,inqg)
       CASE('X')
          EXIT
       CASE default
          WRITE(6,*) ' ERROR : Unsupported graph ID'
          CYCLE
       END SELECT
       EXIT
    END DO

    RETURN
  END SUBROUTINE tr_gout


  SUBROUTINE tr_gr_save
! ***********************************************************
!           SAVE GRAPHIC DATA
! ***********************************************************
!    USE trcomm, ONLY : GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT
    IMPLICIT NONE

    INTEGER(ikind)    :: iosts
    CHARACTER(LEN=32) :: tr_gfile_name
    CHARACTER(LEN=1)  :: kid
    LOGICAL           :: lex

!    INTEGER(ikind) :: GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT

    DO
       WRITE(6,*) '# INPUT : Graphic save file name (CR to CANCEL)'
       READ(5,'(A32)',IOSTAT=iosts) tr_gfile_name
       IF(iosts /= 0) CYCLE

       IF(tr_gfile_name.EQ.'                                ') RETURN
       INQUIRE(FILE=tr_gfile_name, EXIST=lex)

       IF(lex) THEN
          DO
             WRITE(6,*) '# Old file is going to be OVERWRITTEN. ', &
                        'Are you sure {y/n}?'
             READ(5,'(A1)',IOSTAT=iosts) kid
             IF(iosts /= 0)THEN
                WRITE(6,*) ' XX Invalid Input.'
                CYCLE
             END IF
             EXIT
          END DO
          
          CALL GUCPTL(kid)
          IF(kid.NE.'Y') CYCLE

          OPEN(22,FILE=tr_gfile_name,IOSTAT=iosts,STATUS='OLD', &
                                                 FORM='UNFORMATTED')
          IF(iosts /= 0) THEN
             WRITE(6,*) '# XX Old file open ERROR : IOSTAT=',iosts
             CYCLE
          END IF

          WRITE(6,*) '# Old file (',tr_gfile_name,') is assigned for output.'

       ELSE
          OPEN(22,FILE=tr_gfile_name,IOSTAT=iosts,STATUS='NEW', &
                                                 FORM='UNFORMATTED')
          IF(iosts /= 0)THEN
             WRITE(6,*) 'XX New file open ERROR : IOSTAT=',iosts
             CYCLE
          END IF

          WRITE(6,*) '# New file (',tr_gfile_name,') is created for output.'
       ENDIF
       EXIT
    END DO

!    WRITE(22) GVR,GRM,GRG,GTR,NGR
!    WRITE(22) GVT,GT,NGT
!    CLOSE(22)

    WRITE(6,*) '# Data was SUCCESSFULLY SAVED to the file.'

    RETURN
  END SUBROUTINE tr_gr_save


  SUBROUTINE tr_gr_load
! ***********************************************************
!           LOAD GRAPHIC DATA
! ***********************************************************
!    USE trcomm, ONLY : GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT
    IMPLICIT NONE

    INTEGER(ikind)    :: iosts
    CHARACTER(LEN=32) :: tr_gfile_name
    LOGICAL           :: lex

    INTEGER(ikind) :: GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT

    DO
       WRITE(6,*) '# INPUT : Graphic load file name (CR to CANCEL)'
       READ(5,'(A32)',IOSTAT=iosts) tr_gfile_name
       IF(iosts /= 0) CYCLE

       IF(tr_gfile_name.EQ.'                                ') RETURN
       INQUIRE(FILE=tr_gfile_name,EXIST=lex)

       IF(lex) THEN
          OPEN(22,FILE=tr_gfile_name,IOSTAT=iosts,STATUS='OLD', &
                                                 FORM='UNFORMATTED')
          IF(iosts /= 0)THEN
             WRITE(6,*) '# XX OLD FILE OPEN ERROR : IOSTAT=',IOSTS
             CYCLE
          END IF

          WRITE(6,*) '# File (',tr_gfile_name,') is assigned for input.'

       ELSE
          WRITE(6,*) 'XX File (',tr_gfile_name,') not found.'
          CYCLE
       ENDIF
       EXIT
    END DO

    READ(22) GVR,GRM,GRG,GTR,NGR
    READ(22) GVT,GT,NGT
    CLOSE(22)

    WRITE(6,*) '# Data was SUCCESSFULLY LOADED from the file.'

    RETURN
  END SUBROUTINE tr_gr_load


  SUBROUTINE tr_gr_time
! ***********************************************************************
!            Write time on figure
! ***********************************************************************

    USE trcomm, ONLY : t
    IMPLICIT NONE

    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.0)
    CALL SETFNT(32)
    CALL MOVE(11.8,18.0)
    CALL TEXT('t =',2)
    CALL NUMBD(t,'(1F7.3)',7)
    CALL TEXT(' sec.',4)
    RETURN
  END SUBROUTINE tr_gr_time

END MODULE trgout
