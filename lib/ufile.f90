module ufile_integer
  integer(4) :: MODEL
end module ufile_integer

!***********************************************************
!
!      LOCATING UFILE DIRECTORY
!
!***********************************************************

  SUBROUTINE UFILE_INTERFACE(KDIRX,KUFDIR,KUFDEV,KUFDCG,MODE)

    use ufile_integer
    implicit none
    integer(4),        intent(in)  :: MODE
    character(len=80), intent(in)  :: KUFDIR, KUFDEV, KUFDCG
    character(len=80), intent(out) :: KDIRX
    integer(4)         :: IKNDEV, IKNDCG, IKDIRX, IST
    character(len=100) :: KFILE
!    COMMON /UFMODE/ MODEL
    LOGICAL :: DIR, FILE

!     MODE is the parameter which determines how to handle exp. files.
!     MODE = 0 : Binary files are loaded if available, or ASCII files
!                are loaded and aftermath binary files are created.
!            1 : Only binary files are loaded.
!         else : Only ASCII files are loaded and binary files are NOT
!                created.

    MODEL = MODE

    IKNDEV = len_trim(KUFDEV)
    IKNDCG = len_trim(KUFDCG)

    KDIRX=TRIM(KUFDIR)//'/'//KUFDEV(1:IKNDEV)//'/' &
         &                          //KUFDCG(1:IKNDCG)//'/in/'

    IKDIRX = len_trim(KDIRX)
    KFILE=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'2d'//KUFDCG(1:IKNDCG)//'.NE'
    INQUIRE(FILE=KDIRX,EXIST=DIR,IOSTAT=IST)
    IF(IST > 0) RETURN
    IF(DIR.EQV..FALSE.) THEN

!     ifort compiler cannot recognize whether the directory exists or
!     not (always says "FALSE"), so that NE ufile must be checked in
!     case of using ifort because it is the most possible existence file
!     in a set of ufiles.

       INQUIRE(FILE=KFILE,EXIST=FILE,IOSTAT=IST)
       IF(IST > 0) RETURN
       IF(FILE.EQV..FALSE.) THEN
          WRITE(6,'(A)') '## DESIGNATED DIRECTORY DOES NOT EXIST!'
          WRITE(6,'(A7,A80)') '   ==> ',KDIRX
          STOP
       ENDIF
    ENDIF

  END SUBROUTINE UFILE_INTERFACE

!**************************************************
!**    UFILE read for TR (Time Evolution UFILE)  **
!**************************************************

!  input:

!  KDIRX    : Directory
!  KUFDEV   : Device name in UFILEs
!  KUFDCG   : Shot number in UFILEs
!  KFID     : Variable name in UFILEs

!  output:

!  RUF(NRMU)      : Equally Spaced Normalized Radial Data
!  TMU(NTUM)      : Total Time Data (The Number of DT)
!  F1(NTUM)       : Functional Values
!  F2(NTUM,NRMU)  : Functional Values
!  NRFMAX         : Maximum Number of the Radial Mesh
!  NTXMAX         : Maximum Number of the Time Mesh
!  NRMU           : Array size (radial)
!  NTUM           : Array size (time)
!  MDCHK          : Loop Check Value
!  IERR           : Error Indicator

!   ***************************************************************

  SUBROUTINE UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TT,F1,NTXMAX,NTUM,MDCHK,IERR)

    use ufile_integer
    implicit none
!    COMMON /UFMODE/ MODEL
    character(len=80),        intent(in)  :: KDIRX, KUFDEV, KUFDCG
    character(len=10),        intent(in)  :: KFID
    integer(4),               intent(in)  :: NTUM
    real(8), dimension(NTUM), intent(out) :: TT
    real(8), dimension(NTUM), intent(out) :: F1
    integer(4),               intent(out) :: NTXMAX
    integer(4),               intent(out) :: MDCHK, IERR
    integer(4) :: IKNDEV, IKNDCG, IKDIRX, KL1, IST
    character(len=100) :: KDIRR1, KFILE
    LOGICAL :: LEX

    IKNDEV = len_trim(KUFDEV)
    IKNDCG = len_trim(KUFDCG)

!    IF(MDCHK /= 0) RETURN

    IKDIRX = len_trim(KDIRX)
    KDIRR1=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'1d'//KUFDCG(1:IKNDCG)//'.'

    KL1 = len_trim(KDIRR1)
    KFILE=KDIRR1(1:KL1)//KFID

    INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=IST)
    IF(IST > 0) RETURN
    IF(LEX) THEN
       CALL TRXR1D(KDIRR1,KFID,TT,F1,NTUM,NTXMAX,MODEL)
       MDCHK=1
       IERR=0
    ELSE
       F1(1:NTUM)=0.D0
       IERR=1
    ENDIF

  END SUBROUTINE UFREAD_TIME

!*****

  SUBROUTINE UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,F2, &
       &                        NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)

    use ufile_integer
    implicit none
!    COMMON /UFMODE/ MODEL
    integer(4), intent(in)  :: NRMU, NTUM
    character(len=80), intent(in) :: KDIRX, KUFDEV, KUFDCG
    character(len=10), intent(in) :: KFID
    integer(4),                    intent(out) :: MDCHK, IERR
    real(8), dimension(NRMU),      intent(out) :: RUF
    real(8), dimension(NTUM),      intent(out) :: TMU
    real(8), dimension(NTUM,NRMU), intent(out) :: F2
    integer(4) :: IKNDEV, IKNDCG, IKDIRX, KL2, NRFMAX, NTXMAX, NTX, NRF, MD, IST
    real(8) :: FCTR, AITKEN2P
    real(8), dimension(NTUM) :: F2CTR, F2EDG
    real(8), dimension(NRMU,NTUM) :: F2I
    character(len=100) :: KDIRR2, KFILE
    LOGICAL :: LEX

    IKNDEV = len_trim(KUFDEV)
    IKNDCG = len_trim(KUFDCG)

!    IF(MDCHK /= 0) RETURN

    IKDIRX = len_trim(KDIRX)
    KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'2d'//KUFDCG(1:IKNDCG)//'.'

    KL2 = len_trim(KDIRR2)
    KFILE=KDIRR2(1:KL2)//KFID

    INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=IST)
    IF(IST > 0) RETURN
    IF(LEX) THEN
       NRFMAX = NRMU
       CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2I,NRMU,NTUM,NRFMAX,NTXMAX,MODEL)
       IF(NRFMAX == 1) THEN
          WRITE(6,*) '## ',KFID,'HAS NOT BEEN READ DUE TO A SINGLE RADIAL POINT.'
          F2(1:NTUM,1:NRMU)=0.D0
          IERR=2
          RETURN
       ENDIF
       MDCHK=1
       IERR=0
       F2(1:NTXMAX,1:NRFMAX) = transpose(F2I(1:NRFMAX,1:NTXMAX))
    ELSE
       F2(1:NTUM,1:NRMU)=0.D0
       IERR=1
       RETURN
    ENDIF

!*****
!  For interpolation, center and/or edge value is added
!  if necessary. The variables with center value having to
!  be zero such as BPOL, RMINOR, SURF and VOLUME were beforehand 
!  manupulated and thus in this section MD=4 will be selected.
!*****

    MD=0
    IF(RUF(1) /= 0.D0.AND.RUF(NRFMAX) /= 1.D0) THEN
       DO NTX=1,NTXMAX
          F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
          F2EDG(NTX)=AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1), &
               &                   F2(NTX,NRFMAX-2),RUF(NRFMAX), &
               &                   RUF(NRFMAX-1),RUF(NRFMAX-2))
       ENDDO
       MD=1
    ELSEIF(RUF(1) /= 0.D0.AND.RUF(NRFMAX) == 1.D0) THEN
       DO NTX=1,NTXMAX
          F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
       ENDDO
       MD=2
    ELSEIF(RUF(1) == 0.D0.AND.RUF(NRFMAX) /= 1.D0) THEN
       DO NTX=1,NTXMAX
          F2EDG(NTX)=AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1), &
     &                             F2(NTX,NRFMAX-2),RUF(NRFMAX), &
     &                             RUF(NRFMAX-1),RUF(NRFMAX-2))
       ENDDO
       MD=3
    ELSEIF(RUF(1) == 0.D0.AND.RUF(NRFMAX) == 1.D0) THEN
       MD=4
    ELSE
       STOP 'XX UFREAD2_TIME: TRXR2D: ERROR'
    ENDIF

    CALL DATA_ERROR_CORRECT(KUFDEV,KUFDCG,KFID,RUF,F2,NTXMAX,NRMU,NTUM)

    select case(MD)
    case(1)
       NRFMAX=NRFMAX+2
       RUF(2:NRFMAX-1) = RUF(1:NRFMAX-2)
       RUF(1)      = 0.D0
       RUF(NRFMAX) = 1.D0
       DO NTX=1,NTXMAX
          F2(NTX,2:NRFMAX-1) = F2(NTX,1:NRFMAX-2)
          F2(NTX,1)=F2CTR(NTX)
          F2(NTX,NRFMAX)=F2EDG(NTX)
       ENDDO
    case(2)
       NRFMAX=NRFMAX+1
       RUF(2:NRFMAX) = RUF(1:NRFMAX-1)
       RUF(1)=0.D0
       DO NTX=1,NTXMAX
          F2(NTX,2:NRFMAX) = F2(NTX,1:NRFMAX-1)
          F2(NTX,1)=F2CTR(NTX)
       ENDDO
    case(3)
       NRFMAX=NRFMAX+1
       RUF(NRFMAX)=1.D0
       F2(1:NTXMAX,NRFMAX)=F2EDG(1:NTXMAX)
    end select

  END SUBROUTINE UFREAD2_TIME

!***********************************************************
!
!     WRONG DATA CORRECT
!
!***********************************************************

  SUBROUTINE DATA_ERROR_CORRECT(KUFDEV,KUFDCG,KFID,RUF,F2,NTXMAX,NRMU,NTUM)

    implicit none
    integer(4) :: NTXMAX, NRMU, NTUM
    character(len=80), intent(in) :: KUFDEV, KUFDCG
    character(len=10), intent(in) :: KFID
    real(8), dimension(NRMU) :: RUF
    real(8), dimension(NTUM,NRMU) :: F2
    integer(4) :: NTX
    real(8) :: FCTR

    IF(KUFDEV == 'jet' .AND. KFID == 'GRHO1') THEN
       IF(KUFDCG == '57987' .OR. KUFDCG == '58159' .OR. KUFDCG == '58323') THEN
          DO NTX=1,NTXMAX
             F2(NTX,1)=FCTR(RUF(2),RUF(3),F2(NTX,2),F2(NTX,3))
          ENDDO
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE DATA_ERROR_CORRECT

!*****

  SUBROUTINE UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,F2, &
       &                         NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)

    use ufile_integer
    implicit none
!    COMMON /UFMODE/ MODEL
    integer(4),        intent(in) :: NRMU, NTUM
    character(len=80), intent(in) :: KDIRX, KUFDEV, KUFDCG
    character(len=10), intent(in) :: KFID
    integer(4),                    intent(out) :: MDCHK, IERR
    real(8), dimension(NRMU),      intent(out) :: RUF
    real(8), dimension(NTUM),      intent(out) :: TMU
    real(8), dimension(NTUM,NRMU), intent(out) :: F2
    integer(4) :: IKNDEV, IKNDCG, IKDIRX, KL2, NRFMAX, NTXMAX, NRF, NTX, IST
    real(8), dimension(NRMU,NTUM) :: F2I
    character(len=100) :: KDIRR2, KFILE
    LOGICAL :: LEX

    IKNDEV = len_trim(KUFDEV)
    IKNDCG = len_trim(KUFDCG)

!    IF(MDCHK /= 0) RETURN

    IKDIRX = len_trim(KDIRX)
    KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'2d'//KUFDCG(1:IKNDCG)//'.'

    KL2 = len_trim(KDIRR2)
    KFILE=KDIRR2(1:KL2)//KFID

    INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=IST)
    IF(IST > 0) RETURN
    IF(LEX) THEN
       NRFMAX=NRMU           ! equal to NRMU
       CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2I,NRMU,NTUM,NRFMAX,NTXMAX,MODEL)
       MDCHK=1
       IERR=0
       F2(1:NTXMAX,1:NRFMAX) = transpose(F2I(1:NRFMAX,1:NTXMAX))
    ELSE
       F2(1:NTUM,1:NRMU)=0.D0
       IERR=1
       RETURN
    ENDIF

    IF(RUF(NRFMAX) > 1.D0) THEN
       DO NRF=1,NRFMAX
          IF(RUF(NRF) > 1.D0) THEN
             NRFMAX=NRF-1
             EXIT
          ENDIF
       ENDDO
    ENDIF
    IF(KFID == 'NEXP'.OR.KFID == 'NEXPEB') THEN
       F2(1:NTXMAX,1:NRMU)=F2(1:NTXMAX,1:NRMU)*1.D-20
    ELSE
       F2(1:NTXMAX,1:NRMU)=F2(1:NTXMAX,1:NRMU)*1.D-3
    ENDIF

  END SUBROUTINE UFREAD2_ERROR

!***** INITIALIZATION FOR READING UFILE *****

  SUBROUTINE NDINIT(NDPOS)

    implicit none
    integer(4), intent(out) :: NDPOS
    integer(4), parameter :: NLENM=80

    NDPOS=NLENM
    RETURN
  END SUBROUTINE NDINIT

!***** READING UFILE *****

  SUBROUTINE NDREAD(NDPOS,IRD,XD,IERR)

    implicit none
    integer(4), intent(inout) :: NDPOS
    integer(4), intent(in)    :: IRD
    integer(4), intent(out)   :: IERR
    real(8),    intent(out)   :: XD
    integer(4), parameter :: NLENM=80
    integer(4) :: NR, NT, NDPOS1, NDPOS2, IST
    character(len=NLENM) :: LINE
    character(len=40)    :: NUMD

    OUT:DO 
       IF(NDPOS == NLENM) THEN
          READ(IRD,'(A80)',IOSTAT=IST) LINE
          IF(IST > 0) THEN
             IERR = 1
             RETURN
          ELSEIF(IST < 0) THEN
             IERR = 2
             RETURN
          END IF
          NDPOS=0
       ENDIF

       IN1:DO 
          NDPOS=NDPOS+1
          IF(LINE(NDPOS:NDPOS) == ' '.AND.NDPOS < NLENM) THEN
             CYCLE IN1
          ELSE
             EXIT IN1
          END IF
       END DO IN1
    
       IF(NDPOS < NLENM) THEN
          NDPOS1=NDPOS
          IN2:DO 
             NDPOS=NDPOS+1
             IF(LINE(NDPOS:NDPOS) /= 'E'.AND.LINE(NDPOS:NDPOS) /= 'e'.AND. &
                  &      NDPOS < NLENM) THEN
                CYCLE IN2
             ELSE
                EXIT IN2
             END IF
          END DO IN2

          IF(NDPOS < NLENM) THEN
             NDPOS=NDPOS+1
             IN3:DO 
                NDPOS=NDPOS+1
                IF(LINE(NDPOS:NDPOS) /= ' '.AND.LINE(NDPOS:NDPOS) /= '-'.AND. &
                     &         NDPOS < NLENM) THEN
                   CYCLE IN3
                ELSE
                   EXIT IN3
                END IF
             END DO IN3

             IF(NDPOS /= NLENM) NDPOS=NDPOS-1
             NDPOS2=NDPOS
             NUMD=LINE(NDPOS1:NDPOS2)
             READ(NUMD,*,IOSTAT=IST) XD
             IF(IST > 0) THEN
                NDPOS=NLENM
                CYCLE OUT
             END IF
             IERR=0
             RETURN
          ENDIF
       ENDIF
    END DO OUT

  END SUBROUTINE NDREAD

!***** READ 1D VARIABLE *****

  SUBROUTINE TRXR1D(KDIRR1,KFID,T,F1,NTM,NTXMAX,MODEL)

    implicit none
    integer(4),              intent(in)  :: NTM, MODEL
    character(len=80),       intent(in)  :: KDIRR1
    character(len=10),       intent(in)  :: KFID
    integer(4),              intent(out) :: NTXMAX
    real(8), dimension(NTM), intent(out) :: T, F1
    integer(4) :: KL1, KL2, MODE, IST, IMATCH, NTX, KL, INCL, NDPOS, IERR, NT
    character(len=80) :: KFILE, KFILEB, KMATCH1, KMATCH2, KMATCH3, KKLINE, KKLINE1, KKLINE2
    character(len=15) :: KFIDB
    LOGICAL :: KMATCH, LEX

    MODE=MODEL
    IMATCH=0
    KMATCH1='# OF X PTS'
    KMATCH2='# OF Y PTS'
    KMATCH3='# OF PTS'

    KL1 = len_trim(KDIRR1)
    KFILE=KDIRR1(1:KL1)//KFID
    KL2 = len_trim(KFILE)
    KFILEB=KFILE(1:KL2)//'.bin'

    INQUIRE(FILE=KFILEB,EXIST=LEX,IOSTAT=IST)
    IF(IST > 0) RETURN
    IF(LEX.AND.MODE == 0) MODE=1

    IF(MODE == 1) THEN
       OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='OLD',ERR=8)
       READ(15) NTXMAX
       READ(15) (T(NTX),NTX=1,NTXMAX)
       READ(15) (F1(NTX),NTX=1,NTXMAX)
       CLOSE(15)
       KFIDB=KFID//'.bin'
       WRITE(6,"(' ',3A,I3,A)") ' - READ FILE :',KFIDB,'(',NTXMAX,')'
       RETURN
    ENDIF

8   OPEN(15,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',STATUS='OLD')
    IF(IST > 0) THEN
       WRITE(6,*) 'XX NOT FOUND :',KFILE(1:45)
       NTXMAX=0
       WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
       RETURN
    ENDIF

    DO 
       READ(15,'(A80)',IOSTAT=IST) KKLINE
       IF(IST > 0) THEN
          WRITE(6,*) 'XX FILE READ ERROR'
          RETURN
       ELSEIF(IST < 0) THEN
          RETURN
       END IF
       KL = len_trim(KKLINE)
       IF(KL == 0) CYCLE
       INCL=INDEX(KKLINE,';')
       IF(INCL == 0) THEN
          IF(IMATCH >= 1) EXIT
          CYCLE
       ENDIF
       CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
       KL = len_trim(KKLINE1)
       CALL KEXTR(KKLINE2)
       IF(KMATCH(KKLINE2,KMATCH1)) THEN
          READ(KKLINE1,*,ERR=8000) NTXMAX
          IMATCH=IMATCH+1
       ELSE IF(KMATCH(KKLINE2,KMATCH2)) THEN
          READ(KKLINE1,*,ERR=8000) NTXMAX
          IMATCH=IMATCH+1
       ELSE IF(KMATCH(KKLINE2,KMATCH3)) THEN
          READ(KKLINE1,*,ERR=8000) NTXMAX
          IMATCH=IMATCH+1
       ENDIF
    END DO

    BACKSPACE 15

    CALL NDINIT(NDPOS)
    DO NTX=1,NTXMAX
       CALL NDREAD(NDPOS,15,T(NTX),IERR)
       IF(IERR /= 0) GOTO 8000
    ENDDO
    DO NTX=1,NTXMAX
       CALL NDREAD(NDPOS,15,F1(NTX),IERR)
       IF(IERR /= 0) GOTO 8000
    ENDDO

    CLOSE(15)
    WRITE(6,"(' ',3A,I3,A)") ' - READ ASCII FILE :',KFID,'(',NTXMAX,')'

    IF(MODE == 0) THEN
       OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='NEW')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
          RETURN
       ENDIF
       WRITE(15) NTXMAX
       WRITE(15) (T(NT),NT=1,NTXMAX)
       WRITE(15) (F1(NT),NT=1,NTXMAX)
       CLOSE(15)
       WRITE(6,*) ' - WRITE FILE:',KFILEB(1:60)
    ENDIF
    RETURN

8000 WRITE(6,*) 'XX FILE READ ERROR'
    RETURN
  END SUBROUTINE TRXR1D

!***** READ 2D VARIABLE *****

  SUBROUTINE TRXR2D(KDIRR2,KFID,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,MODEL)

    implicit none
    integer(4),                  intent(in)  :: NRM, NTM, MODEL
    character(len=80),           intent(in)  :: KDIRR2
    character(len=10),           intent(in)  :: KFID
    integer(4),                  intent(out) :: NRXMAX, NTXMAX
    real(8), dimension(NTM),     intent(out) :: T
    real(8), dimension(NRM),     intent(out) :: R
    real(8), dimension(NRM,NTM), intent(out) :: F2
    integer(4) :: KL2, MODE, IST, IMATCH, NRX, NTX, KL, INCL, NDPOS, IERR, NT
    character(len=80) :: KFILE, KFILEB, KMATCH1, KMATCH2, KMATCH3, KKLINE, KKLINE1, KKLINE2
    character(len=15) :: KFIDB
    LOGICAL :: KMATCH, LEX

    MODE=MODEL
    IMATCH=0
    KMATCH1='# OF X PTS'
    KMATCH2='# OF Y PTS'
    KMATCH3='# OF PTS'

    KL2 = len_trim(KDIRR2)
    KFILE=KDIRR2(1:KL2)//KFID
    KL2 = len_trim(KFILE)
    KFILEB=KFILE(1:KL2)//'.bin'

    INQUIRE(FILE=KFILEB,EXIST=LEX,IOSTAT=IST)
    IF(IST > 0) GOTO 9000
    IF(LEX.AND.MODE == 0) MODE=1

    IF(MODE == 1) THEN
       OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='OLD',ERR=8)
       READ(15) NRXMAX,NTXMAX
       READ(15) (R(NRX),NRX=1,NRXMAX)
       READ(15) (T(NTX),NTX=1,NTXMAX)
       READ(15) ((F2(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)
       CLOSE(15)
       KFIDB=KFID//'.bin'
       WRITE(6,"(' ',3A,2(I3,A))") ' - READ FILE :',KFIDB,'(',NRXMAX,',',NTXMAX,')'
       GOTO 9000
    ENDIF

8   OPEN(15,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',STATUS='OLD')
    IF(IST > 0) THEN
       WRITE(6,*) 'XX NOT FOUND :',KFILE(1:55)
       NTXMAX=0
       WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
       RETURN
    ENDIF

    DO
       READ(15,'(A80)',IOSTAT=IST) KKLINE
       IF(IST > 0) THEN
          WRITE(6,*) 'XX FILE READ ERROR'
          RETURN
       ELSEIF(IST < 0) THEN
          GOTO 9000
       END IF
       KL = len_trim(KKLINE)
       IF(KL == 0) CYCLE
       INCL=INDEX(KKLINE,';')
       IF(INCL == 0) THEN
          IF(IMATCH >= 2) EXIT
          CYCLE
       ENDIF
       CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
       KL = len_trim(KKLINE1)
       CALL KEXTR(KKLINE2)
       IF(KMATCH(KKLINE2,KMATCH1)) THEN
          READ(KKLINE1,*,ERR=8000) NRXMAX
          IMATCH=IMATCH+1
       ELSE IF(KMATCH(KKLINE2,KMATCH2)) THEN
          READ(KKLINE1,*,ERR=8000) NTXMAX
          IMATCH=IMATCH+1
       ELSE IF(KMATCH(KKLINE2,KMATCH3)) THEN
          READ(KKLINE1,*,ERR=8000) NRXMAX
          IMATCH=IMATCH+1
       ENDIF
    END DO

    BACKSPACE 15

    CALL NDINIT(NDPOS)
    DO NRX=1,NRXMAX
       CALL NDREAD(NDPOS,15,R(NRX),IERR)
       IF(IERR /= 0) GOTO 8000
    ENDDO
    DO NTX=1,NTXMAX
       CALL NDREAD(NDPOS,15,T(NTX),IERR)
       IF(IERR /= 0) GOTO 8000
    ENDDO
    DO NTX=1,NTXMAX
       DO NRX=1,NRXMAX
          CALL NDREAD(NDPOS,15,F2(NRX,NTX),IERR)
          IF(IERR /= 0) GOTO 8000
       ENDDO
    ENDDO

    CLOSE(15)
    WRITE(6,"(' ',3A,2(I3,A))") ' - READ ASCII FILE :',KFID,'(',NRXMAX,',',NTXMAX,')'

    IF(MODE == 0) THEN
       OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='NEW')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
          RETURN
       ENDIF
       WRITE(15) NRXMAX,NTXMAX
       WRITE(15) (R(NRX),NRX=1,NRXMAX)
       WRITE(15) (T(NTX),NTX=1,NTXMAX)
       WRITE(15) ((F2(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)
       WRITE(6,*) ' - WRITE FILE:',KFILEB(1:60)
       CLOSE(15)
    ENDIF

9000 IF(KFID(1:6) == 'VOLUME'.OR.KFID(1:4) == 'BPOL'.OR. &
    &   KFID(1:4) == 'SURF'  .OR.KFID(1:6) == 'RMINOR') THEN
       IF(R(1) /= 0.D0) THEN
          NRXMAX=NRXMAX+1
          R(2:NRXMAX) = R(1:NRXMAX-1)
          R(1) = 0.D0
          DO NTX=1,NTXMAX
             F2(2:NRXMAX,NTX) = F2(1:NRXMAX-1,NTX)
             F2(1,NTX)=0.D0
          ENDDO
       ENDIF
    ENDIF
    RETURN

8000 WRITE(6,*) 'XX FILE READ ERROR'
    RETURN
  END SUBROUTINE TRXR2D
