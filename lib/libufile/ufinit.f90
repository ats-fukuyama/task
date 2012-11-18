MODULE ufinit

  PUBLIC

  ! global variables for ufile modules
  INTEGER, PARAMETER :: rkind = SELECTED_REAL_KIND(12,100)
  INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

  INTEGER, PARAMETER :: N0DMAX = 180 ! 2008 Public Release

  INTEGER(ikind),DIMENSION(N0DMAX) :: NID_ND

  CHARACTER(LEN=80) :: KUFDIR ! database directory path
  CHARACTER(LEN=80) :: KUFDEV ! device name
  CHARACTER(LEN=80) :: KUFDCG ! shot number
  CHARACTER(LEN=80) :: KDIRX  ! data storage directory path (***/in/)

CONTAINS

  SUBROUTINE ufile_init(KUF_DIR,KUF_DEV,KUF_DCG,IERR)
! --------------------------------------------------------------------
!  ***  Initialization routine.                                 ***
!  ***   - Locating ufile directory and checking its existence  ***
!  ***   - Substitution of module variables (Initialization)    ***
!
!  ***  This subrouitne must be called firstly every time you   ***
!  ***  refer to a data set of UFILE (i.e. every shot).         ***
!   
!   (This subroutine searches 'KUFDIR/DEVICE/SHOT/in'.)
!
!   KUF_DIR  : the directory path containing UFILE database
!   KUF_DEV  : the device name
!   KUF_DCG  : the discharge number
!   IERR     : error identifier
! --------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=80), INTENT(IN)  :: KUF_DIR, KUF_DEV, KUF_DCG
    INTEGER(ikind),    INTENT(OUT) :: IERR

    LOGICAL            :: DIR, FILE
    CHARACTER(LEN=100) :: KFILE
    INTEGER(ikind)     :: IKNDEV, IKNDCG, IKDIRX, IST

    ! substitution to module variables -----
    KUFDIR = KUF_DIR
    KUFDEV = KUF_DEV
    KUFDCG = KUF_DCG
    ! --------------------------------------

    IERR   = 0
    IKNDEV = LEN_TRIM(KUFDEV)
    IKNDCG = LEN_TRIM(KUFDCG)

    KDIRX = TRIM(KUFDIR)//'/'//KUFDEV(1:IKNDEV)//'/' &
         &                   //KUFDCG(1:IKNDCG)//'/in/'

    INQUIRE(FILE=KDIRX,EXIST=DIR,IOSTAT=IST)

    IF(IST > 0) THEN
       WRITE(6,*) ' XX Directory inquiring error. (ufile_dir_inquire)'
       WRITE(6,*) ' XX IOSTAT= ',IST
       IERR = 2
       RETURN
    ENDIF

    ! Some compilers cannot recognize whether the directory exists or not
    !  (always says "FALSE"), so that NE ufile should be checked
    !  in that case because it is the most possible existence file
    !  in a set of ufiles
    IKDIRX = LEN_TRIM(KDIRX)
    IF(DIR.EQV. .FALSE.) THEN
       KFILE=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'2d'//KUFDCG(1:IKNDCG)//'.NE'
       INQUIRE(FILE=KFILE,EXIST=FILE,IOSTAT=IST)
       IF(IST > 0) THEN
          WRITE(6,*) ' XX ufile_dir_inquire: File inquiring error.'
          WRITE(6,*) ' XX IOSTAT= ',IST
          IERR = 2
          RETURN
       END IF
       IF(FILE.EQV..FALSE.) THEN
          WRITE(6,'(A)') &
               ' XX ufile_dir_inquire: Designated directory does not exist!'
          WRITE(6,'(A7,A80)') '   ==> ',KDIRX
          IERR = 1
          RETURN
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE ufile_init


  SUBROUTINE ufile_inquire(KUFDIM,KFID,KFILE,KFILEB,IDBIN,ERROUT,IERR)
! -------------------------------------------------------------------------
!   *** inquire ufile and check its existence ***
!   < input >
!   KUFDIM : the dimension of ufile data; 1d or 2d
!   KFID   : the name of variables; AMIN, PNBI, NE, ...etc.
!   ERROUT : Error message output switch
!          : = 0  write error message to standard output 
!          : = 1  suppress error message output
!   
!   < output >
!   KFILE  : file name of a varible data including its path
!   KFILEB : binary file name of a varible data including its path
!   IDBIN  : the parameter which determines how to handle exp. files.
!   IDBIN =   0 : Binary files are loaded if available, or ASCII files
!                 are loaded and aftermath binary files are created.
!             1 : Only binary files are loaded.
!             2 : Only ASCII files are loaded and binary files are NOT
!                  created.
!   IERR   : error identifier
! -------------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=10), INTENT(IN)    :: KUFDIM, KFID
    INTEGER(ikind),    INTENT(IN)    :: ERROUT
    CHARACTER(LEN=100),INTENT(OUT)   :: KFILE,KFILEB
    INTEGER(ikind),    INTENT(INOUT) :: IDBIN
    INTEGER(ikind),    INTENT(OUT)   :: IERR

    LOGICAL            :: LEX
    CHARACTER(LEN=100) :: KDIRR1
    INTEGER(ikind)         :: IKNDEV,IKNDCG,IKNDIM,IKDIRX,KL,IST

    IERR   = 0
    IKNDEV = len_trim(KUFDEV)
    IKNDCG = len_trim(KUFDCG)
    IKNDIM = len_trim(KUFDIM)
    IKDIRX = len_trim(KDIRX)

    KDIRR1 = KDIRX(1:IKDIRX) // KUFDEV(1:IKNDEV) //  &
             KUFDIM(1:IKNDIM)// KUFDCG(1:IKNDCG)

    KL     = len_trim(KDIRR1)
    KFILE  = KDIRR1(1:KL)//'.'//KFID
    KL     = len_trim(KFILE)
    KFILEB = KFILE(1:KL)//'.bin'

    ! ASCII file
    INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=IST)
    IF(LEX .EQV. .FALSE.)THEN
       IF(ERROUT == 0)THEN
          WRITE(6,*) ' XX ufile_inquire: KFILE does not exist.'
          WRITE(6,*) ' XX KFILE: ',TRIM(KFILE)
       END IF
       IERR = -1
       RETURN
    ELSE IF(IST /= 0) THEN
       IF(ERROUT == 0)THEN
          WRITE(6,*) ' XX ufile_inquire: KFILE inquiring error. IOSTAT= ',IST
          WRITE(6,*) ' XX KFILE: ',TRIM(KFILE)
       END IF
       IERR = -1
       RETURN
    ENDIF

    IF(IDBIN < 2)THEN
       ! Binary file
       INQUIRE(FILE=KFILEB,EXIST=LEX,IOSTAT=IST)
       IF(LEX .EQV. .FALSE.)THEN
          IF(IDBIN==1)THEN
             IF(ERROUT == 0)THEN
                WRITE(6,*) ' XX ufile_inquire: KFILEB does not exist.'
                WRITE(6,*) ' XX KFILEB: ',TRIM(KFILEB)
             END IF
             IERR = -1
             RETURN
          END IF
          IERR = 0
          RETURN
       ELSE IF(IST /= 0) THEN
          IF(ERROUT == 0)THEN
             WRITE(6,*) ' XX ufile_inqurie: KFILEB inquiring error. IOSTAT= ',IST
             WRITE(6,*) ' XX KFILEB: ',TRIM(KFILEB)
          END IF
          IF(IDBIN==1)THEN
             IERR = -1
             RETURN
          END IF
          IERR = 1
          RETURN
       ENDIF
       IDBIN = 1 ! read binary file because of the existence of it
    END IF
    IERR = 0

    RETURN
  END SUBROUTINE ufile_inquire

END MODULE ufinit
