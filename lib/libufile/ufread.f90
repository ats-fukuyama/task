MODULE ufread

  USE ufinit,ONLY: rkind,ikind,KUFDIR,KUFDEV,KUFDCG

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC ufread_0d,ufread_1d_time,ufread_2d_time,ufread_2d_errbar
  
CONTAINS

  SUBROUTINE ufread_0d(NDMAX,IERR)
! -----------------------------------------------------------------------
!   ***  0D UFILE reader ***
!
!   This subroutine reads the ***_0d.dat file of UFILE datasets.
!
!   The subroutine 'uflist_init' initialize the structure 'uf0d'
!    to contain variable names and its data type.
!
!   The subroutine 'uflist_set' contains the variable data 
!    into the structure 'ud0d'
!
!   < output >
!   NDMAX  : the number of data included in the _0d.dat file
!   IERR   : error identifier
! -----------------------------------------------------------------------
    USE uflist, ONLY: uflist_init,uflist_set
    USE ufinit, ONLY: kufdir,kufdev,kufdcg,n0dmax
    IMPLICIT NONE

    INTEGER(ikind)             ,INTENT(OUT) :: NDMAX, IERR

    CHARACTER(LEN=2000)                     :: KKLINE
    CHARACTER(LEN=100)                      :: KFILE
    CHARACTER(LEN=80)                       :: KFDAT
    CHARACTER(LEN=10)                       :: KFVER,KUFDIM
    LOGICAL,              DIMENSION(N0DMAX) :: LEX0
    REAL(rkind),          DIMENSION(N0DMAX) :: FR0
    INTEGER(ikind),       DIMENSION(N0DMAX) :: FI0
    CHARACTER(LEN=15),    DIMENSION(N0DMAX) :: FC0
    INTEGER(ikind) :: IKNDEV,IKNDCG,NDTMAX,IST,UNIT
    LOGICAL    :: LEX

    IERR   = 0
    NDMAX  = 0
    UNIT   = 20

    KUFDIM = '0d'
    KFVER  = 'pr08'

    IKNDEV = LEN_TRIM(KUFDEV)
    IKNDCG = LEN_TRIM(KUFDCG)

    KFDAT = TRIM(KFVER)//'_'//TRIM(KUFDEV)//'_'// &
            TRIM(KUFDCG)//'_'//TRIM(KUFDIM)//'.dat'
    KFILE = TRIM(KUFDIR)//'/'//KUFDEV(1:IKNDEV)//'/'//KUFDCG(1:IKNDCG)// &
            '/'//TRIM(KFDAT)

    INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=IST)
    IF(LEX .EQV. .FALSE.)THEN
       WRITE(6,*) 'XX ufread_0d: _0d.dat file does not exist.'
       WRITE(6,*) 'XX KFILE: ',TRIM(KFILE)
       IERR = 1
       RETURN
    ELSE IF(IST /= 0) THEN
       WRITE(6,*) 'XX ufread_0d: _0d.dat inquire error. IOSTAT= ',IST
       WRITE(6,*) 'XX KFILE: ',TRIM(KFILE)
       IERR = 2
       RETURN
    ENDIF

    OPEN(UNIT,FILE=KFILE,IOSTAT=IST)
    IF(IST /= 0)THEN
       WRITE(6,*) 'XX ufread_0d: _0d.dat file open error. IOSTAT= ',IST
       WRITE(6,*) 'XX KFILE: ',TRIM(KFILE)
       IERR = 2
       RETURN
    END IF

    CALL uflist_init

    ! read labels 
    READ(UNIT,'(A2000)',IOSTAT=IST) KKLINE
    IF(IST /= 0) THEN
       WRITE(6,*) 'XX ufread_0d: 0d.dat read errors.'
       WRITE(6,*) 'XX KFILE: ',TRIM(KFILE)
       RETURN
    END IF
    IF(SCAN(KKLINE,',') == 0)THEN
       ! traditional format (TAB format)
       CALL ufread_0d_pr98(UNIT,NDTMAX,NDMAX,LEX0,FR0,FI0,FC0,IERR)
    ELSE
       ! up-to-date format (csv format)
       CALL ufread_0d_pr08(UNIT,NDTMAX,NDMAX,LEX0,FR0,FI0,FC0,IERR)
    END IF

    IF(IERR /= 0)THEN
       WRITE(6,*) 'XX ufread_0d: 0d.dat read errors. IERR= ',IERR
       WRITE(6,*) 'XX KFILE: ',TRIM(KFILE)
       RETURN
    END IF

    CLOSE(UNIT)

    CALL uflist_set(LEX0,FR0,FI0,FC0)

    WRITE(6,'(1X,3A,I2,A,I3,A)') ' - Read 0D data: ',TRIM(KFDAT), &
                                 '   (',NDTMAX,',',NDMAX,')'

    RETURN
  END SUBROUTINE ufread_0d


  SUBROUTINE ufread_0d_pr08(UNIT,NDTMAX,NDMAX,LEX0,FR0,FI0,FC0,IERR)
! -------------------------------------------------------------------------
!  0-D UFILE reader for the up-to-date format (csv format)
! -------------------------------------------------------------------------
    USE uflist, ONLY: uf0d
    USE ufinit, ONLY: n0dmax, nid_nd
    IMPLICIT NONE

    INTEGER(ikind), INTENT(IN)  :: UNIT
    INTEGER(ikind), INTENT(OUT) :: NDMAX,NDTMAX,IERR

    LOGICAL,          DIMENSION(N0DMAX),INTENT(OUT) :: LEX0
    REAL(rkind),      DIMENSION(N0DMAX),INTENT(OUT) :: FR0
    INTEGER(ikind),   DIMENSION(N0DMAX),INTENT(OUT) :: FI0
    CHARACTER(LEN=15),DIMENSION(N0DMAX),INTENT(OUT) :: FC0

    CHARACTER(LEN=2000) :: KKLINE
    INTEGER(ikind)      :: ND,NID,NDT,NDMAX_LABEL,IKNFID,EOL,IST

    IERR = 0
    REWIND UNIT

    ! --- read labels ---
    READ(UNIT,'(A2000)',IOSTAT=IST) KKLINE
    IF(IST /= 0) THEN
       IERR = 2008
       RETURN
    END IF
    
    NID_ND(1:N0DMAX) = 0
    ND  = 0
    EOL = 0
    label: DO
       KKLINE = TRIM(KKLINE)
       IKNFID = INDEX(KKLINE,',')
       IF(IKNFID == 0) THEN ! the last item of the list
          IKNFID = 15 + 1
          EOL = 1
       END IF

       DO NID = 1, N0DMAX
          IF(KKLINE(1:IKNFID-1)==uf0d(NID)%kfid)THEN
             ND = ND + 1

             LEX0(NID)  = .TRUE.
             NID_ND(ND) = NID
             EXIT
          END IF
       END DO

       IF(EOL == 1)THEN ! the end of line
          NDMAX_LABEL = ND
          EXIT
       END IF
       
       KKLINE(1:) = KKLINE(IKNFID+1:)
    END DO label

    ! --- initialize ---
    FR0(1:N0DMAX) = 0.d0
    FI0(1:N0DMAX) = 0
    FC0(1:N0DMAX) = '---'

    ! --- read values ---
    NDT = 0
    time_slice: DO
       READ(UNIT,'(A2000)',IOSTAT=IST) KKLINE
       IF(IST > 0) THEN
          IERR = 2008
          RETURN
       ELSE IF(IST < 0)THEN ! the end of file
          NDTMAX = NDT
          EXIT
       END IF

       ND  = 0
       EOL = 0
       value: DO
          ND = ND + 1
          KKLINE = TRIM(KKLINE)
          IKNFID = INDEX(KKLINE,',')
          IF(IKNFID == 0)THEN ! the last item of the list
             IKNFID = 15 + 1
             EOL = 1
          END IF

          NID = NID_ND(ND)
          IF(NID /= 0)THEN
             SELECT CASE(uf0d(NID)%id_type)
             CASE(1) ! REAL
                READ(KKLINE(1:IKNFID-1),*) FR0(NID)
             CASE(2) ! INTEGER
                READ(KKLINE(1:IKNFID-1),*) FI0(NID)
             CASE(3) ! CHARACTER
                FC0(NID) = KKLINE(1:IKNFID-1)
             END SELECT
          END IF

          IF(EOL == 1)THEN ! the end of line
             NDMAX = ND
             EXIT
          END IF
          KKLINE(1:) = KKLINE(IKNFID+1:)
       END DO value

       NDT = NDT + 1

       IF(NDMAX /= 0)THEN
          IF(NDMAX /= NDMAX_LABEL)THEN
             WRITE(6,*) 'XX ufread_0d_pr08: 0d.dat has some error.'
             WRITE(6,'(1X,A,I3,A,I4,A,I4)') 'XX NDT= ',NDT,   &
              &             'NDMAX= ',NDMAX,'NDMAX_LABEL= ',NDMAX_LABEL
             IERR = 2008
             RETURN
          END IF
       END IF
    END DO time_slice

    RETURN
  END SUBROUTINE ufread_0d_pr08


  SUBROUTINE ufread_0d_pr98(UNIT,NDTMAX,NDMAX,LEX0,FR0,FI0,FC0,IERR)
! -----------------------------------------------------------------------
!  0-D UFILE reader for the traditional format
! -----------------------------------------------------------------------
    USE uflist, ONLY: uf0d
    USE ufinit, ONLY: n0dmax, nid_nd
    IMPLICIT NONE

    INTEGER(ikind), INTENT(IN)  :: UNIT
    INTEGER(ikind), INTENT(OUT) :: NDMAX,NDTMAX,IERR

    LOGICAL,          DIMENSION(N0DMAX),INTENT(OUT) :: LEX0
    REAL(rkind),      DIMENSION(N0DMAX),INTENT(OUT) :: FR0
    INTEGER(ikind),   DIMENSION(N0DMAX),INTENT(OUT) :: FI0
    CHARACTER(LEN=15),DIMENSION(N0DMAX),INTENT(OUT) :: FC0

    CHARACTER(LEN=80) :: KKLINE
    CHARACTER(LEN=10) :: KKVAL
    INTEGER(ikind)      :: NL,NL1,ND,NID,NDT,NDMAX_LABEL,NLMAX, &
                           IKNFID,IKNVAL,EOL,EOF,IST

    IERR = 0
    REWIND UNIT

    IKNFID = 10 ! fixed field length

    NID_ND(1:N0DMAX) = 0
    ND  = 0
    NL  = 0
    EOL = 0

    ! --- read labels ---
    label: DO
       READ(UNIT,'(A80)',IOSTAT=IST) KKLINE
       IF(IST /= 0) THEN
          IERR = 1998
          RETURN
       END IF
       NL = NL + 1 ! the number of lines in a data block

       line1: DO NL1 = 1, 7
          ND = ND + 1 ! the number of signals
          KKLINE = TRIM(KKLINE)

          DO NID = 1, N0DMAX
             IF(TRIM(ADJUSTL(KKLINE(1:IKNFID)))==uf0d(NID)%kfid)THEN
                LEX0(NID)  = .TRUE.
                NID_ND(ND) = NID
                EXIT
             END IF
          END DO
       
          KKLINE(1:) = KKLINE(IKNFID+2:)
          IF(NL1 /= 7 .AND. LEN_TRIM(KKLINE) == 0)THEN ! the end of block
             EOL = 1
             NDMAX_LABEL = ND
             NLMAX       = NL
             EXIT
          END IF
       END DO line1

       IF(EOL==1) EXIT
    END DO label

    ! --- initialize ---
    FR0(1:N0DMAX) = 0.d0
    FI0(1:N0DMAX) = 0
    FC0(1:N0DMAX) = '---'


    ! --- read value ---
    NDT = 0
    EOF = 0
    time_slice: DO

       ND  = 0
       EOL = 0
       value: DO
          READ(UNIT,'(A80)',IOSTAT=IST) KKLINE
          IF(IST > 0) THEN
             IERR = 1998
             RETURN
          ELSE IF(IST < 0)THEN ! the end of file
             NDTMAX = NDT
             EOF = 1
             EXIT
          END IF

          line2: DO NL1 = 1, 7
             ND = ND + 1
             KKLINE = TRIM(KKLINE)

             NID = NID_ND(ND)
             IF(NID /= 0)THEN
                KKVAL = ADJUSTL(KKLINE(1:IKNFID))
                IKNVAL = LEN_TRIM(KKVAL)
                SELECT CASE(uf0d(NID)%id_type)
                CASE(1) ! REAL
                   READ(KKVAL(1:IKNVAL),*) FR0(NID)
                CASE(2) ! INTEGER
                   READ(KKVAL(1:IKNVAL),*) FI0(NID)
                CASE(3) ! CHARACTER
                   FC0(NID) = ADJUSTL(KKLINE(1:IKNFID))
                END SELECT
             END IF

             KKLINE(1:) = KKLINE(IKNFID+2:)
             IF(NL1 /= 7 .AND. LEN_TRIM(KKLINE) == 0)THEN ! the end of block
                EOL   = 1
                NDMAX = ND
                EXIT
             END IF

          END DO line2

          IF(EOL==1) EXIT
       END DO value

       NDT = NDT + 1

       IF(NDMAX /= 0)THEN
          IF(NDMAX /= NDMAX_LABEL)THEN
             WRITE(6,*) 'XX ufread_0d_pr98: 0d.dat has some error.'
             WRITE(6,'(1X,A,I3,A,I4,A,I4)') 'XX NDT= ',NDT,   &
              &             'NDMAX= ',NDMAX,'NDMAX_LABEL= ',NDMAX_LABEL
             IERR = 1998
             RETURN
          END IF
       END IF

       ! Skip lines of the labels
       DO NL = 1, NLMAX
          READ(UNIT,'()',IOSTAT=IST)
          IF(IST < 0)THEN ! the end of file
             NDTMAX = NDT
             EOF    = 1
             EXIT
          END IF
       END DO

       IF(EOF==1) EXIT
    END DO time_slice

    RETURN
  END SUBROUTINE ufread_0d_pr98

! =======================================================================
!  UFILE reader subroutines
!     - ufread_1d_time
!     - ufread_2d_time
!     - ufread_2d_errbar
!
!  input:
!
!  KDIRX    : Directory
!  KUFDEV   : Device name in UFILEs
!  KUFDCG   : Shot number in UFILEs
!  KFID     : Variable name in UFILEs
!  NRUM     : Array size (radial)
!  NTUM     : Array size (time)
!  NTXMAX   : Maximum Number of the Time Mesh
!  ERROUT   : Switch for writing inquire error message to standard outoput
!              = 0: write, = 1(else): suppress
!
!  output:
!
!  RUF(NRUM)      : Equally Spaced Normalized Radial Data
!  TMU(NTUM)      : Total Time Data (The Number of DT)
!  F1(NTUM)       : Functional Values
!  F2(NTUM,NRUM)  : Functional Values
!  NRFMAX         : Maximum Number of the Radial Mesh
!  IERR           : Error Indicator
! =======================================================================

  SUBROUTINE ufread_1d_time(KFID,TMU,F1,NTXMAX,NTUM,ID_BIN,ERROUT,IERR)
! ----------------------------------------------------------------------
!    1D UFILE reader (Time Evolution UFILE)
! ----------------------------------------------------------------------
    USE ufinit, ONLY: ufile_inquire
    IMPLICIT NONE

    CHARACTER(LEN=10),            INTENT(IN)    :: KFID
    INTEGER(ikind),               INTENT(IN)    :: NTUM, ID_BIN, ERROUT
    REAL(rkind), DIMENSION(NTUM), INTENT(OUT)   :: TMU
    REAL(rkind), DIMENSION(NTUM), INTENT(OUT)   :: F1
    INTEGER(ikind),               INTENT(OUT)   :: NTXMAX,IERR

    CHARACTER(LEN=100) :: KFILE,KFILEB
    CHARACTER(LEN=10)  :: KUFDIM
    INTEGER(ikind)     :: IDBIN

    KUFDIM = '1d'
    IERR   = 0

    IDBIN  = ID_BIN

    CALL ufile_inquire(KUFDIM,KFID,KFILE,KFILEB,IDBIN,ERROUT,IERR)
    SELECT CASE(IERR)
    CASE(-1) ! fatal error: file not found or broken
       F1(1:NTUM) = 0.D0
       IERR = 1
       RETURN
    CASE(0)
       CONTINUE
    CASE(1) ! binary file is broken. read ASCII file. (if IDBIN=0)
       WRITE(6,*) '## Binary file has some error.'       
       WRITE(6,*) '## read ASCII file and re-create the binary file.'
       IERR = 0
    END SELECT

    CALL uf_bin_1d_read(KFID,KFILE,KFILEB,TMU,F1,NTUM,NTXMAX,IDBIN)

    RETURN
  END SUBROUTINE ufread_1d_time


  SUBROUTINE ufread_2d_time(KFID,RUF,TMU,F2,NRFMAX,NTXMAX,NRUM,NTUM, &
                            ID_BIN,ERROUT,IERR)
! ----------------------------------------------------------------------
!   2D UFILE reader (Time Evolution UFILE)
! ----------------------------------------------------------------------
    USE ufinit, ONLY: ufile_inquire
    USE libitp
    IMPLICIT NONE

    CHARACTER(LEN=10),                 INTENT(IN)    :: KFID
    INTEGER(ikind),                    INTENT(IN)    :: NRUM,NTUM,ID_BIN,ERROUT
    INTEGER(ikind),                    INTENT(OUT)   :: NRFMAX,NTXMAX,IERR
    REAL(rkind), DIMENSION(NRUM),      INTENT(OUT)   :: RUF
    REAL(rkind), DIMENSION(NTUM),      INTENT(OUT)   :: TMU
    REAL(rkind), DIMENSION(NTUM,NRUM), INTENT(OUT)   :: F2

    INTEGER(ikind)                    :: IDBIN, NTX, NRF, IST
    REAL(rkind), DIMENSION(NTUM)      :: F2CTR, F2EDG
    REAL(rkind), DIMENSION(NRUM,NTUM) :: F2I
    CHARACTER(LEN=100)                :: KFILE,KFILEB
    CHARACTER(LEN=10)                 :: KUFDIM

    KUFDIM = '2d'
    IERR   = 0

    IDBIN  = ID_BIN

    CALL ufile_inquire(KUFDIM,KFID,KFILE,KFILEB,IDBIN,ERROUT,IERR)
    SELECT CASE(IERR)
    CASE(-1) ! fatal error: file not found or broken
       F2(1:NTUM,1:NRUM) = 0.D0
       IERR = 1
       RETURN
    CASE(0)
       CONTINUE
    CASE(1) ! binary file is broken. read ASCII file. (if IDBIN=0)
       WRITE(6,*) '## Binary file has some error.'
       WRITE(6,*) '## read ASCII file and re-create the binary file.'
       IERR = 0
    END SELECT


    CALL uf_bin_2d_read(KFID,KFILE,KFILEB,TMU,RUF,F2I,NRUM,NTUM,NRFMAX,NTXMAX,IDBIN)
    IF(NRFMAX == 1) THEN
       WRITE(6,*) '## ufread_2d_time: '
       WRITE(6,*) &
            '## ',KFID,'has not been read due to lacking radial point data.'
       F2(1:NTUM,1:NRUM)=0.D0
       IERR=2
       RETURN
    ENDIF

    F2(1:NTXMAX,1:NRFMAX) = transpose(F2I(1:NRFMAX,1:NTXMAX))


    ! -------------------------------------------------------------------
    !  For interpolation, center and/or edge value is added if necessary. 
    !  The variables with center value having to be zero such as 
    !   BPOL, RMINOR, SURF and VOLUME have been manupulated beforehand
    !   in 'uf_bin_2d_read' subroutine.
    ! -------------------------------------------------------------------

    IF(RUF(1) /= 0.D0 .AND. RUF(NRFMAX) /= 1.D0) THEN
       DO NTX = 1, NTXMAX
          F2CTR(NTX) = FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
          F2EDG(NTX) = AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1), &
               &                     F2(NTX,NRFMAX-2),RUF(NRFMAX),    &
               &                     RUF(NRFMAX-1),RUF(NRFMAX-2))
       ENDDO

       ! substitution
       NRFMAX = NRFMAX + 2
       RUF(2:NRFMAX-1) = RUF(1:NRFMAX-2)
       RUF(1)          = 0.D0
       RUF(NRFMAX)     = 1.D0
       F2(1:NTXMAX,2:NRFMAX-1) = F2(1:NTXMAX,1:NRFMAX-2)
       F2(1:NTXMAX,1)          = F2CTR(1:NTXMAX)
       F2(1:NTXMAX,NRFMAX)     = F2EDG(1:NTXMAX)

    ELSE IF(RUF(1) /= 0.D0 .AND. RUF(NRFMAX) == 1.D0) THEN
       DO NTX=1,NTXMAX
          F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
       ENDDO

       ! substitution
       NRFMAX = NRFMAX + 1
       RUF(2:NRFMAX) = RUF(1:NRFMAX-1)
       RUF(1)        = 0.D0
       F2(1:NTXMAX,2:NRFMAX) = F2(1:NTXMAX,1:NRFMAX-1)
       F2(1:NTXMAX,1)        = F2CTR(1:NTXMAX)


    ELSEIF(RUF(1) == 0.D0 .AND. RUF(NRFMAX) /= 1.D0) THEN
       DO NTX=1,NTXMAX
          F2EDG(NTX)=AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1),     &
               &                          F2(NTX,NRFMAX-2),RUF(NRFMAX), &
               &                          RUF(NRFMAX-1),RUF(NRFMAX-2))
       ENDDO

       ! substitution
       NRFMAX = NRFMAX + 1
       RUF(NRFMAX)         = 1.D0
       F2(1:NTXMAX,NRFMAX) = F2EDG(1:NTXMAX)

    ELSEIF(RUF(1) == 0.D0 .AND. RUF(NRFMAX) == 1.D0) THEN
       CONTINUE
    ELSE
       STOP 'XX ufread_2d_time: Radial point data has some error.'
    ENDIF

    CALL data_error_correct(KFID,RUF,F2,NTXMAX,NRUM,NTUM)

    RETURN
  END SUBROUTINE ufread_2d_time


  SUBROUTINE ufread_2d_errbar(KFID,RUF,TMU,F2,NRFMAX,NTXMAX, &
                              NRUM,NTUM,ID_BIN,ERROUT,IERR)
! -----------------------------------------------------------------------
!   error bar
! ----------------------------------------------------------------------
    USE ufinit, ONLY: ufile_inquire
    IMPLICIT NONE

    CHARACTER(LEN=10),     INTENT(IN) :: KFID
    INTEGER(ikind),        INTENT(IN) :: NRUM, NTUM, ID_BIN, ERROUT
    INTEGER(ikind),                    INTENT(OUT)   :: NRFMAX,NTXMAX,IERR
    REAL(rkind), DIMENSION(NRUM),      INTENT(OUT)   :: RUF
    REAL(rkind), DIMENSION(NTUM),      INTENT(OUT)   :: TMU
    REAL(rkind), DIMENSION(NTUM,NRUM), INTENT(OUT)   :: F2

    INTEGER(ikind) :: IKNDEV, IKNDCG, IKDIRX, IDBIN, NRF, IST
    REAL(rkind), DIMENSION(NRUM,NTUM) :: F2I
    CHARACTER(LEN=100)                :: KFILE,KFILEB
    CHARACTER(LEN=10)                 :: KUFDIM
    LOGICAL :: LEX

    KUFDIM = '2d'
    IERR   = 0

    IDBIN  = ID_BIN

    CALL ufile_inquire(KUFDIM,KFID,KFILE,KFILEB,IDBIN,ERROUT,IERR)
    SELECT CASE(IERR)
    CASE(-1) ! fatal error; file not found or broken
       F2(1:NTUM,1:NRUM) = 0.D0
       IERR = 1
       RETURN
    CASE(0)
       CONTINUE
    CASE(1) ! binary file is broken. read ASCII file. (if IDBIN=0)
       WRITE(6,*) '## Binary file has some error.'
       WRITE(6,*) '## read ASCII file and re-create the binary file.'
       IERR = 0
    END SELECT


    CALL uf_bin_2d_read(KFID,KFILE,KFILEB,TMU,RUF,F2I,NRUM,NTUM,NRFMAX,NTXMAX,IDBIN)
    IF(NRFMAX == 1) THEN
       WRITE(6,*) '## ufread_2d_errbar: '
       WRITE(6,*) &
            '## ',KFID,'has not been read due to lacking radial point data.'
       F2(1:NTUM,1:NRUM)=0.D0
       IERR=2
       RETURN
    ENDIF

    F2(1:NTXMAX,1:NRFMAX) = transpose(F2I(1:NRFMAX,1:NTXMAX))


    !  ----- for the time being -----
    IF(RUF(NRFMAX) > 1.D0) THEN
       DO NRF=1,NRFMAX
          IF(RUF(NRF) > 1.D0) THEN
             NRFMAX=NRF-1
             EXIT
          ENDIF
       ENDDO
    ENDIF
    IF(KFID == 'NEXP'.OR.KFID == 'NEXPEB') THEN
       F2(1:NTXMAX,1:NRUM)=F2(1:NTXMAX,1:NRUM)*1.D-20
    ELSE
       F2(1:NTXMAX,1:NRUM)=F2(1:NTXMAX,1:NRUM)*1.D-3
    ENDIF

    RETURN
  END SUBROUTINE ufread_2d_errbar


! ************************************************************************
! ************************************************************************

  SUBROUTINE uf_bin_1d_read(KFID,KFILE,KFILEB,T,F1,NTM,NTXMAX,IDBIN)
! ----------------------------------------------------------------------
!   Read 1D variable from ASCII file or binary file of UFILE
!   Make binary file in order to faster reading of data if necessary
!
!   This subroutine is called only from the subroutine 'ufread_1d_time'
! ----------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=10),           INTENT(IN)    :: KFID
    CHARACTER(LEN=100),          INTENT(IN)    :: KFILE,KFILEB
    INTEGER(ikind),              INTENT(IN)    :: NTM
    INTEGER(ikind),              INTENT(INOUT) :: IDBIN
    INTEGER(ikind),              INTENT(OUT)   :: NTXMAX
    REAL(rkind), DIMENSION(NTM), INTENT(OUT)   :: T, F1

    CHARACTER(LEN=15) :: KFIDB
    INTEGER(ikind)    :: UNIT,IST,IERR
    INTEGER(ikind)    :: NTX, NDPOS, NT
    LOGICAL           :: LEX

    IERR = 0
    UNIT = 15

    bin_ascii: SELECT CASE(IDBIN)
    CASE(1) ! read binary file
       OPEN(UNIT,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='OLD')
       IF(IST > 0)THEN
          WRITE(6,*) 'XX uf_bin_1d_read: binary file open error. IOSTAT= ',IST
          WRITE(6,*) 'XX  failed to read the binary file'
          WRITE(6,*) 'XX KFILEB: ',TRIM(KFILEB)
          NTXMAX = 0
          RETURN
       END IF
       READ(UNIT) NTXMAX
       READ(UNIT) (T(NTX) ,NTX = 1, NTXMAX)
       READ(UNIT) (F1(NTX),NTX = 1, NTXMAX)
       CLOSE(UNIT)
       KFIDB=KFID//'.bin'
       WRITE(6,"(' ',3A,I3,A)") ' - Read Binary file :',TRIM(KFIDB),' (',NTXMAX,')'
       RETURN

    CASE DEFAULT ! read ASCII file
       OPEN(UNIT,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',STATUS='OLD')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX uf_bin_1d_read: KFILE open error. IOSTAT= ',IST
          WRITE(6,*) 'XX  failed to read the ASCII file'
          WRITE(6,*) 'XX KFILE= ',TRIM(KFILE)
          NTXMAX=0
          RETURN
       ENDIF

       CALL uf_head_search(UNIT,NTXMAX,NTXMAX,NTXMAX,IERR)
       IF(IERR /= 0)THEN
          WRITE(6,*) 'XX uf_bin_1d_read: Read error of the header. IERR= ',IERR
          WRITE(6,*) 'XX KFILE= ',TRIM(KFILE)
       END IF

       CALL NDINIT(NDPOS)
       DO NTX=1,NTXMAX
          CALL NDREAD(NDPOS,UNIT,T(NTX),IERR)
          IF(IERR /= 0)THEN
             WRITE(6,*) 'XX uf_bin_1d_read: Time data read error. NTX= ',NTX
             RETURN
          END IF
       ENDDO
       DO NTX=1,NTXMAX
          CALL NDREAD(NDPOS,UNIT,F1(NTX),IERR)
          IF(IERR /= 0)THEN
             WRITE(6,*) 'XX uf_bin_1d_read: Function data read error. NTX= ',NTX
             RETURN
          END IF
       ENDDO

       CLOSE(UNIT)
       WRITE(6,"(' ',3A,I3,A)") ' - Read ASCII file :',TRIM(KFID),' (',NTXMAX,')'

    END SELECT bin_ascii

    ! write data into binary file (.bin)
    IF(IDBIN == 0) THEN
       OPEN(UNIT,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='REPLACE')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX uf_bin_1d_read: binary file open error. IOSTAT= ',IST
          WRITE(6,*) 'XX   failed to write the binary file'
          WRITE(6,*) 'XX KFILEB: ',TRIM(KFILEB)
          RETURN
       ENDIF
       WRITE(UNIT) NTXMAX
       WRITE(UNIT) (T(NT),NT=1,NTXMAX)
       WRITE(UNIT) (F1(NT),NT=1,NTXMAX)
       CLOSE(UNIT)
       WRITE(6,*) ' - Write binary file:',TRIM(KFILEB)
    ENDIF

    RETURN
  END SUBROUTINE uf_bin_1d_read


  SUBROUTINE uf_bin_2d_read(KFID,KFILE,KFILEB,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,IDBIN)
! ----------------------------------------------------------------------
!   Read 2D variable from ASCII file or binary file
!   Make binary file in order to faster reading of data if necessary
!
!   This subroutine is called only from the subroutine 'ufread_2d_time'
! ----------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=10),              INTENT(IN)    :: KFID
    CHARACTER(LEN=100),             INTENT(IN)    :: KFILE,KFILEB
    INTEGER(ikind),                 INTENT(IN)    :: NRM, NTM
    INTEGER(ikind),                 INTENT(INOUT) :: IDBIN
    INTEGER(ikind),                 INTENT(OUT)   :: NRXMAX, NTXMAX
    REAL(rkind), DIMENSION(NTM),    INTENT(OUT)   :: T
    REAL(rkind), DIMENSION(NRM),    INTENT(OUT)   :: R
    REAL(rkind), DIMENSION(NRM,NTM),INTENT(OUT)   :: F2

    CHARACTER(LEN=15) :: KFIDB
    INTEGER(ikind)    :: UNIT, IST, IERR
    INTEGER(ikind)    :: NRX, NTX, NDPOS, NT
    LOGICAL :: LEX

    IERR = 0
    UNIT = 15

    bin_ascii: SELECT CASE(IDBIN)
    CASE(1) ! read binary file
       OPEN(UNIT,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='OLD')
       IF(IST > 0)THEN
          WRITE(6,*) 'XX uf_bin_2d_read: binary file open error. IOSTAT= ',IST
          WRITE(6,*) 'XX KFILEB: ',TRIM(KFILEB)
          NTXMAX = 0
          RETURN
       END IF
       READ(UNIT) NRXMAX, NTXMAX
       READ(UNIT) (R(NRX), NRX = 1, NRXMAX)
       READ(UNIT) (T(NTX), NTX = 1, NTXMAX)
       READ(UNIT) ((F2(NRX,NTX), NRX = 1, NRXMAX), NTX = 1, NTXMAX)
       CLOSE(UNIT)
       KFIDB=KFID//'.bin'
       WRITE(6,"(' ',3A,2(I3,A))") ' - Read binary file :',TRIM(KFIDB),' (',NRXMAX,',',NTXMAX,')'

    CASE DEFAULT ! read ASCII file
       OPEN(UNIT,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',STATUS='OLD')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX uf_bin_2d_read: KFILE open error. IOSTAT= ',IST
          WRITE(6,*) 'XX KFILE= ',TRIM(KFILE)
          NTXMAX = 0
          RETURN
       ENDIF

       CALL uf_head_search(UNIT,NRXMAX,NTXMAX,NTXMAX,IERR)
       IF(IERR /= 0)THEN
          WRITE(6,*) 'XX uf_bin_2d_read: Read error of the header. IERR= ',IERR
          WRITE(6,*) 'XX KFILE= ',TRIM(KFILE)
       END IF

       CALL NDINIT(NDPOS)
       DO NRX=1,NRXMAX
          CALL NDREAD(NDPOS,UNIT,R(NRX),IERR)
          IF(IERR /= 0)THEN
             WRITE(6,*) 'XX uf_bin_2d_read: Radial data read error. NRX= ',NRX
             RETURN
          END IF
       ENDDO
       DO NTX=1,NTXMAX
          CALL NDREAD(NDPOS,UNIT,T(NTX),IERR)
          IF(IERR /= 0)THEN
             WRITE(6,*) 'XX uf_bin_1d_read: Time data read error. NTX= ',NTX
             RETURN
          END IF
       ENDDO
       DO NTX=1,NTXMAX
          DO NRX=1,NRXMAX
             CALL NDREAD(NDPOS,UNIT,F2(NRX,NTX),IERR)
             IF(IERR /= 0)THEN
                WRITE(6,*) &
      & 'XX uf_bin_1d_read: Function data read error. NTX= ',NTX,'NRX= ',NRX
                RETURN
             END IF
          ENDDO
       ENDDO
       
       CLOSE(UNIT)
       WRITE(6,"(' ',3A,2(I3,A))") ' - Read ASCII file :',TRIM(KFID),' (',NRXMAX,',',NTXMAX,')'

    END SELECT bin_ascii

    ! write data into binary file (.bin)
    IF(IDBIN == 0) THEN
       OPEN(UNIT,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',STATUS='REPLACE')
       IF(IST > 0) THEN
          WRITE(6,*) 'XX uf_bin_2d_read: binary file open error. IOSTAT= ',IST
          WRITE(6,*) 'XX KFILEB: ',TRIM(KFILEB)
          RETURN
       ENDIF
       WRITE(UNIT) NRXMAX, NTXMAX
       WRITE(UNIT) (R(NRX), NRX = 1, NRXMAX)
       WRITE(UNIT) (T(NTX), NTX = 1, NTXMAX)
       WRITE(UNIT) ((F2(NRX,NTX), NRX = 1, NRXMAX), NTX = 1, NTXMAX)
       WRITE(6,*) ' - Write binary file:',TRIM(KFILEB)
       CLOSE(UNIT)
    ENDIF

    ! insert the value at the axis
    IF(KFID(1:6) == 'VOLUME' .OR. KFID(1:4) == 'BPOL' .OR. &
    &  KFID(1:4) == 'SURF'   .OR. KFID(1:6) == 'RMINOR') THEN
       IF(R(1) /= 0.D0) THEN
          NRXMAX = NRXMAX + 1
          R(2:NRXMAX) = R(1:NRXMAX-1)
          R(1) = 0.D0
          DO NTX = 1, NTXMAX
             F2(2:NRXMAX,NTX) = F2(1:NRXMAX-1,NTX)
             F2(1,NTX) = 0.D0
          ENDDO
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE uf_bin_2d_read


  SUBROUTINE data_error_correct(KFID,RUF,F2,NTXMAX,NRUM,NTUM)
! ----------------------------------------------------------------------
!     WRONG DATA CORRECT
! ----------------------------------------------------------------------
    USE libitp
    IMPLICIT NONE

    CHARACTER(LEN=10),     INTENT(IN) :: KFID
    INTEGER(ikind),        INTENT(IN) :: NTXMAX, NRUM, NTUM
    REAL(rkind), DIMENSION(NRUM)      :: RUF
    REAL(rkind), DIMENSION(NTUM,NRUM) :: F2

    INTEGER(ikind) :: NTX

    IF(KUFDEV == 'jet' .AND. KFID == 'GRHO1') THEN
       IF(KUFDCG=='57987' .OR. KUFDCG=='58159' .OR. KUFDCG=='58323') THEN
          DO NTX=1,NTXMAX
             F2(NTX,1)=FCTR(RUF(2),RUF(3),F2(NTX,2),F2(NTX,3))
          ENDDO
       ENDIF

       IF(KUFDCG == '35156' .OR. KUFDCG == '35171') THEN
          DO NTX=1,NTXMAX
             F2(NTX,2)=FCTR(RUF(3),RUF(4),F2(NTX,3),F2(NTX,4))
             F2(NTX,1)=FCTR(RUF(2),RUF(3),F2(NTX,2),F2(NTX,3))
          ENDDO
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE data_error_correct

! ************************************************************************
! ************************************************************************

  SUBROUTINE uf_head_search(UNIT,MAXVAL1,MAXVAL2,MAXVAL3,IERR)
! -------------------------------------------------------------------------
!   *** search the phrase matching key sentences and get the certain values
!   
!  < input >
!   UNIT  : the unit number for input
!
!  < output >
!   the maximum number of data ...
!   MAXVAL1 : ... corresponding to the phrase '# OF X PTS' (e.g. NRXMAX)
!   MAXVAL2 : ... corresponding to the phrase '# OF Y PTS' (e.g. NTXMAX)
!   MAXVAL3 : ... corresponding to the phrase '# OF PTS'   (e.g. NTXMAX)
!
!   IERR    : error identifier
! -------------------------------------------------------------------------
    USE libchar
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: UNIT
    INTEGER(ikind),INTENT(OUT) :: MAXVAL1, MAXVAL2, MAXVAL3, IERR
    
    CHARACTER(LEN=80) :: KMATCH1,KMATCH2,KMATCH3, KKLINE,KKLINE1,KKLINE2
    INTEGER(ikind)    :: IMATCH,KL,INCL,IST,LINES

    IERR    = 0
    LINES   = 0
    IMATCH  = 0
    KMATCH1 = '# OF X PTS'
    KMATCH2 = '# OF Y PTS'
    KMATCH3 = '# OF PTS'

    DO 
       LINES = LINES + 1
       READ(UNIT,'(A80)',IOSTAT=IST) KKLINE
       IF(IST > 0) THEN
          WRITE(6,*) 'XX uf_head_search: Line read error.'
          WRITE(6,*) 'XX Lines: ',LINES
          IERR = 1
          RETURN
       ELSEIF(IST < 0) THEN ! end of file
          IERR = 0
          RETURN
       END IF

       KL = len_trim(KKLINE)
       IF(KL == 0)THEN
          CYCLE
       ELSE
          INCL = INDEX(KKLINE,';')
          IF(INCL == 0) THEN
             IF(IMATCH >= 1) EXIT
             CYCLE
          ENDIF
       END IF

       ! KSPLIT, KEXTR and KMATCH are defined in TASK/lib libchar.f90
       CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
!       KL = len_trim(KKLINE1)
       CALL KEXTR(KKLINE2)

       IF(KMATCH(KKLINE2,KMATCH1)) THEN
          READ(KKLINE1,*,IOSTAT=IST) MAXVAL1
          IF(IST /= 0)THEN
             WRITE(6,*) 'XX uf_head_search: Read error. IOSTAT= ',IST
             WRITE(6,*) &
       'XX Error in reading the value corresponding to the phrase: ', &
       &                                                   TRIM(KMATCH1)
             IERR = 2
             RETURN
          END IF
          IMATCH = IMATCH + 1

       ELSE IF(KMATCH(KKLINE2,KMATCH2)) THEN
          READ(KKLINE1,*,IOSTAT=IST) MAXVAL2
          IF(IST /= 0)THEN
             WRITE(6,*) 'XX uf_head_search: Read error. IOSTAT= ',IST
             WRITE(6,*) &
       'XX Error in reading the value corresponding to the phrase: ', &
       &                                                   TRIM(KMATCH2)
             IERR = 2
             RETURN
          END IF
          IMATCH = IMATCH + 1

       ELSE IF(KMATCH(KKLINE2,KMATCH3)) THEN
          READ(KKLINE1,*,IOSTAT=IST) MAXVAL3
          IF(IST /= 0)THEN
             WRITE(6,*) 'XX uf_head_search: Read error. IOSTAT= ',IST
             WRITE(6,*) &
       'XX Error in reading the value corresponding to the phrase: ', &
       &                                                   TRIM(KMATCH3)
             IERR = 2
             RETURN
          END IF
          IMATCH = IMATCH + 1
       ENDIF
    END DO

    BACKSPACE UNIT

    RETURN
  END SUBROUTINE uf_head_search


  SUBROUTINE NDINIT(NDPOS)
! ------------------------------------------------------------------------
!    INITIALIZATION FOR READING LINES OF UFILE
! ------------------------------------------------------------------------
    IMPLICIT NONE

    integer(ikind), intent(out) :: NDPOS
    integer(ikind), parameter :: NLENM=80

    NDPOS=NLENM

    RETURN
  END SUBROUTINE NDINIT


  SUBROUTINE NDREAD(NDPOS,IRD,XD,IERR)
! ------------------------------------------------------------------------
!    READING LINES OF UFILE
! ------------------------------------------------------------------------
    IMPLICIT NONE

    integer(ikind), intent(inout) :: NDPOS
    integer(ikind), intent(in)    :: IRD
    integer(ikind), intent(out)   :: IERR
    real(rkind),    intent(out)   :: XD
    integer(ikind), parameter :: NLENM=80
    integer(ikind) :: NR, NT, NDPOS1, NDPOS2, IST
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

    RETURN
  END SUBROUTINE NDREAD

END MODULE ufread
