!     ***********************************************************

!          FILE Output PROFILE DATA

!     ***********************************************************

      SUBROUTINE trfout

      USE TRCOMM,ONLY : GT,GRM,GVT,GVR,GRT,NGT,NCTM,NCRTM,NRMAX,NGR,NCGM
      IMPLICIT NONE
      CHARACTER(LEN=80):: LINE
      CHARACTER(LEN=80),SAVE:: FILENAME='stdout'
      INTEGER,SAVE:: NFL=6
      CHARACTER(LEN=2):: KID
      CHARACTER(LEN=8):: KNCTM,KNCGM,KNCRTM
      INTEGER:: I,ID,IERR

      WRITE(6,'(A,I3,A,A)') 'NFL=',NFL, &
                            '  FILENAME=',FILENAME(1:LEN_TRIM(FILENAME))
1     CALL KKINT(NCTM,KNCTM)
      CALL KKINT(NCGM,KNCGM)
      CALL KKINT(NCRTM,KNCRTM)
      WRITE(6,'(12A)') &
           '# FOUT: GT[0-',KNCTM(1:LEN_TRIM(KNCTM)),'], ', &
                   'RT[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RN[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RG[0-',KNCGM(1:LEN_TRIM(KNCGM)),'], F:filename, ?:help, X:exit'
      READ(5,'(A80)',END=9000,ERR=1) LINE
      KID=LINE(1:2)
      CALL GUCPTL(KID(1:1))
      CALL GUCPTL(KID(2:2))
      IF(KID(1:1).EQ.'X') GOTO 9000
      IF(KID(1:1).EQ.'F') THEN
3        WRITE(6,'(A)') '# output file name?'
         READ(5,*,END=1,ERR=3) FILENAME
         IF(FILENAME(1:6).EQ.'stdout') THEN
            NFL=6
         ELSE
            NFL=21
            CALL FWOPEN(NFL,FILENAME,1,1,'FOUT',IERR)
            IF(IERR.NE.0) GOTO 1
         ENDIF
         GOTO 1
      ENDIF
      IF(KID(1:1).EQ.'?') THEN
         CALL VIEWGTLIST(NCTM)
         CALL VIEWRTLIST(NCRTM)
         GOTO 1
      ENDIF

      READ(LINE(3:),*,END=1,ERR=1) ID
      SELECT CASE(KID)
!     ----- history of global data -----
      CASE('GT')
         IF(ID.EQ.0) then
            DO I=1,NCTM
               CALL TRF1DGT(NFL,GT,GVT,NGT,I)
            ENDDO
         ELSE
            CALL TRF1DGT(NFL,GT,GVT,NGT,ID)
         ENDIF
!     ----- history of profile data -----
      CASE('RT')
         IF(ID.EQ.0) then
            DO I=1,NCTM
               CALL TRF2DRT(NFL,GRM,GT,GRT,NRMAX,NGT,I)
            ENDDO
         ELSE
            CALL TRF2DRT(NFL,GRM,GT,GRT,NRMAX,NGT,ID)
         ENDIF
!     ----- snapshot of profile data -----
      CASE('RN')
         IF(ID.EQ.0) then
            DO I=1,NCTM
               CALL TRF2DRN(NFL,GRM,GRT,NRMAX,NGT,I)
            ENDDO
         ELSE
            CALL TRF2DRN(NFL,GRM,GRT,NRMAX,NGT,ID)
         ENDIF
!     ----- history of profile data -----
      CASE('RG')
         IF(ID.EQ.0) then
            DO I=1,NCGM
               CALL TRF2DRG(NFL,GRM,GT,GVR,NRMAX,NGR,I)
            ENDDO
         ELSE
            CALL TRF2DRG(NFL,GRM,GT,GVR,NRMAX,NGR,ID)
         ENDIF
      CASE('FN')
      END SELECT
      GOTO 1

9000  RETURN
      END SUBROUTINE TRFOUT

!     ===== 1D file output ====

      SUBROUTINE TRF1DGT(NFL,GT,GF,NTMAX,ID)
      USE TRCOMM,ONLY: NCTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NTMAX,ID
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NTMAX,NCTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NT

      CALL GETKGT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(t)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NTMAX
      WRITE(NFL,'(1P2E15.7)') (GT(NT),GF(NT,ID),NT=1,NTMAX)
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I6,A)') ID,':',KSTR,'(',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF1DGT

!     ===== 2D file output ====

      SUBROUTINE TRF2DRT(NFL,GR,GT,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCRTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NRMAX,NTMAX,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTMAX
      WRITE(NFL,'(1P5E15.7)') (GR(NR),NR=1,NRMAX)
      WRITE(NFL,'(1P5E15.7)') (GT(NT),NT=1,NTMAX)
      DO NT=1,NTMAX
         WRITE(NFL,'(1P5E15.7)') (GF(NR,NT,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRT

!     ===== 2D file last data output ====

      SUBROUTINE TRF2DRN(NFL,GR,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCRTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NRMAX,NTMAX,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR

      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,tmax)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NRMAX
      DO NR=1,NRMAX
         WRITE(NFL,'(1P2E15.7)') GR(NR),GF(NR,NTMAX,ID)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A)') &
           ID,':',KSTR,'(',NRMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRN

!     ===== 2D file output ====

      SUBROUTINE TRF2DRG(NFL,GR,GT,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCGM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NRMAX,NTMAX,NCGM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTMAX
      WRITE(NFL,'(1P5E15.7)') (GR(NR),NR=1,NRMAX)
      WRITE(NFL,'(1P5E15.7)') (GT(NT),NT=1,NTMAX)
      DO NT=1,NTMAX
         WRITE(NFL,'(1P5E15.7)') (GF(NR,NT,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRG

