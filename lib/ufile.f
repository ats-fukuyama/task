C
C     ***********************************************************
C
C           LOCATING UFILE DIRECTORY
C
C     ***********************************************************
C
      SUBROUTINE UFILE_INTERFACE(KDIRX,KUFDEV,KUFDCG,IKNDEV,IKNDCG,MODE)
C
      INTEGER MODE
      CHARACTER KDIRX*80
      CHARACTER KUFDEV*80,KUFDCG*80
      COMMON /UFMODE/ MODEL
      LOGICAL DIR
C
C     MODE is the parameter which determines how to handle exp. files.
C     MODE = 0 : Binary files are loaded if available, or ASCII files
C                are loaded and aftermath binary files are created.
C            1 : Only binary files are loaded.
C         else : Only ASCII files are loaded and binary files are NOT
C                created.
C
      MODEL = MODE
C
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
      KDIRX='../../../profiledb/profile_data/'//KUFDEV(1:IKNDEV)//'/'
     &                          //KUFDCG(1:IKNDCG)//'/in/'
      INQUIRE(FILE=KDIRX,EXIST=DIR,ERR=9000)
      IF(DIR.NEQV..TRUE.) THEN
         WRITE(6,'(A25,A34,A17)') 
     &        '## DESIGNATED DIRECTORY( ',KDIRX,' ) DOES NOT EXIST!'
         STOP
      ENDIF
C
 9000 RETURN
      END
C
C   **************************************************
C   **    UFILE read for TR (Time Evolution UFILE)  **
C   **************************************************
C
C     input:
C
C     KDIRX    : Directory
C     KUFDEV   : Device name in UFILEs
C     KUFDCG   : Shot number in UFILEs
C     KFID     : Variable name in UFILEs
C
C     output:
C
C     RUF(NRMU)      : Equally Spaced Normalized Radial Data
C     TMU(NTUM)      : Total Time Data (The Number of DT)
C     F1(NTUM)       : Functional Values
C     F2(NTUM,NRMU)  : Functional Values
C     NRFMAX         : Maximum Number of the Radial Mesh
C     NTXMAX         : Maximum Number of the Time Mesh
C     NRMU           : Array size (radial)
C     NTUM           : Array size (time)
C     MDCHK          : Loop Check Value
C     IERR           : Error Indicator
C
C   ***************************************************************
C
      SUBROUTINE UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,
     &                       TT,F1,NTXMAX,NTUM,MDCHK,IERR)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /UFMODE/ MODEL
      DIMENSION TT(NTUM),F1(NTUM)
      CHARACTER KDIRX*80
      CHARACTER KDIRR1*80
      CHARACTER KFILE*80,KFID*10
      CHARACTER KUFDEV*80,KUFDCG*80
      LOGICAL LEX
C
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR1=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'1d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR1,KL1)
      KFILE=KDIRR1(1:KL1)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000) 
      IF(LEX) THEN
         CALL TRXR1D(KDIRR1,KFID,TT,F1,NTUM,NTXMAX,MODEL)
         MDCHK=1
         IERR=0
      ELSE
         DO NTA=1,NTUM
            F1(NTA)=0.D0
         ENDDO
         IERR=1
      ENDIF
C
 9000 RETURN
      END
C
C     *****
C
      SUBROUTINE UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,F2,
     &                        NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /UFMODE/ MODEL
      DIMENSION RUF(NRMU),TMU(NTUM),F2I(NRMU,NTUM),F2(NTUM,NRMU)
      DIMENSION F2CTR(NTUM),F2EDG(NTUM)
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*10
      CHARACTER KUFDEV*80,KUFDCG*80
      LOGICAL LEX
C
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'2d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         NRFMAX=NRMU
         CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2I,NRMU,NTUM,NRFMAX,NTXMAX, 
     &               MODEL)
         IF(NRFMAX.EQ.1) THEN
            WRITE(6,*) '## ',KFID,
     &           'HAS NOT BEEN READ DUE TO A SINGLE RADIAL POINT.'
            DO NTA=1,NTUM
               DO NRF=1,NRMU
                  F2(NTA,NRF)=0.D0
               ENDDO
            ENDDO
            IERR=1
            GOTO 9000
         ENDIF
         MDCHK=1
         IERR=0
         DO NTX=1,NTXMAX
            DO NRF=1,NRFMAX
               F2(NTX,NRF)=F2I(NRF,NTX)
            ENDDO
         ENDDO
      ELSE
         DO NTA=1,NTUM
            DO NRF=1,NRMU
               F2(NTA,NRF)=0.D0
            ENDDO
         ENDDO
         IERR=1
         GOTO 9000
      ENDIF
C
C     *****
C       For interpolation, center and/or edge value is added
C       if necessary. The variables with center value having to
C       be zero such as BPOL, RMINOR, SURF and VOLUME were beforehand 
C       manupulated and thus in this section MD=4 will be selected.
C     *****
C
      MD=0
      IF(RUF(1).NE.0.D0.AND.RUF(NRFMAX).NE.1.D0) THEN
         DO NTX=1,NTXMAX
            F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
            F2EDG(NTX)=AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1),
     &                          F2(NTX,NRFMAX-2),RUF(NRFMAX),
     &                          RUF(NRFMAX-1),RUF(NRFMAX-2))
         ENDDO
         MD=1
      ELSEIF(RUF(1).NE.0.D0.AND.RUF(NRFMAX).EQ.1.D0) THEN
         DO NTX=1,NTXMAX
            F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(NTX,1),F2(NTX,2))
         ENDDO
         MD=2
      ELSEIF(RUF(1).EQ.0.D0.AND.RUF(NRFMAX).NE.1.D0) THEN
         DO NTX=1,NTXMAX
            F2EDG(NTX)=AITKEN2P(1.D0,F2(NTX,NRFMAX),F2(NTX,NRFMAX-1),
     &                          F2(NTX,NRFMAX-2),RUF(NRFMAX),
     &                          RUF(NRFMAX-1),RUF(NRFMAX-2))
         ENDDO
         MD=3
      ELSEIF(RUF(1).EQ.0.D0.AND.RUF(NRFMAX).EQ.1.D0) THEN
         MD=4
      ELSE
         STOP 'XX UFILE: TRXR2D: ERROR'
      ENDIF
C
      CALL DATA_ERROR_CORRECT(KUFDEV,KUFDCG,KFID,RUF,F2,
     &                        NTXMAX,NRMU,NTUM)
C
      IF(MD.EQ.1) THEN
         NRFMAX=NRFMAX+2
         DO NRF=NRFMAX-1,2,-1
            RUF(NRF)=RUF(NRF-1)
         ENDDO
         RUF(1)=0.D0
         RUF(NRFMAX)=1.D0
         DO NTX=1,NTXMAX
            DO NRF=NRFMAX-1,2,-1
               F2(NTX,NRF)=F2(NTX,NRF-1)
            ENDDO
            F2(NTX,1)=F2CTR(NTX)
            F2(NTX,NRFMAX)=F2EDG(NTX)
         ENDDO
      ELSEIF(MD.EQ.2) THEN
         NRFMAX=NRFMAX+1
         DO NRF=NRFMAX,2,-1
            RUF(NRF)=RUF(NRF-1)
         ENDDO
         RUF(1)=0.D0
         DO NTX=1,NTXMAX
            DO NRF=NRFMAX,2,-1
               F2(NTX,NRF)=F2(NTX,NRF-1)
            ENDDO
            F2(NTX,1)=F2CTR(NTX)
         ENDDO
      ELSEIF(MD.EQ.3) THEN
         NRFMAX=NRFMAX+1
         RUF(NRFMAX)=1.D0
         DO NTX=1,NTXMAX
            F2(NTX,NRFMAX)=F2EDG(NTX)
         ENDDO
      ENDIF
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           WRONG DATA CORRECT
C
C     ***********************************************************
C
      SUBROUTINE DATA_ERROR_CORRECT(KUFDEV,KUFDCG,KFID,RUF,F2,
     &                              NTXMAX,NRMU,NTUM)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RUF(NRMU),F2(NTUM,NRMU)
      CHARACTER KUFDEV*80,KUFDCG*80
      CHARACTER KFID*10
C
      IF(KUFDEV.EQ.'jet') THEN
         IF(KFID.EQ.'GRHO1') THEN
            IF(KUFDCG.EQ.'57987'.OR.KUFDCG.EQ.'58159'.OR.
     &         KUFDCG.EQ.'58323') THEN
               DO NTX=1,NTXMAX
                  F2(NTX,1)=FCTR(RUF(2),RUF(3),F2(NTX,2),F2(NTX,3))
               ENDDO
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,F2,
     &                         NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /UFMODE/ MODEL
      DIMENSION RUF(NRMU),TMU(NTUM),F2I(NRMU,NTUM),F2(NTUM,NRMU)
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*10
      CHARACTER KUFDEV*80,KUFDCG*80
      LOGICAL LEX
C
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'2d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         NRFMAX=NRMU           ! equal to NRMU
         CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2I,NRMU,NTUM,NRFMAX,NTXMAX,
     &               MODEL)
         MDCHK=1
         IERR=0
         DO NTX=1,NTXMAX
            DO NRF=1,NRFMAX
               F2(NTX,NRF)=F2I(NRF,NTX)
            ENDDO
         ENDDO
      ELSE
         DO NTX=1,NTUM
            DO NRF=1,NRMU
               F2(NTX,NRF)=0.D0
            ENDDO
         ENDDO
         IERR=1
         GOTO 9000
      ENDIF
C
      IF(RUF(NRFMAX).GT.1.D0) THEN
         DO NRF=1,NRFMAX
            IF(RUF(NRF).GT.1.D0) THEN
               NRFMAX=NRF-1
               GOTO 1000
            ENDIF
         ENDDO
      ENDIF
 1000 CONTINUE
      IF(KFID.EQ.'NEXP'.OR.KFID.EQ.'NEXPEB') THEN
         DO NTX=1,NTXMAX
            DO NRF=1,NRMU
               F2(NTX,NRF)=F2(NTX,NRF)*1.D-20
            ENDDO
         ENDDO
      ELSE
         DO NTX=1,NTXMAX
            DO NRF=1,NRMU
               F2(NTX,NRF)=F2(NTX,NRF)*1.D-3
            ENDDO
         ENDDO
      ENDIF
C
 9000 RETURN
      END
C
C     ***** INITIALIZATION FOR READING UFILE *****
C
      SUBROUTINE NDINIT(NDPOS)
C
      PARAMETER (NLENM=80)
C
      NDPOS=NLENM
      RETURN
      END
C
C     ***** READING UFILE *****
C
      SUBROUTINE NDREAD(NDPOS,IRD,XD,IERR,ID)
C
      PARAMETER (NLENM=80)
      REAL*8 XD
      CHARACTER LINE*(NLENM),NUMD*40
C
   10 IF(NDPOS.EQ.NLENM) THEN
         READ(IRD,'(A80)',ERR=8000,END=9000) LINE
         NDPOS=0
      ENDIF
C
   20 CONTINUE
         NDPOS=NDPOS+1
      IF(LINE(NDPOS:NDPOS).EQ.' '.AND.NDPOS.LT.NLENM) GOTO 20
C
      IF(NDPOS.LT.NLENM) THEN
         NDPOS1=NDPOS
   30    CONTINUE
            NDPOS=NDPOS+1
         IF(LINE(NDPOS:NDPOS).NE.'E'.AND.
     &      LINE(NDPOS:NDPOS).NE.'e'.AND.
     &      NDPOS.LT.NLENM) GOTO 30
C
         IF(NDPOS.LT.NLENM) THEN
            NDPOS=NDPOS+1
   40       CONTINUE
               NDPOS=NDPOS+1
            IF(LINE(NDPOS:NDPOS).NE.' '.AND.
     &         LINE(NDPOS:NDPOS).NE.'-'.AND.
     &         NDPOS.LT.NLENM) GOTO 40
C
            IF(NDPOS.NE.NLENM) NDPOS=NDPOS-1
            NDPOS2=NDPOS
            NUMD=LINE(NDPOS1:NDPOS2)
            READ(NUMD,*,ERR=100) XD
            IF(ID.EQ.1) WRITE(6,*) NDPOS1,NDPOS2,XD
            IERR=0
            RETURN
         ENDIF
      ENDIF
      GOTO 10
C
  100 NDPOS=NLENM
      GOTO 10
C
 8000 IERR=1
      RETURN
 9000 IERR=2
      RETURN
      END
C
C     ***** READ 1D VARIABLE *****
C
      SUBROUTINE TRXR1D(KDIRR1,KFID,T,F1,NTM,NTXMAX,MODEL)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(NTM),F1(NTM)
      CHARACTER KDIRR1*80
      CHARACTER KFID*10,KFIDB*15,KFILE*80,KFILEB*80
      CHARACTER KMATCH1*80,KMATCH2*80,KMATCH3*80
      CHARACTER KKLINE*80,KKLINE1*80,KKLINE2*80
      LOGICAL KMATCH,LEX
C
      MODE=MODEL
      IMATCH=0
      KMATCH1='# OF X PTS'
      KMATCH2='# OF Y PTS'
      KMATCH3='# OF PTS'
C
      CALL KTRIM(KDIRR1,KL1)
      KFILE=KDIRR1(1:KL1)//KFID
      CALL KTRIM(KFILE,KL2)
      KFILEB=KFILE(1:KL2)//'.bin'
C
      INQUIRE(FILE=KFILEB,EXIST=LEX,ERR=9000)
      IF(LEX.AND.MODE.EQ.0) THEN
         MODE=1
      ENDIF
C
      IF(MODE.EQ.1) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='OLD',ERR=8)
         READ(15) NTXMAX
         READ(15) (T(NTX),NTX=1,NTXMAX)
         READ(15) (F1(NTX),NTX=1,NTXMAX)
         CLOSE(15)
         KFIDB=KFID//'.bin'
         WRITE(6,"(' ',3A,I3,A)") ' - READ FILE :',KFIDB,'(',NTXMAX,')'
         GOTO 9000
      ENDIF
C
    8 OPEN(15,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',
     &     STATUS='OLD')
      IF(IST.GT.0) THEN
         WRITE(6,*) 'XX NOT FOUND :',KFILE(1:45)
         NTXMAX=0
         WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
         RETURN
      ENDIF
C
  100 CONTINUE
         READ(15,'(A80)',ERR=8000,END=9000) KKLINE
         CALL KTRIM(KKLINE,KL)
         IF(KL.EQ.0) GOTO 100
         INCL=KKINDEX(KKLINE,';')
         IF(INCL.EQ.0) THEN
            IF(IMATCH.GE.1) GOTO 1000
            GOTO 100
         ENDIF
         CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
         CALL KTRIM(KKLINE1,KL)
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
      GOTO 100
C
 1000 CONTINUE
C
      BACKSPACE 15
C
      CALL NDINIT(NDPOS)
      DO NTX=1,NTXMAX
         CALL NDREAD(NDPOS,15,T(NTX),IERR,0)
         IF(IERR.NE.0) GOTO 8000
      ENDDO
      DO NTX=1,NTXMAX
         CALL NDREAD(NDPOS,15,F1(NTX),IERR,0)
         IF(IERR.NE.0) GOTO 8000
      ENDDO
C
      CLOSE(15)
      WRITE(6,"(' ',3A,I3,A)") ' - READ ASCII FILE :',
     &                         KFID,'(',NTXMAX,')'
C
      IF(MODE.EQ.0) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='NEW')
         IF(IST.GT.0) THEN
            WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
            RETURN
         ENDIF     
         WRITE(15) NTXMAX
         WRITE(15) (T(NT),NT=1,NTXMAX)
         WRITE(15) (F1(NT),NT=1,NTXMAX)
         CLOSE(15)
         WRITE(6,*) ' - WRITE FILE:',KFILEB(1:60)
      ENDIF

 9000 RETURN
C
 8000 WRITE(6,*) 'XX FILE READ ERROR'
      RETURN
      END
C
C     ***** READ 2D VARIABLE *****
C
      SUBROUTINE TRXR2D(KDIRR2,KFID,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,MODEL)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(NTM),R(NRM),F2(NRM,NTM)
      CHARACTER KDIRR2*80
      CHARACTER KFID*10,KFIDB*15,KFILE*80,KFILEB*80
      CHARACTER KMATCH1*80,KMATCH2*80,KMATCH3*80
      CHARACTER KKLINE*80,KKLINE1*80,KKLINE2*80
      LOGICAL KMATCH,LEX
C
      MODE=MODEL
      IMATCH=0
      KMATCH1='# OF X PTS'
      KMATCH2='# OF Y PTS'
      KMATCH3='# OF PTS'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
      CALL KTRIM(KFILE,KL2)
      KFILEB=KFILE(1:KL2)//'.bin'
C
      INQUIRE(FILE=KFILEB,EXIST=LEX,ERR=9000)
      IF(LEX.AND.MODE.EQ.0) THEN
         MODE=1
      ENDIF
C
      IF(MODE.EQ.1) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='OLD',ERR=8)
         READ(15) NRXMAX,NTXMAX
         READ(15) (R(NRX),NRX=1,NRXMAX)
         READ(15) (T(NTX),NTX=1,NTXMAX)
         READ(15) ((F2(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)
         CLOSE(15)
         KFIDB=KFID//'.bin'
         WRITE(6,"(' ',3A,2(I3,A))") ' - READ FILE :',KFIDB,
     &                '(',NRXMAX,',',NTXMAX,')'
         GOTO 9000
      ENDIF
C
    8 OPEN(15,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',
     &     STATUS='OLD')
      IF(IST.GT.0) THEN
         WRITE(6,*) 'XX NOT FOUND :',KFILE(1:55)
         NTXMAX=0
         WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
         RETURN
      ENDIF
C
  100 CONTINUE
         READ(15,'(A80)',ERR=8000,END=9000) KKLINE
         CALL KTRIM(KKLINE,KL)
         IF(KL.EQ.0) GOTO 100
         INCL=KKINDEX(KKLINE,';')
         IF(INCL.EQ.0) THEN
            IF(IMATCH.GE.2) GOTO 1000
            GOTO 100
         ENDIF
         CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
         CALL KTRIM(KKLINE1,KL)
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
      GOTO 100
C
 1000 CONTINUE
C
      BACKSPACE 15
C
      CALL NDINIT(NDPOS)
      DO NRX=1,NRXMAX
         IF(KFID.EQ.'NE'.OR.KFID.EQ.'TE'.OR.KFID.EQ.'TI'.OR.
     &        KFID.EQ.'NFAST'.OR.KFID.EQ.'ZEFFR') THEN
            CALL NDREAD(NDPOS,15,R(NRX),IERR,1)
         ELSE
            CALL NDREAD(NDPOS,15,R(NRX),IERR,0)
         ENDIF
         IF(IERR.NE.0) GOTO 8000
      ENDDO
      DO NTX=1,NTXMAX
         CALL NDREAD(NDPOS,15,T(NTX),IERR,0)
!         IF(KFID.EQ.'NE') write(6,*) NTX,T(NTX)
         IF(IERR.NE.0) GOTO 8000
      ENDDO
      DO NTX=1,NTXMAX
      DO NRX=1,NRXMAX
         CALL NDREAD(NDPOS,15,F2(NRX,NTX),IERR,0)
!         IF(KFID.EQ.'NE'.AND.NTX.LE.3) write(6,*) NTX,NRX,F2(NRX,NTX)
         IF(IERR.NE.0) GOTO 8000
      ENDDO
      ENDDO
C
      CLOSE(15)
      WRITE(6,"(' ',3A,2(I3,A))") ' - READ ASCII FILE :',KFID,
     &             '(',NRXMAX,',',NTXMAX,')'
C
      IF(MODE.EQ.0) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='NEW')
         IF(IST.GT.0) THEN
            WRITE(6,*) 'XX FILE OPEN ERROR: STATUS = ',IST
            RETURN
         ENDIF
         write(6,*) 
         WRITE(15) NRXMAX,NTXMAX
         WRITE(15) (R(NRX),NRX=1,NRXMAX)
         WRITE(15) (T(NTX),NTX=1,NTXMAX)
         WRITE(15) ((F2(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)
         WRITE(6,*) ' - WRITE FILE:',KFILEB(1:60)
         CLOSE(15)
      ENDIF
C
 9000 IF(KFID(1:6).EQ.'VOLUME'.OR.
     &   KFID(1:4).EQ.'BPOL'.OR.
     &   KFID(1:4).EQ.'SURF'.OR.
     &   KFID(1:6).EQ.'RMINOR') THEN
         IF(R(1).NE.0.D0) THEN
            NRXMAX=NRXMAX+1
            DO NRX=NRXMAX,2,-1
               R(NRX)=R(NRX-1)
            ENDDO
            R(1)=0.D0
            DO NTX=1,NTXMAX
               DO NRX=NRXMAX,2,-1
                  F2(NRX,NTX)=F2(NRX-1,NTX)
               ENDDO
               F2(1,NTX)=0.D0
            ENDDO
         ENDIF
      ENDIF
      RETURN
C
 8000 WRITE(6,*) 'XX FILE READ ERROR'
      RETURN
      END
