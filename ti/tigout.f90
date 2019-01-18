! MODULE tigout

MODULE tigout

  PRIVATE
  PUBLIC ti_gout

CONTAINS

  SUBROUTINE ti_gout

    USE ticomm
    IMPLICIT NONE
    CHARACTER(LEN=2):: KIG
    CHARACTER(LEN=1):: K1,K2
    INTEGER:: IST


    DO
       WRITE(6,'(A)') '# SELECT : R1-8, R9, T1-3, G1-3, X/EXIT'
       READ(5,'(A2)',IOSTAT=IST) KIG
       IF(IST.GT.0) CYCLE
       IF(IST.LT.0) EXIT
       K1=KIG(1:1)
       K2=KIG(2:2)
       CALL GUCPTL(K1)
       CALL GUCPTL(K2)

       SELECT CASE(K1)
       CASE('R')
          CALL TI_GOUT_R(K2)
       CASE('G')
          CALL TI_GOUT_G(K2)
       CASE('T')
          CALL TI_GOUT_T(K2)
       CASE('X')
          EXIT
       END SELECT
    END DO
    RETURN
  END SUBROUTINE ti_gout

! --- R Graphics: present radial profile --- 

  SUBROUTINE ti_gout_r(K2)

    USE ticomm
    USE libgrf
    USE ADPOST
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID,NSA,NSA1X,NSA2X
    INTEGER,SAVE:: NSA1=1
    INTEGER,SAVE:: NSA2=0

    READ(K2,'(I1)',ERR=9000,END=9000) ID

1   WRITE(6,'(A,2I5)') &
         '## INPUT: NSA1,NSA2 (NSA1=0 for end,=NSA2 for single) ',NSA1,NSA2
    NSA1X=NSA1
    NSA2X=NSA2
    READ(5,*,END=8000,ERR=1) NSA1X,NSA2X
    IF(NSA1X.LE.0) THEN
       GOTO 8000
    ELSE
       NSA1=NSA1X
    ENDIF
    IF(NSA2X.LE.0) THEN
       NSA2=NSA1
    ELSE
       NSA2=NSA2X
    ENDIF
    IF(NSA1.GT.nsa_max) NSA1=nsa_max
    IF(NSA2.GT.nsa_max) NSA2=nsa_max

    NXMAX=NRMAX
    SELECT CASE(ID)
    CASE(1:5)
       NFMAX=NSA2-NSA1+1
    CASE DEFAULT
       NFMAX=1
    END SELECT

    ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
    DO NX=1,NXMAX
       XDATA(NX)=RM(NX)
    END DO
    DO NX=1,NXMAX
       DO NF=1,NFMAX
          NSA=NSA1+NF-1
          SELECT CASE(ID)
          CASE(1)
             FDATA(NX,NF)=MAX(RNA(NSA,NX),1.D-6)
          CASE(2)
             FDATA(NX,NF)=MAX(RTA(NSA,NX),1.D-6)
          CASE(3)
             FDATA(NX,NF)=RUA(NSA,NX)
          CASE(4)
             FDATA(NX,NF)=LOG10(MAX(RNA(NSA,NX),1.D-6))
          CASE(5)
             FDATA(NX,NF)=LOG10(MAX(RTA(NSA,NX),1.D-6))
          CASE(6)
             FDATA(NX,NF)=func_adpost(6,1,RTA(1,NX))
          CASE(7)
             FDATA(NX,NF)=func_adpost(18,1,RTA(1,NX))
          CASE(8)
             FDATA(NX,NF)=func_adpost(26,1,RTA(1,NX))
          CASE(9)
             FDATA(NX,NF)=func_adpost(74,1,RTA(1,NX))
          CASE DEFAULT
             FDATA(NX,NF)=0.D0
          END SELECT
       END DO
    END DO

    CALL PAGES
    SELECT CASE(ID)
    CASE(1)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@N vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(2)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@T vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(3)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@U vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(4)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln N vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(5)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln T vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(6)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z C vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(7)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z Ar vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(8)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z Fe vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(9)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z W vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    END SELECT
    CALL PAGEE
    DEALLOCATE(XDATA,FDATA)
8000 CONTINUE
    RETURN

9000 WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_r
    
! --- T Graphics: time evolution --- 

  SUBROUTINE ti_gout_t(K2)

    USE ticomm
    USE tirecord
    USE libgrf
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID,NSA,NSA1X,NSA2X
    INTEGER,SAVE:: NSA1=1
    INTEGER,SAVE:: NSA2=0

    READ(K2,'(I1)',ERR=9000,END=9000) ID

1   WRITE(6,'(A,2I5)') &
         '## INPUT: NSA1,NSA2 (NSA1=0 for end,NSA2= for single) ',NSA1,NSA2
    NSA1X=NSA1
    NSA2X=NSA2
    READ(5,*,END=8000,ERR=1) NSA1X,NSA2X
    IF(NSA1X.LE.0) THEN
       GOTO 8000
    ELSE
       NSA1=NSA1X
    ENDIF
    IF(NSA2X.LE.0) THEN
       NSA2=NSA1
    ELSE
       NSA2=NSA2X
    ENDIF
    IF(NSA1.GT.nsa_max) NSA1=nsa_max
    IF(NSA2.GT.nsa_max) NSA2=nsa_max

    SELECT CASE(ID)
    CASE(1:3)
       NXMAX=ngt_max
       NFMAX=NSA2-NSA1+1
       IF(NFMAX.EQ.1) NFMAX=3

       WRITE(6,*) 'NXMAX,NFMAX=',NXMAX,NFMAX

       ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
       DO NX=1,NXMAX
          XDATA(NX)=gt(NX)
       END DO
       DO NX=1,NXMAX
          IF(NSA1.EQ.NSA2) THEN
             FDATA(NX,1)=gvta(NX,NSA1,ID+6)
             FDATA(NX,2)=gvta(NX,NSA1,ID)
             FDATA(NX,3)=gvta(NX,NSA1,ID+3)
          ELSE
             DO NSA=NSA1,NSA2
                NF=NSA-NSA1+1
                FDATA(NX,NF)=gvta(NX,NSA,ID)
             END DO
          END IF
       END DO

       CALL PAGES
       SELECT CASE(ID)
       CASE(1)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@N(0) vs t@',0, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       CASE(2)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@N(a) vs t@',2, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       CASE(3)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@U(0) vs t@',0, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       END SELECT
       DEALLOCATE(XDATA,FDATA)
    END SELECT
    CALL PAGEE
    GO TO 1

8000 CONTINUE
    RETURN

9000 CONTINUE
    WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_t
    
! --- G Graphics: time evolution of radial profile --- 

  SUBROUTINE ti_gout_g(K2)

    USE ticomm
    USE tirecord
    USE libgrf
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID,NSAX
    INTEGER,SAVE:: NSA=1

    READ(K2,'(I1)',ERR=9000,END=9000) ID

1   WRITE(6,'(A,2I5)') &
         '## INPUT: NSA (NSA=0 for end) ',NSA
    NSAX=NSA
    READ(5,*,END=8000,ERR=1) NSAX
    IF(NSAX.LE.0) GO TO 8000
    NSA=NSAX
    IF(NSA.GT.nsa_max) NSA=nsa_max

    NXMAX=NRMAX
    NFMAX=ngr_max
    ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
    DO NX=1,NXMAX
       XDATA(NX)=RM(NX)
    END DO
    DO NX=1,NXMAX
       DO NF=1,NFMAX
          SELECT CASE(ID)
          CASE(1)
             FDATA(NX,NF)=MAX(gvrta(NX,NF,NSA,1),1.D-6)
          CASE(2)
             FDATA(NX,NF)=MAX(gvrta(NX,NF,NSA,2),1.D-6)
          CASE(3)
             FDATA(NX,NF)=MAX(gvrta(NX,NF,NSA,3),1.D-6)
          CASE DEFAULT
             FDATA(NX,NF)=0.D0
          END SELECT
       END DO
    END DO

    CALL PAGES
    SELECT CASE(ID)
    CASE(1)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@N vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(2)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@U vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(3)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@T vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    END SELECT
    CALL PAGEE
    DEALLOCATE(XDATA,FDATA)
    GO TO 1

8000 CONTINUE
    RETURN

9000 CONTINUE
    WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_g
    
END MODULE tigout
