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
       WRITE(6,'(A)') '# SELECT : R1-5, P1-6, T1-5, G1-5, ?/HELP X/EXIT'
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
       CASE('P')
          CALL TI_GOUT_P(K2)
       CASE('G')
          CALL TI_GOUT_G(K2)
       CASE('T')
          CALL TI_GOUT_T(K2)
       CASE('?')
          WRITE(6,*) ' R: radial profile'
          WRITE(6,*) '    1:n  2:T  3:u  4:ln n  5:;n T'
          WRITE(6,*) ' P: ADPOST Z profile'
          WRITE(6,*) '    1:Be 2:C  3:Ar 4:Ne 5:Fe 6:W'
          WRITE(6,*) ' T: time evolution'
          WRITE(6,*) '    1:n  2:T  3:u  4:ln n  5:;n T'
          WRITE(6,*) ' G: time evolution of profile'
          WRITE(6,*) '    1:n  2:T  3:u  4:ln n  5:;n T'
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
    NFMAX=NSA2-NSA1+1

    ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
    DO NX=1,NXMAX
       XDATA(NX)=RM(NX)
    END DO
    DO NX=1,NXMAX
       DO NF=1,NFMAX
          NSA=NSA1+NF-1
          SELECT CASE(ID)
          CASE(1)
             FDATA(NX,NF)=RNA(NSA,NX)
          CASE(2)
             FDATA(NX,NF)=RTA(NSA,NX)
          CASE(3)
             FDATA(NX,NF)=RUA(NSA,NX)
          CASE(4)
             FDATA(NX,NF)=LOG10(MAX(RNA(NSA,NX),glog_min))
          CASE(5)
             FDATA(NX,NF)=LOG10(MAX(RTA(NSA,NX),glog_min))
          CASE DEFAULT
             GoTO 1
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
    END SELECT
    CALL PAGEE
    DEALLOCATE(XDATA,FDATA)
    GO TO 1

8000 CONTINUE
    RETURN

9000 WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_r
    
! --- P Graphics: Adpost profile --- 

  SUBROUTINE ti_gout_p(K2)

    USE ticomm
    USE libgrf
    USE ADPOST
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID

    READ(K2,'(I1)',ERR=9000,END=9000) ID
    GO TO 2

1   WRITE(6,'(A)') &
         '## INPUT: ID: 1:Be 2:C 3:Ne 4:Ar 5:Fe 6:W, 0 for end'
    READ(5,*,END=8000,ERR=1) ID

2   CONTINUE
    IF(ID.LE.0) THEN
       GOTO 8000
    ELSE IF(ID.GT.6) THEN
       GOTO 1
    ENDIF

    NXMAX=NRMAX
    NFMAX=1

    ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
    DO NX=1,NXMAX
       XDATA(NX)=RM(NX)
    END DO
    DO NX=1,NXMAX
       DO NF=1,NFMAX
          SELECT CASE(ID)
          CASE(1)
             FDATA(NX,NF)=func_adpost(4,1,RTA(1,NX))
          CASE(2)
             FDATA(NX,NF)=func_adpost(6,1,RTA(1,NX))
          CASE(3)
             FDATA(NX,NF)=func_adpost(10,1,RTA(1,NX))
          CASE(4)
             FDATA(NX,NF)=func_adpost(18,1,RTA(1,NX))
          CASE(5)
             FDATA(NX,NF)=func_adpost(26,1,RTA(1,NX))
          CASE(6)
             FDATA(NX,NF)=func_adpost(74,1,RTA(1,NX))
          CASE DEFAULT
             GO TO 1
          END SELECT
       END DO
    END DO

    CALL PAGES
    SELECT CASE(ID)
    CASE(1)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Be: Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(2)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@C:  Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(3)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Ne: Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(4)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Ar: Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(5)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Fe: Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(6)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@W:  Z vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    END SELECT
    CALL PAGEE
    DEALLOCATE(XDATA,FDATA)
    GOTO 1

8000 CONTINUE
    RETURN

9000 WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_p
    
! --- T Graphics: time evolution --- 

  SUBROUTINE ti_gout_t(K2)

    USE ticomm
    USE tirecord
    USE libgrf
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID,NSA,NSA1X,NSA2X,IDX
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
    CASE(1:5)
       NXMAX=ngt_max
       NFMAX=NSA2-NSA1+1
       IF(NFMAX.EQ.1) NFMAX=3

       WRITE(6,*) 'NXMAX,NFMAX=',NXMAX,NFMAX

       ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
       DO NX=1,NXMAX
          XDATA(NX)=gt(NX)
       END DO
       IDX=MOD(ID-1,3)+1
       SELECT CASE(ID)
       CASE(1:3)
          DO NX=1,NXMAX
             IF(NSA1.EQ.NSA2) THEN
                FDATA(NX,1)=gvta(NX,NSA1,IDX+6)
                FDATA(NX,2)=gvta(NX,NSA1,IDX)
                FDATA(NX,3)=gvta(NX,NSA1,IDX+3)
             ELSE
                DO NSA=NSA1,NSA2
                   NF=NSA-NSA1+1
                   FDATA(NX,NF)=gvta(NX,NSA,IDX+6)
                END DO
             END IF
          END DO
       CASE(4:5)
          DO NX=1,NXMAX
             IF(NSA1.EQ.NSA2) THEN
                FDATA(NX,1)=LOG10(MAX(gvta(NX,NSA1,IDX+6),glog_min))
                FDATA(NX,2)=LOG10(MAX(gvta(NX,NSA1,IDX),glog_min))
                FDATA(NX,3)=LOG10(MAX(gvta(NX,NSA1,IDX+3),glog_min))
             ELSE
                DO NSA=NSA1,NSA2
                   NF=NSA-NSA1+1
                   FDATA(NX,NF)=LOG10(MAX(gvta(NX,NSA,IDX+6),glog_min))
                END DO
             END IF
          END DO
       END SELECT

       CALL PAGES
       SELECT CASE(ID)
       CASE(1)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@n vs t@',0, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       CASE(2)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@T vs t@',0, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       CASE(3)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@u vs t@',0, &
                     XMIN=0.D0,XMAX=gt(ngt_max),FMIN=0.D0)
       CASE(4)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln n vs t@',2, &
                     XMIN=0.D0,XMAX=gt(ngt_max))
       CASE(5)
          CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln T vs t@',2, &
                     XMIN=0.D0,XMAX=gt(ngt_max))
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
             FDATA(NX,NF)=gvrta(NX,NF,NSA,1)
          CASE(2)
             FDATA(NX,NF)=gvrta(NX,NF,NSA,2)
          CASE(3)
             FDATA(NX,NF)=gvrta(NX,NF,NSA,3)
          CASE(4)
             FDATA(NX,NF)=MAX(gvrta(NX,NF,NSA,1),glog_min)
          CASE(5)
             FDATA(NX,NF)=MAX(gvrta(NX,NF,NSA,2),glog_min)
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
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@T vs r@',2, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(3)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@u vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(4)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln n vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(5)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln T vs r@',2, &
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
