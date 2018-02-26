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
       WRITE(6,'(A)') '# SELECT : R0-R5, R9, T1, X/EXIT'
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
       CASE('T')
!          CALL TI_GOUT_T(K2)
       CASE('X')
          EXIT
       END SELECT
    END DO
    RETURN
  END SUBROUTINE ti_gout

  SUBROUTINE ti_gout_r(K2)

    USE ticomm
    USE libgrf
    USE ADPOST
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA
    INTEGER:: NXMAX,NFMAX,NX,NF,ID

    READ(K2,'(I1)',ERR=9000,END=9000) ID

    NXMAX=NRMAX
    NFMAX=NSAMAX
    ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,NFMAX))
    DO NX=1,NXMAX
       XDATA(NX)=RM(NX)
    END DO
    DO NX=1,NXMAX
       DO NF=1,NFMAX
          SELECT CASE(ID)
          CASE(0)
             FDATA(NX,NF)=MAX(RNA(NF,NX),1.D-6)
          CASE(1)
             FDATA(NX,NF)=LOG10(MAX(RNA(NF,NX),1.D-6))
          CASE(2)
             FDATA(NX,NF)=RUA(NF,NX)
          CASE(3)
             FDATA(NX,NF)=MIN(MAX(RUA(NF,NX),-1.D6),1.D6)
          CASE(4)
             FDATA(NX,NF)=MAX(RTA(NF,NX),1.D-6)
          CASE(5)
             FDATA(NX,NF)=LOG10(MAX(RTA(NF,NX),1.D-6))
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
    CASE(0)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@N vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(1)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln N vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(2)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@U vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(3)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@U vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(4)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@T vs r@',0, &
                  XMIN=0.D0,XMAX=RA,FMIN=0.D0)
    CASE(5)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@ln T vs r@',2, &
                  XMIN=0.D0,XMAX=RA)
    CASE(7)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z Ar vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(8)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z FE vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    CASE(9)
       CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,NFMAX,'@Z W vs r@',0, &
                  XMIN=0.D0,XMAX=RA)
    END SELECT
    CALL PAGEE
    DEALLOCATE(XDATA,FDATA)
    RETURN

9000 WRITE(6,'(A,A1)') 'XX ti_gout_r: K1 error',K2
    RETURN
  END SUBROUTINE ti_gout_r
    
END MODULE tigout
