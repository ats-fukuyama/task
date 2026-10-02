! plprof_TOTAL

MODULE plprof_TOTAL

  USE plcomm,ONLY: rkind
  INTEGER:: nrmax_profm_TOTAL,nrmax_profg_TOTAL
  REAL(rkind),ALLOCATABLE:: spline_rn_TOTAL(:,:,:),spline_rt_TOTAL(:,:,:)
  REAL(rkind),ALLOCATABLE:: spline_qp_TOTAL(:,:)
  REAL(rkind),ALLOCATABLE:: data_rm(:),data_rg(:)
  REAL(rkind),ALLOCATABLE:: derivm(:),derivg(:)

  PUBLIC pl_load_TOTAL
  PUBLIC pl_read_prof_TOTAL
  PUBLIC pl_read_qp_TOTAL

CONTAINS
  
  ! *** load TOTAL profile data ***
    
  SUBROUTINE pl_load_TOTAL(ierr)

    USE plcomm
    USE libfio
    USE libspl1d
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(rkind),ALLOCATABLE:: data_profm(:,:),data_profg(:,:)
    REAL(rkind),ALLOCATABLE:: data_rn(:,:),data_rt(:,:),data_qp(:)
    INTEGER:: i,j,k,nr,ns
    CHARACTER(LEN=256):: line
    
    CALL FROPEN(21,knam_profm_TOTAL,1,0,'PN',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,A,I8)') &
            'XX file open error in tr_load_TOTAL: knam_profm_TOTAL,ierr=', &
            knam_profm_TOTAL,ierr
       STOP
    END IF

    ! --- count number of data ---
         
    i=0
    READ(21,'(A)') ! skip time data
    READ(21,'(A)') ! skip label line
100 CONTINUE
    READ(21,*,ERR=190,END=200) j
    i=i+1
    GO TO 100
190 CONTINUE
    WRITE(6,*) 'XX file read error in tr_profm_TOTAL: i=',i
    STOP
200 CONTINUE
    nrmax_profm_TOTAL=i+1 ! r=0 is added
    ALLOCATE(data_profm(21,nrmax_profm_TOTAL))
    ALLOCATE(data_rm(nrmax_profm_TOTAL))
    ALLOCATE(data_rn(nrmax_profm_TOTAL,6))
    ALLOCATE(data_rt(nrmax_profm_TOTAL,6))
    ALLOCATE(spline_rn_TOTAL(4,nrmax_profm_TOTAL,6))
    ALLOCATE(spline_rt_TOTAL(4,nrmax_profm_TOTAL,6))
    ALLOCATE(derivm(nrmax_profm_TOTAL))
    nsmax=4
    nszmax=2
    nsnmax=0
    nstmax=nsmax+nszmax+nsnmax
         
    ! --- read electron density profile data ---
         
    REWIND(21)
    READ(21,'(A)') ! skip time data
    READ(21,'(A)') ! skip label line
    DO nr=1,nrmax_profm_TOTAL-1
       READ(21,'(I4,21E12.4)',ERR=390,END=400) &
            j,(data_profm(k,nr),k=1,21)
    END DO
    
    NS_e=1
    NPA(NS_e)= 0  
    PA(NS_e) = AME/AMP
    PZ(NS_e) =-1.D0
         
    NS_D=2
    NPA(NS_D)= 1
    PA(NS_D) = 2.D0
    PZ(NS_D) = 1.D0

    NS_T=3
    NPA(NS_T)= 1
    PA(NS_T) = 3.D0
    PZ(NS_T) = 1.D0

    NS_He4=4
    NPA(NS_He4)= 2
    PA(NS_He4) = 4.D0
    PZ(NS_He4) = 2.D0

    DO nr=2,nrmax_profm_TOTAL
       data_rm(nr)=data_profm(1,nr-1)
       data_rn(nr,1)=data_profm( 6,nr-1)*1.D-20  ! n_e
       data_rn(nr,2)=data_profm( 7,nr-1)*1.D-20  ! n_D
       data_rn(nr,3)=data_profm( 8,nr-1)*1.D-20  ! n_T
       data_rn(nr,4)=data_profm(10,nr-1)*1.D-20  ! n_He4
       data_rt(nr,1)=data_profm( 4,nr-1)  ! T_e
       data_rt(nr,2)=data_profm( 5,nr-1)  ! T_D
       data_rt(nr,3)=data_profm( 5,nr-1)  ! T_T
       data_rt(nr,4)=data_profm( 5,nr-1)  ! T_He4
       data_rn(nr,5)=data_profm(12,nr-1)*1.D-20  ! n_imp1
       data_rn(nr,6)=data_profm(13,nr-1)*1.D-20  ! n_imp2
       data_rt(nr,5)=data_profm( 5,nr-1)  ! T_imp1
       data_rt(nr,6)=data_profm( 5,nr-1)  ! T_imp2
    END DO
    data_rm(1)=0.D0
    !  f(2)=f(1)+f''*X(2)**2
    !  f(3)=f(1)+F''*X(3)**2
    !  f''=(f(2)-f(1))/X(2)**2=(f(3)-f(1))/X(3)**2
    !  (f(2)-f(1))*X(3)**2=(f(3)-f(1))*X(2)**2
    !  (X(2)**2-X(3)**2)*f(1)=(f(3)*X(2)**2-f(2)*X(3)**2)
    !  f(1)=(f(2)*X(3)**2-f(3)*X(2)**2)/(X(3)**2-X(2)**2)
    DO ns=1,6
       data_rn(1,ns)=(data_rn(2,ns)*data_rm(3)**2 &
                     -data_rn(3,ns)*data_rm(2)**2) &
                     /(data_rm(3)**2-data_rm(2)**2)
       data_rt(1,ns)=(data_rt(2,ns)*data_rm(3)**2 &
                     -data_rt(3,ns)*data_rm(2)**2) &
                     /(data_rm(3)**2-data_rm(2)**2)
    END DO
    WRITE(6,*) 'nrmax_profm_TOTAL=',nrmax_profm_TOTAL
    DO ns=1,6
       CALL SPL1D(data_rm,data_rn(:,ns),derivm,spline_rn_TOTAL(:,:,ns), &
            nrmax_profm_TOTAL,1,ierr)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,I6)') 'XX pl_load_TOTAL: SPL1D error: rn: ns=',ns
          STOP
       END IF
       CALL SPL1D(data_rm,data_rt(:,ns),derivm,spline_rt_TOTAL(:,:,ns), &
            nrmax_profm_TOTAL,1,ierr)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,I6)') 'XX pl_load_TOTAL: SPL1D error: rt: ns=',ns
          STOP
       END IF
    END DO
    GOTO 400
390 WRITE(6,*) 'XX profm data error: nr=',nr
    STOP
               
400 CONTINUE

    CALL FROPEN(21,knam_profg_TOTAL,1,0,'PN',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,A,I8)') &
            'XX file open error in tr_prof_TOTAL: knam_profg_TOTAL,ierr=', &
            knam_profg_TOTAL,ierr
       STOP
    END IF

    ! --- count number of data ---
         
    i=0
    READ(21,'(A)') ! skip time data
    READ(21,'(A)') ! skip label line
500 CONTINUE
    READ(21,*,ERR=590,END=600)
    i=i+1
    GO TO 500
590 CONTINUE
    WRITE(6,*) 'XX file read error in data_profg_TOTAL: i=',i
    STOP
600 CONTINUE
    nrmax_profg_TOTAL=i
    WRITE(6,*) 'nrmax_profg_TOTAL=',nrmax_profg_TOTAL
    ALLOCATE(data_profg(27,nrmax_profg_TOTAL))
    ALLOCATE(data_rg(nrmax_profg_TOTAL))
    ALLOCATE(data_qp(nrmax_profg_TOTAL))
    ALLOCATE(spline_qp_TOTAL(4,nrmax_profg_TOTAL))
    ALLOCATE(derivg(nrmax_profg_TOTAL))
         
    ! --- read electron density profile data ---
         
    REWIND(21)
    READ(21,'(A)') line ! skip time data
    READ(21,'(A)') line! skip label line
    DO nr=1,nrmax_profg_TOTAL
       READ(21,'(E11.4,26(1X,E11.4))',ERR=790,END=800) &
            (data_profg(k,nr),k=1,27)
    END DO
    DO nr=1,nrmax_profg_TOTAL
       data_rg(nr)=data_profg(1,nr)
       data_qp(nr)=data_profg(14,nr)  ! q_p
    END DO
    CALL SPL1D(data_rg,data_qp,derivg,spline_qp_TOTAL(:,:), &
               nrmax_profg_TOTAL,0,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I6)') 'XX pl_load_TOTAL: SPL1D error: qp:'
       STOP
    END IF
    GOTO 800
790 WRITE(6,*) 'XX profg data error: nr,k=',nr,k
    STOP
               
800 CONTINUE
    RETURN
  END SUBROUTINE pl_load_TOTAL

  ! *** read TOTAL profile data for rhon ***
    
  SUBROUTINE pl_read_prof_TOTAL(rhon,nsmaxl,pnl,ptl)

    USE plcomm
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    INTEGER,INTENT(IN):: nsmaxl
    REAL(rkind),INTENT(OUT):: pnl(nsmaxl),ptl(nsmaxl)
    INTEGER:: ns,ierr
    
    DO ns=1,nsmaxl
       CALL SPL1DF(rhon,pnl(ns),data_rm,spline_rn_TOTAL(:,:,ns), &
            nrmax_profm_TOTAL,ierr)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,ES12.4,I6,I6)') &
               'XX pl_read_prof_TOTAL: SPL1DF error: rn: rhon,ns,ierr=', &
               rhon,ns,ierr
          STOP
       END IF
       CALL SPL1DF(rhon,ptl(ns),data_rm,spline_rt_TOTAL(:,:,ns), &
            nrmax_profm_TOTAL,ierr)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,ES12.4,I6,I6)') &
               'XX pl_read_prof_TOTAL: SPL1DF error: rt: rhon,ns,ierr=', &
               rhon,ns,ierr
          STOP
       END IF
    END DO
    RETURN
  END SUBROUTINE pl_read_prof_TOTAL
  
  ! *** read TOTAL qp data for rhon ***
    
  SUBROUTINE pl_read_qp_TOTAL(rhon,qpl)

    USE plcomm
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    REAL(rkind),INTENT(OUT):: qpl
    INTEGER:: ierr
    
    CALL SPL1DF(rhon,qpl,data_rg,spline_qp_TOTAL, &
                nrmax_profg_TOTAL,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,ES12.4,I6)') &
            'XX pl_read_qp_TOTAL: SPL1DF error: rhon,ierr:',rhon,ierr
       STOP
    END IF
    RETURN
  END SUBROUTINE pl_read_qp_TOTAL
  
END MODULE plprof_TOTAL
