! trfixed.f90

MODULE trfixed

  !  *** Define fixed profile of density and temperature ***

  USE task_kinds,ONLY: dp

  PRIVATE
                                            ! density profile
  INTEGER:: ntime_nfixed_max                ! number of time points
  INTEGER:: ndata_nfixed_max                ! number of coef data
  REAL(dp),PUBLIC:: rho_min_nfixed,rho_max_nfixed  ! range of fixed profile
  REAL(dp),PUBLIC,ALLOCATABLE:: time_nfixed(:)     ! time points t_i
  REAL(dp),ALLOCATABLE:: coef_nfixed(:,:)   ! coef data for t_i<= t <t_{i+1}
                                            ! temperature profile
  INTEGER:: ntime_tfixed_max                ! number of time points
  INTEGER:: ndata_tfixed_max                ! number of coef data
  REAL(dp),PUBLIC :: rho_min_tfixed,rho_max_tfixed  ! range of fixed profile
  REAL(dp),PUBLIC,ALLOCATABLE:: time_tfixed(:)     ! time points t_i
  REAL(dp),ALLOCATABLE:: coef_tfixed(:,:)   ! coef data for t_i<= t <t_{i+1}

  PUBLIC tr_prof_nfixed ! set fixed density profile
  PUBLIC tr_prof_tfixed ! set fixed temperature profile
  PUBLIC tr_prep_nfixed ! read fixed density pfofile parameters
  PUBLIC tr_prep_tfixed ! read fixed temperature profile parameters

CONTAINS

  ! *** set fixed density profile ***
  
  SUBROUTINE tr_prof_nfixed(rho,time,rn)
  
    USE task_kinds,ONLY: dp
    IMPLICIT NONE    
    REAL(dp),INTENT(IN):: rho,time
    REAL(dp),INTENT(OUT):: rn
    REAL(dp):: tr_func_nfixed
    REAL(dp),ALLOCATABLE:: coef(:)
    REAL(dp):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LT.time_nfixed(1)) THEN
       RETURN
    ELSE IF (time.GE.time_nfixed(ntime_nfixed_max)) THEN
       id=ntime_nfixed_max
    ELSE
       DO ntime=1,ntime_nfixed_max-1
          IF(time.GE.time_nfixed(ntime).AND. &
               time.LT.time_nfixed(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_nfixed_max))
    IF(id.EQ.ntime_nfixed_max) THEN ! after time_nfixed(ntime_nfixed_max)
       DO i=0,ndata_nfixed_max
          coef(i)=coef_nfixed(i,ntime_nfixed_max)
       END DO
    ELSE ! between time_nfixed(id) and time_nfixed(id+1)
       factor=(time-time_nfixed(id)) &
             /(time_nfixed(id+1)-time_nfixed(id))
       DO i=0,ndata_nfixed_max
          coef(i)=(1.D0-factor)*coef_nfixed(i,id) &
                        +factor*coef_nfixed(i,id+1)
       END DO
    END IF

    ! --- set local density profile ---
    
    rn=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    RETURN
  END SUBROUTINE tr_prof_nfixed

  ! *** set fixed temperature profile ***
  
  SUBROUTINE tr_prof_tfixed(rho,time,rt)
  
    USE task_kinds,ONLY: dp
    IMPLICIT NONE    
    REAL(dp),INTENT(IN):: rho,time
    REAL(dp),INTENT(OUT):: rt
    REAL(dp):: tr_func_tfixed
    REAL(dp),ALLOCATABLE:: coef(:)
    REAL(dp):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LE.time_tfixed(1)) THEN
       RETURN
    ELSE IF (time.GE.time_tfixed(ntime_tfixed_max)) THEN
       id=ntime_tfixed_max
    ELSE
       DO ntime=1,ntime_tfixed_max-1
          IF(time.GE.time_tfixed(ntime).AND. &
               time.LE.time_tfixed(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_tfixed_max))
    IF(id.EQ.0) THEN ! before time_nfixed(1)
       DO i=0,ndata_tfixed_max
          coef(i)=coef_tfixed(i,1)
       END DO
    ELSE IF(id.EQ.ntime_tfixed_max) THEN ! after time_nfixed(ntime_nfixed_max)
       DO i=0,ndata_tfixed_max
          coef(i)=coef_tfixed(i,ntime_tfixed_max)
       END DO
    ELSE ! between time_nfixed(id) and time_nfixed(id+1)
       factor=(time-time_tfixed(id)) &
             /(time_tfixed(id+1)-time_tfixed(id))
       DO i=0,ndata_tfixed_max
          coef(i)=(1.D0-factor)*coef_tfixed(i,id) &
                        +factor*coef_tfixed(i,id+1)
       END DO
    END IF

    ! --- set temperature profile ---
    
    rt=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    RETURN
  END SUBROUTINE tr_prof_tfixed

  ! *** read density profile data from file ***

  SUBROUTINE tr_prep_nfixed
    USE trcomm,ONLY: knam_nfixed
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_nfixed,1,0,'fn',ierr)
    READ(NFL,*) ntime_nfixed_max,ndata_nfixed_max,rho_min_nfixed,rho_max_nfixed
    ALLOCATE(time_nfixed(ntime_nfixed_max))
    ALLOCATE(coef_nfixed(0:ndata_nfixed_max,ntime_nfixed_max))
    DO ntime=1,ntime_nfixed_max
       READ(NFL,*) time_nfixed(ntime)
       READ(NFL,*) (coef_nfixed(ndata,ntime),ndata=0,ndata_nfixed_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_nfixed
    
  ! *** read temperature profile data from file ***

  SUBROUTINE tr_prep_tfixed
    USE trcomm,ONLY: knam_tfixed
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_tfixed,1,0,'ft',ierr)
    READ(NFL,*) ntime_tfixed_max,ndata_tfixed_max,rho_min_tfixed,rho_max_tfixed
    ALLOCATE(time_tfixed(ntime_tfixed_max))
    ALLOCATE(coef_tfixed(0:ndata_tfixed_max,ntime_tfixed_max))
    DO ntime=1,ntime_tfixed_max
       READ(NFL,*) time_tfixed(ntime)
       READ(NFL,*) (coef_tfixed(ndata,ntime),ndata=0,ndata_tfixed_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_tfixed
END MODULE trfixed
