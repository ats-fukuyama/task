! w1file.f90
  
MODULE w1file
  PRIVATE
  PUBLIC w1_save
  PUBLIC w1_load

CONTAINS

  ! *** save data ***
  
  SUBROUTINE w1_save(ierr)
    USE w1comm
    USE w1parm,ONLY: w1_save_parm
    USE libfio,ONLY: FWOPEN
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    LOGICAL:: l_opened
    CHARACTER(LEN=90):: filename_namelist,filename_data
    
    filename_namelist=TRIM(job_id)//'.namelist'
    filename_data=TRIM(job_id)//'.data'
    
    INQUIRE(nfile_data,OPENED=l_opened)
    IF(l_opened) CLOSE(nfile_data)
    CALL FWOPEN(nfile_data,filename_namelist,1,0,'w1',ierr)
    IF(ierr.EQ.0) THEN
       CALL w1_save_parm(nfile_data)
       CLOSE(nfile_data)
       CALL FWOPEN(nfile_data,filename_data,0,0,'w1',ierr)
       IF(ierr.EQ.0) THEN
          WRITE(nfile_data) NXMAX,NZMAX,NSMAX
          WRITE(nfile_data) XA,XAM,ZA
          WRITE(nfile_data) CE2DA
          WRITE(nfile_data) PABS2D
          WRITE(nfile_data) PROFB,PROFPN,PROFPU,PROFTR,PROFTP
          CLOSE(nfile_data)
       ELSE
          WRITE(6,*) 'XX w1_save: data file open error: ierr=',ierr
       END IF
    ELSE
       WRITE(6,*) 'XX w1_save: namelist file open error: ierr=',ierr
    END IF
  END SUBROUTINE w1_save

! load data file

  SUBROUTINE w1_load(ierr)

    USE w1comm
    USE w1parm,ONLY: w1_load_parm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    LOGICAL:: l_opened
    CHARACTER(LEN=90):: filename_namelist,filename_data
    
    filename_namelist=TRIM(job_id)//'.namelist'
    filename_data=TRIM(job_id)//'.data'

    ! --- open namelist file  ---

    ierr=0

    INQUIRE(nfile_data,OPENED=l_opened)
    IF(l_opened) CLOSE(nfile_data)
       
    CALL FROPEN(nfile_data,filename_namelist,1,0,'w1',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX w1_load: namelist file open error: ierr=',ierr
       RETURN
    END IF
    CALL w1_load_parm(nfile_data,ierr)
    CLOSE(nfile_data)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX w1_load: namelist file read error: ierr=',ierr
       RETURN
    END IF

    CALL FROPEN(nfile_data,filename_data,0,0,'w1',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX w1_load: data file open error: ierr=',ierr
       RETURN
    END IF
    
    READ(nfile_data) NXMAX,NZMAX,NSMAX

    IF(ALLOCATED(XA)) DEALLOCATE(XA)
    IF(ALLOCATED(XAM)) DEALLOCATE(XAM)
    IF(ALLOCATED(ZA)) DEALLOCATE(ZA)
    ALLOCATE(XA(nxmax+1),XAM(nxmax),ZA(nzmax))
    READ(nfile_data) XA,XAM,ZA

    IF(ALLOCATED(CE2DA)) DEALLOCATE(CE2DA)
    ALLOCATE(CE2DA(nzmax,nxmax+1,3))
    READ(nfile_data) CE2DA

    IF(ALLOCATED(PABS2D)) DEALLOCATE(PABS2D)
    ALLOCATE(PABS2D(nzmax,nxmax,nsmax))
    READ(nfile_data) PABS2D

    IF(ALLOCATED(PROFB)) DEALLOCATE(PROFB)
    IF(ALLOCATED(PROFPN)) DEALLOCATE(PROFPN)
    IF(ALLOCATED(PROFPU)) DEALLOCATE(PROFPU)
    IF(ALLOCATED(PROFTR)) DEALLOCATE(PROFTR)
    IF(ALLOCATED(PROFTP)) DEALLOCATE(PROFTP)
    ALLOCATE(PROFB(nxmax+1))
    ALLOCATE(PROFPN(nxmax+1,nsmax),PROFPU(nxmax+1,nsmax))
    ALLOCATE(PROFTR(nxmax+1,nsmax),PROFTP(nxmax+1,nsmax))
    READ(nfile_data) PROFB,PROFPN,PROFPU,PROFTR,PROFTP

    CLOSE(nfile_data)
 
    RETURN
    
    ierr=9000
    RETURN
  END SUBROUTINE w1_load
    
END MODULE w1file

