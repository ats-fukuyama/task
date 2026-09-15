!  equ_convert.f90

MODULE equ_convert

  PUBLIC equ_read_gfile_csv

CONTAINS
  
  SUBROUTINE equ_read_gfile_csv(knam_gfile_csv)
    USE plcomm
    USE libfio
    USE libgrf
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: knam_gfile_csv
    CHARACTER(LEN=256):: line
    INTEGER:: id_f,ierr,n,n_comma,ndata
    CHARACTER(LEN=8),ALLOCATABLE:: knam_label(:)
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
         z_data,r_data,dz_data,dr_data,current_limit_data
    INTEGER,ALLOCATABLE,DIMENSION(:):: &
         turns_data
    REAL(rkind):: rmin,rmax,zmin,zmax,rlen,zlen
    REAL(rkind):: gpzmin,gpzmax,gprmin,gprmax,gpzlen,gprlen
    REAL(rkind):: gxmin,gxmax,gymin,gymax
    REAL:: guclip

    id_f=21
    CALL FROPEN(id_f,knam_gfile_csv,1,0,'GF',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,A,I8)') &
            'XX file open error in pl_read_gfile_csv:',knam_gfile_csv,ierr
       STOP
    END IF
    
    ! read first line

    READ(id_f,'(A)') line

    ! count comma
    
    n_comma=0
    DO n=n,LEN(line)
       IF(line(n:n).EQ.',') n_comma=n_comma+1
    END DO
    WRITE(6,*) n_comma
    ndata=n_comma
    
    ! read label

    ALLOCATE(knam_label(ndata))
    ALLOCATE(z_data(ndata),r_data(ndata),dz_data(ndata),dr_data(ndata))
    ALLOCATE(turns_data(ndata),current_limit_data(ndata))

    READ(line(5:LEN(line)),*) (knam_label(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(3:LEN(line)),*) (z_data(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(3:LEN(line)),*) (r_data(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(4:LEN(line)),*) (dz_data(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(4:LEN(line)),*) (dr_data(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(7:LEN(line)),*) (turns_data(n),n=1,ndata)

    READ(id_f,'(A)') line
    READ(line(23:LEN(line)),*) (current_limit_data(n),n=1,ndata)

    DO n=1,ndata
       WRITE(6,'(A8,4ES12.4,I8,ES12.4)') &
            knam_label(n),z_data(n),r_data(n),dz_data(n),dr_data(n), &
            turns_data(n),current_limit_data(n)
    END DO

    rmin=r_data(1)-0.5D0*dr_data(1)
    rmax=r_data(1)+0.5D0*dr_data(1)
    zmin=z_data(1)-0.5D0*dz_data(1)
    zmax=z_data(1)+0.5D0*dz_data(1)
    DO n=2,ndata
       rmin=MIN(rmin,r_data(n)-0.5D0*dr_data(n))
       rmax=MAX(rmax,r_data(n)+0.5D0*dr_data(n))
       zmin=MIN(zmin,z_data(n)-0.5D0*dz_data(n))
       zmax=MAX(zmax,z_data(n)+0.5D0*dz_data(n))
    END DO
    rmin=rmin-0.1D0
    rmax=rmax+0.1D0
    zmin=zmin-0.1D0
    zmax=zmax+0.1D0

    rlen=rmax-rmin
    zlen=zmax-zmin
    IF(rlen.GT.zlen) THEN
       gprmin= 1.D0+2.5D0
       gprmax=14.D0+2.5D0
       gprlen=13.D0
       gpzlen=13.D0*zlen/rlen
       gpzmin= 6.5D0-0.5D0*gpzlen
       gpzmax= 6.5D0+0.5D0*gpzlen
    ELSE
       gpzmin= 1.D0
       gpzmax=14.D0
       gpzlen=13.D0
       gprlen=13.D0*rlen/zlen
       gprmin= 6.5D0-0.5D0*gprlen+2.5D0
       gprmax= 6.5D0+0.5D0*gprlen+2.5D0
    ENDIF

!    WRITE(6,'(A,3ES12.4)') 'r  ',rmin,rmax,rlen
!    WRITE(6,'(A,3ES12.4)') 'z  ',zmin,zmax,zlen
!    WRITE(6,'(A,3ES12.4)') 'gpr',gprmin,gprmax,gprlen
!    WRITE(6,'(A,3ES12.4)') 'gpz',gpzmin,gpzmax,gpzlen
    
    CALL pages
    CALL grd2d_frame_start(0,rmin,rmax,zmin,zmax,'@gfile coil data@', &
         gpxmin=gprmin,gpxmax=gprmax,gpymin=gpzmin,gpymax=gpzmax)
    DO n=1,ndata
       gxmin=r_data(n)-0.5D0*dr_data(n)
       gxmax=r_data(n)+0.5D0*dr_data(n)
       gymin=z_data(n)-0.5D0*dz_data(n)
       gymax=z_data(n)+0.5D0*dz_data(n)
       CALL MOVE2D(GUCLIP(gxmin),GUCLIP(gymin))
       CALL DRAW2D(GUCLIP(gxmax),GUCLIP(gymin))
       CALL DRAW2D(GUCLIP(gxmax),GUCLIP(gymax))
       CALL DRAW2D(GUCLIP(gxmin),GUCLIP(gymax))
       CALL DRAW2D(GUCLIP(gxmin),GUCLIP(gymin))
       CALL GTEXT2D(GUCLIP(r_data(n)+0.5D0*dr_data(n)),GUCLIP(z_data(n)), &
            knam_label(n),LEN(knam_label(n)),0)
    END DO
    CALL grd2d_frame_end
    CALL pagee
  END SUBROUTINE equ_read_gfile_csv
END MODULE equ_convert
