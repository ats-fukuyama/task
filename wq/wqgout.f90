! wqgout.f90

MODULE wqgout

  PRIVATE
  PUBLIC wq_gout

CONTAINS

  SUBROUTINE wq_gout
    USE wqcomm
    USE wqparm,ONLY: wq_parm,wq_broadcast
    USE wqview,ONLY: wq_view
    USE libchar
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=1):: kid
    INTEGER:: mode

1   CONTINUE
    
    IF(nrank.EQ.0) THEN
       WRITE(6,'(A)') &
            '## wqgout menu: A,B,C  parm  X/end'
       CALL TASK_KLIN(line,kid,mode,wq_parm)
    END IF
    CALL ToUpper(kid)
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)

    IF(mode.EQ.2) CALL wq_broadcast
    IF(mode.NE.1) GO TO 1

    SELECT CASE(kid)
    CASE('A')
       CALL wq_gsub(0)
    CASE('B')
       CALL wq_gsub(1)
    CASE('C')
       CALL wq_gsub_plot
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gout

  SUBROUTINE wq_gsub(flag)

  use wqcomm
  USE libgrf
  IMPLICIT NONE
  INTEGER(4)            :: nx,ny,k,NGP,MODE,IPRD,IG2D,flag
  CHARACTER(LEN=80)     :: STR
  REAL(8)               :: FX(nxmax),FY(nymax),FZ(nxmax,nymax)

  DO nx=1,nxmax
!        FX(nx)=RR+dx*DBLE(nx-1-(nxmax-1)/2)
       FX(nx)=dx*DBLE(nx-1)
  ENDDO
  DO ny=1,nymax
!      FY(ny)=dy*DBLE(ny-1-(nymax-1)/2)
       FY(ny)=dy*DBLE(ny-1)
  ENDDO

  SELECT CASE(flag)
  CASE(0)

  call PAGES
  
  MODE=0
  IPRD=0
  IG2D=1

  NGP=5
  STR='/AP/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=pabs(nx,ny)
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=6
  STR='/OUH/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=OUH(nx,ny)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=7
  STR='/R/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=OR(nx,ny)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=8
  STR='/L/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=OL(nx,ny)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=9
  STR='/OCE/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=OCE(nx,ny)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=10
  STR='/OPE/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=OPE(nx,ny)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=11
  STR='/ne/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=ne(nx,ny)
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  call PAGEE

  CASE(1)

  call PAGES

  MODE=0
  IPRD=0
  IG2D=1

  NGP=5
  STR='/Real(EX)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=REAL(EX(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)

  NGP=8
  STR='/Image(EX)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=IMAG(EX(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)
  
  NGP=6
  STR='/Real(EY)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=Real(EY(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)

  NGP=9
  STR='/Image(EY)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=IMAG(EY(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)
  
  NGP=7
  STR='/Real(EZ)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=Real(EZ(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)

  NGP=10
  STR='/Image(EZ)/'
  DO ny=1,nymax
     DO nx=1,nxmax
        FZ(nx,ny)=IMAG(EZ(nx,ny))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,STR,MODE,IPRD,IG2D)
  
  CALL PAGEE

   END SELECT
  
  RETURN
  END SUBROUTINE wq_gsub

  SUBROUTINE wq_gsub_plot

    use wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER            :: nx,ny,nt,NGP,np
    CHARACTER(LEN=80)  :: title
    REAL(rkind)        :: Emax,Eabs
    REAL(rkind)        :: FX(nxmax),FY(nymax),FZ(nxmax,nymax)

    Emax=0.D0
    DO np=1,ntplot
       DO ny=1,nymax
          DO nx=1,nxmax
             Eabs=ABS(EX_save(nx,ny,np))**2 &
                 +ABS(EY_save(nx,ny,np))**2 &
                 +ABS(EZ_save(nx,ny,np))**2
             Emax=MAX(Emax,SQRT(Eabs))
          END DO
       END DO
    END DO
    
    IF(ntplot.LE.0) THEN
       WRITE(6,*) 'XX wq_gsub_plot: ntplot.LE.0'
       RETURN
    END IF
  
    DO nx=1,nxmax
!        FX(nx)=RR+dx*DBLE(nx-1-(nxmax-1)/2)
       FX(nx)=dx*DBLE(nx-1)
    ENDDO
    DO ny=1,nymax
!      FY(ny)=dy*DBLE(ny-1-(nymax-1)/2)
       FY(ny)=dy*DBLE(ny-1)
    ENDDO

    CALL PAGES

    NGP=1
    title='/EX plot/'
    DO ny=1,nymax
       DO nx=1,nxmax
          FZ(nx,ny)=EX_save(nx,ny,1)
       END DO
    END DO
    CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
               ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    DO nt=2,ntplot
       DO ny=1,nymax
          DO nx=1,nxmax
             FZ(nx,ny)=EX_save(nx,ny,nt)
          END DO
       END DO
       CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
            NOFRAME=1,NOTITLE=1,NOINFO=1, &
            ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    END DO
    
    NGP=2
    title='/EY plot/'
    DO ny=1,nymax
       DO nx=1,nxmax
          FZ(nx,ny)=EY_save(nx,ny,1)
       END DO
    END DO
    CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
               ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    DO nt=2,ntplot
       DO ny=1,nymax
          DO nx=1,nxmax
             FZ(nx,ny)=EY_save(nx,ny,nt)
          END DO
       END DO
       CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
            NOFRAME=1,NOTITLE=1,NOINFO=1, &
            ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    END DO

    NGP=3
    title='/EZ plot/'
    DO ny=1,nymax
       DO nx=1,nxmax
          FZ(nx,ny)=EZ_save(nx,ny,1)
       END DO
    END DO
    CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
               ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    DO nt=2,ntplot
       DO ny=1,nymax
          DO nx=1,nxmax
             FZ(nx,ny)=EZ_save(nx,ny,nt)
          END DO
       END DO
       CALL GRD2D(NGP,FX,FY,FZ,nxmax,nxmax,nymax,title, &
            NOFRAME=1,NOTITLE=1,NOINFO=1, &
            ASPECT=0.D0,FMIN=-Emax,FMAX=Emax)
    END DO

    CALL PAGEE

    RETURN
  END SUBROUTINE wq_gsub_plot
END MODULE wqgout
