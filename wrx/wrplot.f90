! wrplot.f90

PROGRAM wrplot

  USE bpsd_kinds
  USE bpsd_constants

  USE libgrf
  IMPLICIT NONE
  INTEGER:: id
  EXTERNAL GSOPEN,GSCLOS

  CALL GSOPEN

1 CONTINUE
  WRITE(6,*) '## mode input: 1, 2 or 9 ?'
  READ(5,*,END=9000,ERR=1) id

  SELECT CASE(id)
  CASE(1)
     CALL wrplot1
  CASE(2)
     CALL wrplot2
  CASE(9)
     GO TO 9000
  CASE DEFAULT
     GO TO 1
  END SELECT
  GO TO 1

9000 CONTINUE
  CALL GSCLOS
  STOP

CONTAINS

  SUBROUTINE wrplot1
    
    USE libgrf
    IMPLICIT NONE
    EXTERNAL PAGES,PAGEE
    INTEGER:: nf,nfmax,nc,ncmax,nfmaxl,np,npmax,line_pat
    REAL(rkind):: fmin,fmax,rnp,rnpmin,rnpmax
    REAL(rkind):: drnp,df,f,p0,pr2,x,y
    REAL(rkind),ALLOCATABLE:: R_ny(:),F_ny(:,:)
    CHARACTER(LEN=20):: title

    title='@w_c/w vs p_para/c@'

    fmin=0.4D0
    fmax=2.4D0
    rnpmin=0.0D0
    rnpmax=0.7D0
    nfmax=301
    ncmax=3
    npmax=8

1   CONTINUE
    WRITE(6,'(A)') &
         '## Input fmin,fmax,rnpmin,rnpmax,nfmax,npmax,ncmax ?'
    READ(5,*,ERR=1,END=9000) fmin,fmax,rnpmin,rnpmax,nfmax,npmax,ncmax

    IF(ALLOCATED(R_ny)) DEALLOCATE(R_ny)
    IF(ALLOCATED(F_ny)) DEALLOCATE(F_ny)
    ALLOCATE(R_ny(nfmax),F_ny(nfmax,3))

    CALL PAGES
    CALL grd2d_frame_start(0,-1.D0,1.D0,fmin,fmax,title,NOINFO=1)
    drnp=(rnpmax-rnpmin)/(npmax-1)
    df=(fmax-fmin)/(nfmax-1)
    DO np=1,npmax
       CALL SETLIN(0,2,7-MOD(np-1,5))
       rnp=rnpmin+drnp*(np-1)
       DO nc=1,ncmax
          SELECT CASE(ABS(nc))
          CASE(1)
             line_pat=0
          CASE(2)
             line_pat=2
          CASE(3)
             line_pat=4
          CASE DEFAULT
             line_pat=6
          END SELECT
          nfmaxl=nfmax
          DO nf=1,nfmax
             f=fmin+df*(nf-1)
             R_ny(nf)=f
             p0=rnp*nc/f/(1.D0-rnp**2)
             pr2=((nc/f)**2-(1.D0-rnp**2))/(1.D0-rnp**2)**2
             IF(pr2.GE.0.D0) THEN
                F_ny(nf,1)=p0
                F_ny(nf,2)=p0-SQRT(pr2)
                F_ny(nf,3)=p0+SQRT(pr2)
                nfmaxl=nf
             END IF
          END DO

          y=R_ny(1)
          x=F_ny(1,2)
          CALL MOVEPT2D(gdclip(x),gdclip(y),line_pat)
          DO nf=2,nfmaxl
             y=R_ny(nf)
             x=F_ny(nf,2)
             CALL DRAWPT2D(gdclip(x),gdclip(y))
          END DO
          DO nf=nfmaxl,1,-1
             y=R_ny(nf)
             x=F_ny(nf,3)
             CALL DRAWPT2D(gdclip(x),gdclip(y))
          END DO
       END DO
    END DO
    CALL grd2d_frame_end
    CALL PAGEE
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wrplot1

  SUBROUTINE wrplot2
    
    USE libgrf
    IMPLICIT NONE
    EXTERNAL PAGES,PAGEE
    INTEGER:: nf,nfmax,nc,ncmax,nfmaxl,np,npmax,line_pat
    REAL(rkind):: fmin,fmax,rnp,rnpmin,rnpmax
    REAL(rkind):: drnp,df,f,p0,pr2,x,y
    REAL(rkind),ALLOCATABLE:: R_ny(:),F_ny(:,:)
    CHARACTER(LEN=20):: title

    title='@w_c/w vs p_para/c@'

    fmin=0.4D0
    fmax=2.4D0
    rnpmin=0.0D0
    rnpmax=0.7D0
    nfmax=301
    ncmax=2
    npmax=8

1   CONTINUE
    WRITE(6,'(A)') &
         '## Input fmin,fmax,rnpmin,rnpmax,nfmax,npmax,ncmax ?'
    READ(5,*,ERR=1,END=9000) fmin,fmax,rnpmin,rnpmax,nfmax,npmax,ncmax

    IF(ALLOCATED(R_ny)) DEALLOCATE(R_ny)
    IF(ALLOCATED(F_ny)) DEALLOCATE(F_ny)
    ALLOCATE(R_ny(nfmax),F_ny(nfmax,3))

    CALL PAGES
    CALL grd2d_frame_start(0,-1.D0,1.D0,fmin,fmax,title,NOINFO=1)
    drnp=(rnpmax-rnpmin)/(npmax-1)
    df=(fmax-fmin)/(nfmax-1)
    DO np=1,npmax
       CALL SETLIN(0,2,7-MOD(np-1,5))
       DO nc=1,ncmax
          SELECT CASE(ABS(nc))
          CASE(1)
             line_pat=0
             rnp=rnpmin+drnp*(np-1)
          CASE(2)
             line_pat=2
             rnp=-rnpmin-drnp*(np-1)
          CASE(3)
             line_pat=4
             rnp=rnpmin+drnp*(np-1)
          CASE DEFAULT
             line_pat=6
             rnp=rnpmin+drnp*(np-1)
          END SELECT
          nfmaxl=nfmax
          DO nf=1,nfmax
             f=fmin+df*(nf-1)
             R_ny(nf)=f
             p0=rnp*nc/f/(1.D0-rnp**2)
             pr2=((nc/f)**2-(1.D0-rnp**2))/(1.D0-rnp**2)**2
             IF(pr2.GE.0.D0) THEN
                F_ny(nf,1)=p0
                F_ny(nf,2)=p0-SQRT(pr2)
                F_ny(nf,3)=p0+SQRT(pr2)
                nfmaxl=nf
             END IF
          END DO

          y=R_ny(1)
          x=F_ny(1,2)
          CALL MOVEPT2D(gdclip(x),gdclip(y),line_pat)
          DO nf=2,nfmaxl
             y=R_ny(nf)
             x=F_ny(nf,2)
             CALL DRAWPT2D(gdclip(x),gdclip(y))
          END DO
          DO nf=nfmaxl,1,-1
             y=R_ny(nf)
             x=F_ny(nf,3)
             CALL DRAWPT2D(gdclip(x),gdclip(y))
          END DO
       END DO
    END DO
    CALL grd2d_frame_end
    CALL PAGEE
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wrplot2

END PROGRAM wrplot
