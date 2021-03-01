!
!
!
MODULE T2COUT

  USE T2COMM, ONLY: &
         ikind,rkind

  PRIVATE
  PUBLIC T2COUT_EXECUTE
 
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: rhonrho,chinchi
  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: nlnrho,nnnrho
  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: nrhonl,nnnl
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: chinl
  INTEGER(ikind):: nlmax_c,nnmax_c
  
CONTAINS
  
  SUBROUTINE T2COUT_EXECUTE

    USE T2COMM,ONLY:&
         & NVMAX,NDMAX,NKMAX,NMMAX,&
         & MassScaCoef,AdveVecCoef,AdveTenCoef,DiffTenCoef,&
         & GradVecCoef,GradTenCoef,ExciScaCoef,ExciVecCoef,&
         & ExciTenCoef,SourScaCoef
 
    INTEGER(ikind)::&
         i_d,i_k,i_v,&
         j_d,j_k
    REAL(   rkind)::d2coef(1:NVMAX,1:NMMAX)
    CHARACTER(LEN=2)::c2coef

    WRITE(6,*)'********** COEFFICIENT CHECK START**********'

    DO
       WRITE(6,*)'i_v'
       READ(5,*)i_v
       
       IF(     i_v.EQ.0)THEN
          EXIT
       ELSEIF((i_v.GT.NVMAX).OR.(i_v.LT.1))THEN
          CYCLE
       ELSE
          DO
             WRITE(6,*)'ms,av,at,dt,gv,gt,es,ev,et xx/exit'
             READ(5,*)c2coef
             SELECT CASE (c2coef)
             CASE ('ms')
                d2coef(           1:NVMAX,1:NMMAX)&
                     =MassScaCoef(i_v,1:NVMAX,1:NMMAX)
                CALL T2_COUT(d2coef)
             CASE ('av')
                DO 
                   WRITE(6,*)'i_d'
                   READ(5,*)  i_d
                   IF((i_d.GT.NDMAX).OR.(i_d.LT.1)) CYCLE
                   
                   d2coef(                   1:NVMAX,1:NMMAX)&
                        =AdveVecCoef(i_d,i_v,1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                   
                ENDDO
             CASE ('at')
                DO 
                   WRITE(6,*)'i_k,i_d,j_d'
                   READ(5,*)  i_k,i_d,j_d
                   IF(  (i_k.GT.NKMAX).OR.(i_k.LT.1).OR.&
                        (i_d.GT.NDMAX).OR.(i_d.LT.1).OR.&
                        (j_d.GT.NDMAX).OR.(j_d.LT.1)) CYCLE
                   
                   d2coef(                1:NVMAX,1:NMMAX)&
                        = AdveTenCoef(i_d,j_d,i_k,i_v,&
                        &                 1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('dt')
                DO 
                   WRITE(6,*)'i_d,j_d'
                   READ(5,*)  i_d,j_d
                   IF(  (i_d.GT.NDMAX).OR.(i_d.LT.1).OR.&
                        (j_d.GT.NDMAX).OR.(j_d.LT.1)) CYCLE
                   d2coef(                       1:NVMAX,1:NMMAX)&
                        =DiffTenCoef(i_d,j_d,i_v,1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('gv')
                DO 
                   WRITE(6,*)'i_d'
                   READ(5,*)  i_d
                   IF((i_d.GT.NDMAX).OR.(i_d.LT.1)) CYCLE
                   
                   d2coef(                    1:NVMAX,1:NMMAX)&
                        = GradVecCoef(i_d,i_v,1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('gt')
                DO 
                   WRITE(6,*)'i_k,i_d,j_d'
                   READ(5,*)  i_k,i_d,j_d
                   IF(  (i_k.GT.NKMAX).OR.(i_k.LT.1).OR.&
                        (i_d.GT.NDMAX).OR.(i_d.LT.1).OR.&
                        (j_d.GT.NDMAX).OR.(j_d.LT.1)) CYCLE
                   
                   d2coef(            1:NVMAX,1:NMMAX) &
                        = GradTenCoef(i_d,j_d,i_k,i_v,&
                        &             1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('es')
                d2coef(                1:NVMAX,1:NMMAX)&
                     = ExciScaCoef(i_v,1:NVMAX,1:NMMAX)
                CALL T2_COUT(d2coef)
             CASE ('ev')
                DO
                   WRITE(6,*)'i_k,i_d'
                   READ(5,*)  i_k,i_d
                   IF(  (i_k.GT.NKMAX).OR.(i_k.LT.1).OR.&
                        (i_d.GT.NDMAX).OR.(i_d.LT.1)) CYCLE
                   d2coef(                1:NVMAX,1:NMMAX) &
                        = ExciVecCoef(i_d,i_k,i_v,&
                        &                 1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('et')
                DO 
                   WRITE(6,*)'i_k,j_k,i_d,j_d'
                   READ(5,*)  i_k,j_k,i_d,j_d
                   IF(  (i_k.GT.NKMAX).OR.(i_k.LT.1).OR.&
                        (j_k.GT.NKMAX).OR.(j_k.LT.1).OR.&
                        (i_d.GT.NDMAX).OR.(i_d.LT.1).OR.&
                        (j_d.GT.NDMAX).OR.(j_d.LT.1)) CYCLE
                   
                   d2coef(                1:NVMAX,1:NMMAX)&
                        = ExciTenCoef(i_d,j_d,i_k,j_k,i_v,&
                        &                 1:NVMAX,1:NMMAX)
                   CALL T2_COUT(d2coef)
                   EXIT
                ENDDO
             CASE ('ss')
                d2coef(                1:NVMAX,1:NMMAX)&
                     = SourScaCoef(i_v,1:NVMAX,1:NMMAX) 
                CALL T2_COUT(d2coef)
             CASE ('xx')
                EXIT
             CASE DEFAULT
                CYCLE
             END SELECT
          ENDDO
       ENDIF
    ENDDO
    
    WRITE(6,*)'********** COEFFICIENT CHECK END  **********'
    RETURN
    
  END SUBROUTINE T2COUT_EXECUTE


  SUBROUTINE T2_COUT(d2cm)
    

    USE T2PARM,ONLY: T2_PARM
    USE T2COMM,ONLY: NMMAX,NVMAX,NXMAX,i2crt
    USE libkio
    USE libgrf
    IMPLICIT NONE
    
    REAL(rkind),DIMENSION(1:NVMAX,1:NMMAX),INTENT(IN)::d2cm
    REAL(rkind),DIMENSION(1:NVMAX,1:NXMAX)::d2cx
 
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER(LEN=80) :: line,kw
    CHARACTER(LEN=1) :: kid,kch
    INTEGER(ikind) :: nwmax,iloc0,nw,i,ich0,nch,ich,j
    INTEGER(ikind) :: i_m,i_x
    CHARACTER(LEN=80),DIMENSION(40):: kword,kwid,knum
    INTEGER(ikind),DIMENSION(40):: inum
    
    DO i_m =1, NMMAX
       i_x = i2crt(2,i_m)
       d2cx(1:NVMAX,i_x)= d2cm(1:NVMAX,i_m)
    ENDDO
    
    CALL T2_GSETUP

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '#### T2 GOUT: Rn An Cn Pn Bn RA AA X/exit'
    CALL TASK_KLIN(line,kid,mode,T2_PARM)
    IF(mode == 2 .OR. mode == 3) GOTO 1

!   --- separate words in a line ---
    ILOC0=1
    NWMAX=0
    DO I=1,80
       IF(ILOC0==0) THEN
          IF(LINE(I:I)/=' ' .AND. LINE(I:I)/=',') THEN
             IF(NWMAX/=0) THEN
                ILOC0=I
             END IF
          END IF
       ELSE
          IF(LINE(I:I)==' ' .OR. LINE(I:I)==',') THEN
             NWMAX=NWMAX+1
             KWORD(NWMAX)=LINE(ILOC0:I)
             ILOC0=0
          END IF
       END IF
    END DO
    IF(ILOC0 /= 0) THEN
       NWMAX=NWMAX+1
       KWORD(NWMAX)=LINE(ILOC0:I)
    END IF

!    DO NW=1,NWMAX
!       WRITE(6,'(A,I5,4X,A)') 'NW,KWORD=',NW,TRIM(KWORD(NW))
!    END DO

!   --- separate id and number in a word ---

    DO NW=1,NWMAX
       KW=KWORD(NW)
       ICH0=0
       DO NCH=1,LEN(KW)
!         --- lower-case char to upper-case char ---
          ICH=ICHAR(KW(NCH:NCH))
          IF(ICH.GE.97.AND.ICH.LE.122) ICH=ICH-32
          KCH=CHAR(ICH)
          KW(NCH:NCH)=KCH

          IF(KCH.GE.'A'.AND.KCH.LE.'Z') THEN
             ICH0=ICH0+1
          ELSE
             EXIT
          ENDIF
       END DO
       IF(ICH0.EQ.0) THEN
          KWID(NW)=' '
          KNUM(NW)=KW
       ELSE
          KWID(NW)=KW(1:ICH0)
          KNUM(NW)=KW(ICH0+1:LEN(KW))
       END IF
       INUM(NW)=0
       READ(KNUM(NW),*,ERR=901,END=901) INUM(NW)
901    CONTINUE

    END DO

    IF(KWID(1)=='X') GO TO 9000

    CALL PAGES
    DO NW=1,NWMAX
!       WRITE(6,'(A,A,A,A,A,I5)') &
!            'KID,KNUM=:',TRIM(KWID(NW)),':',TRIM(KNUM(NW)),':',INUM(NW)
       SELECT CASE(LEN_TRIM(KWID(NW)))
       CASE(1)
          SELECT CASE(KWID(NW))
          CASE('R')
             CALL T2_CR(INUM(NW),1, 0,d2cx)
          CASE('A')
             CALL T2_CR(INUM(NW),11,0,d2cx)
          CASE('C')
             CALL T2_CC(INUM(NW),1, 0,d2cx)
          CASE('D')
             CALL T2_CC(INUM(NW),2, 0,d2cx)
          CASE('P')
             CALL T2_CC(INUM(NW),3, 0,d2cx)
          CASE('B')
             CALL T2_CC(INUM(NW),15,0,d2cx)
          CASE('Y')
             CALL T2_CC(INUM(NW),11,0,d2cx)
          END SELECT
       CASE(2)
          SELECT CASE(KWID(NW))
          CASE('RA')
             DO J=1,5
                DO I=5*(J-1)+1,5*(J-1)+5
                   CALL T2_CR(I, 1, 29+I,d2cx)
                END DO
             END DO
          CASE('AA')
             DO J=1,5
                DO I=5*(J-1)+1,5*(J-1)+5
                   CALL T2_CR(I,11, 29+I,d2cx)
                END DO
             END DO
          CASE('CA')
             DO J=1,5
                DO I=5*(J-1)+1,5*(J-1)+5
                   CALL T2_CC(I, 1, 29+I,d2cx)
                END DO
             END DO
          CASE('PA')
             DO J=1,5
                DO I=5*(J-1)+1,5*(J-1)+5
                   CALL T2_CC(I, 3, 29+I,d2cx)
                END DO
             END DO
          CASE('BA')
             DO J=1,5
                DO I=5*(J-1)+1,5*(J-1)+5
                   CALL T2_CC(I,15, 29+I,d2cx)
                END DO
             END DO
          END SELECT
       END SELECT
    END DO
    CALL PAGEE
    
    GO TO 1

9000 CONTINUE
    CALL T2_GRELEASE

    RETURN

  END SUBROUTINE T2_COUT

  SUBROUTINE T2_GSETUP
    USE T2COMM, ONLY: &
         twopi,NXMAX,NVMAX, &
         NLMAX,NPMIN,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax
    IMPLICIT NONE
    INTEGER(ikind):: nchi,nl,nrho,nchimaxl,nr,ierr
    REAL(rkind):: dchi,drho
    REAL(rkind):: rho_temp!added by H. SETO

    nlmax_c=NLMAX
    nchimax=NPMIN*2**(i1mlvl(nlmax_c)-1)
    ALLOCATE(nrhonl(0:nlmax_c),nnnl(0:nlmax_c))
    ALLOCATE(chinchi(nchimax+1),chinl(nchimax+1,nlmax_c))

    dchi=twopi/nchimax
    DO nchi=1,nchimax+1
       chinchi(nchi)=dchi*(nchi-1)
    END DO

    nrhonl(0)=1
    nnnl(0)=1
    nrhomax=2
    nnmax_c=2
    DO nl=1,nlmax_c
       nrhonl(nl)=nrhomax
       nnnl(nl)=nnmax_c
       nchimaxl=NPMIN*2**(i1mlvl(nl)-1)
       nrhomax=nrhomax+i1rdn2(nl)
       nnmax_c=nnmax_c+nchimaxl*i1rdn2(nl)
    END DO
    nrhomax=nrhomax-1
    nnmax_c=nnmax_c-1
    AlLOCATE(nlnrho(nrhomax),rhonrho(nrhomax),nnnrho(nrhomax))

    nrho=1
    nlnrho(nrho)=0
    rhonrho(nrho)=0.D0
    nnnrho(nrho)=1

    DO nl=1,nlmax_c
       drho=(d1rec(nl)-d1rec(nl-1))/i1rdn2(nl)
       nchimaxl=NPMIN*2**(i1mlvl(nl)-1)
       dchi=twopi/nchimaxl
       DO nchi=1,nchimaxl+1
          chinl(nchi,nl)=dchi*(nchi-1)
       END DO
       DO nr=1,i1rdn2(nl)
          nrho=nrhonl(nl)+(nr-1)
          nlnrho(nrho)=nl
          rhonrho(nrho)=d1rec(nl-1)+drho*nr
          nnnrho(nrho)=nnnl(nl)+nchimaxl*(nr-1)
!          write(6,'(A,4I5,1PE12.4)') &
!               'nl,nr,nrho,nnnrho(nrho),rhonrho(nrho)=', &
!               nl,nr,nrho,nnnrho(nrho),rhonrho(nrho)
       END DO
    END DO

    RETURN

  END SUBROUTINE T2_GSETUP

  SUBROUTINE T2_GRELEASE

    DEALLOCATE(nrhonl,nnnl)
    DEALLOCATE(chinchi,chinl)
    DEAlLOCATE(nlnrho,rhonrho,nnnrho)

    RETURN
  END SUBROUTINE T2_GRELEASE

  SUBROUTINE T2_CC(INUM,ID,NGP,D2CX)
    USE T2COMM, ONLY: & ! changed by 2014-02-05 H.Seto 
         twopi,NXMAX,NVMAX, &
         NLMAX,NPMIN,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax,d0rw
    USE libspl2d
    USE libgrf,ONLY: GRD2D
    IMPLICIT NONE
    INTEGER,PARAMETER:: nxmax_c=41,nymax_c=41
    INTEGER(ikind),INTENT(IN):: inum,id,ngp
    REAL(rkind),DIMENSION(1:NVMAX,1:NXMAX),INTENT(IN)::d2cx
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: gz
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: gzl,dgzl,chig
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: ugzl
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: gx,gy
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: ddx,ddy,ddxy,gxy
    REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: ugz
    INTEGER(ikind):: nchi,nl,nrho,nchimaxl,nr,ierr,nchig,nx,ny
    REAL(rkind):: dchig,xmin,xmax,ymin,ymax,dx,dy,x,y,r,th
    CHARACTER(LEN=80):: LINE

    nchig=MAX(nchimax,72)
    ALLOCATE(gz(nrhomax,nchig+1),chig(nchig+1))
    ALLOCATE(gzl(nchimax+1),dgzl(nchimax+1),ugzl(4,nchimax+1))
    dchig=TWOPI/nchig
    DO nchi=1,nchig+1
       chig(nchi)=dchig*(nchi-1)
       gz(1,nchi)=d2cx(inum,1)! changed by 2014-02-05 H.SETO
       IF(gz(1,nchi).GT. 1.D10) gz(1,nchi)= 1.D10
       IF(gz(1,nchi).LT.-1.D10) gz(1,nchi)=-1.D10
    END DO
    DO nrho=2,nrhomax
       nl=nlnrho(nrho)
       nchimaxl=NPMIN*2**(i1mlvl(nl)-1)
       DO nchi=1,nchimaxl
          gzl(nchi)=d2cx(inum,nnnrho(nrho)+nchi-1)! changed by 2014-02-05 H.SETO
          IF(gzl(nchi).GT. 1.D10) gzl(nchi)= 1.D10
          IF(gzl(nchi).LT.-1.D10) gzl(nchi)=-1.D10
       END DO
       gzl(nchimaxl+1)=gzl(1)
       CALL SPL1D(chinl(1:nchimaxl+1,nl),gzl,dgzl,ugzl,nchimaxl+1,4,ierr)
       DO nchi=1,nchig+1
          CALL SPL1DF(chig(nchi),gz(nrho,nchi), &
                      chinl(1:nchimaxl+1,nl),ugzl,nchimaxl+1,ierr)
       END DO
    END DO

    IF(ID.GT.10) THEN
       ALLOCATE(ddx(nrhomax,nchig),ddy(nrhomax,nchig),ddxy(nrhomax,nchig))
       ALLOCATE(ugz(4,4,nrhomax,nchig))
       ALLOCATE(gx(nxmax_c),gy(nymax_c),gxy(nxmax_c,nymax_c))

       CALL SPL2D(rhonrho,chig,gz,ddx,ddy,ddxy,ugz, &
                  nrhomax,nrhomax,nchig+1,0,0,ierr)
       xmin=-rhonrho(nrhomax)
       xmax= rhonrho(nrhomax)
       ymin=-rhonrho(nrhomax)
       ymax= rhonrho(nrhomax)
       dx=(xmax-xmin)/(nxmax_c-1)
       dy=(ymax-ymin)/(nymax_c-1)
       DO nx=1,nxmax_c
          gx(nx)=xmin+dx*(nx-1)
       END DO
       DO ny=1,nymax_c
          gy(ny)=ymin+dy*(ny-1)
       END DO
       DO nx=1,nxmax_c
          x=gx(nx)
          DO ny=1,nymax_c
             y=gy(ny)
             r=SQRT(x*x+y*y)
             th=ATAN2(y,x)      ! for compatibility with contour
             IF(th.LT.0.D0) th=th+TWOPI
             IF(r >= d0rw) THEN
                gxy(nx,ny)=0.D0
             ELSE
                CALL SPL2DF(r,th,gxy(nx,ny),rhonrho,chig,ugz, &
                            nrhomax,nrhomax,nchig,ierr)
             END IF
!             WRITE(6,'(2I5,1P5E12.4)') nx,ny,x,y,r,th,gxy(nx,ny)
          END DO
       END DO
    END IF
             
    WRITE(LINE,'(A,I3,A)') '@diguv(',inum,')@'

    SELECT CASE(ID)
    CASE(1)
       CALL GRD2D(ngp,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=4,TITLE_SIZE=0.4D0)
    CASE(2)
       CALL GRD2D(ngp,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=1,TITLE_SIZE=0.4D0)
    CASE(3)
       CALL GRD2D(ngp,rhonrho,chig,gz,nrhomax,nrhomax,nchig, &
                  TITLE=LINE,MODE_XY=1,MODE_2D=2,TITLE_SIZE=0.4D0)
    CASE(11)
       CALL GRD2D(ngp,gx,gy,gxy,nxmax_c,nxmax_c,nymax_c, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=4,TITLE_SIZE=0.4D0, &
                  XMIN=-d0rw,XMAX=d0rw,YMIN=-d0rw,YMAX=d0rw, &
                  ASPECT=1.D0)
    CASE(12)
       CALL GRD2D(ngp,gx,gy,gxy,nxmax_c,nxmax_c,nymax_c, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=1,TITLE_SIZE=0.4D0,&
                  XMIN=-d0rw,XMAX=d0rw,YMIN=-d0rw,YMAX=d0rw, &
                  ASPECT=1.D0)
    CASE(13)
       CALL GRD2D(ngp,gx,gy,gxy,nxmax_c,nxmax_c,nymax_c, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=2,TITLE_SIZE=0.4D0, &
                  XMIN=-d0rw,XMAX=d0rw,YMIN=-d0rw,YMAX=d0rw, &
                  ASPECT=1.D0)
    CASE(14)
       CALL GRD2D(ngp,gx,gy,gxy,nxmax_c,nxmax_c,nymax_c, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=11,TITLE_SIZE=0.4D0)
    CASE(15)
       CALL GRD2D(ngp,gx,gy,gxy,nxmax_c,nxmax_c,nymax_c, &
                  TITLE=LINE,MODE_XY=0,MODE_2D=12,TITLE_SIZE=0.4D0)
    END SELECT

    IF(ID.GT.10) THEN
       DEALLOCATE(ddx,ddy,ddxy)
       DEALLOCATE(ugz)
       DEALLOCATE(gx,gy,gxy)
    END IF
    RETURN
  END SUBROUTINE T2_CC

  SUBROUTINE T2_CR(INUM,ID,NGP,D2CX)
    USE libgrf,ONLY: GRD1D
    
    USE T2COMM, ONLY: &! changed by 2014-02-05 H.SETO
         twopi,NXMAX,NVMAX, &
         NLMAX,NPMIN,i1mlvl,i1rdn2,d1rec, &
         nrhomax,nchimax

    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN):: inum,id,ngp
    REAL(rkind),DIMENSION(1:NVMAX,1:NXMAX),INTENT(IN):: d2cx
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: gz
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: gzl,dgzl,ga
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: ugzl
    INTEGER(ikind):: nchi,nl,nrho,nchimaxl,nr,ierr
    CHARACTER(LEN=80):: LINE

    ALLOCATE(gz(nrhomax,nchimax+1))
    ALLOCATE(gzl(nchimax+1),dgzl(nchimax+1),ugzl(4,nchimax+1))
    ALLOCATE(ga(nrhomax))
    DO nchi=1,nchimax+1
       gz(1,nchi)=d2cx(inum,1)! changed by 2014-02-05 H.SETO
       IF(gz(1,nchi).GT. 1.D10) gz(1,nchi)= 1.D10
       IF(gz(1,nchi).LT.-1.D10) gz(1,nchi)=-1.D10
    END DO
    ga(1)=d2cx(inum,1)! changed by 2014-02-05 H.SETO
    DO nrho=2,nrhomax
       nl=nlnrho(nrho)
       nchimaxl=NPMIN*2**(i1mlvl(nl)-1)
       DO nchi=1,nchimaxl
          gzl(nchi)=d2cx(inum,nnnrho(nrho)+nchi-1)
          IF(gzl(nchi).GT. 1.D10) gzl(nchi)= 1.D10
          IF(gzl(nchi).LT.-1.D10) gzl(nchi)=-1.D10
       END DO
       SELECT CASE(inum)
       CASE(1:3)
          ga(nrho)=gzl(nchimaxl)
       CASE DEFAULT
          ga(nrho)=0.D0
          DO nchi=1,nchimaxl
             ga(nrho)=ga(nrho)+gzl(nchi)
          END DO
          ga(nrho)=ga(nrho)/nchimaxl
       END SELECT

       IF(nchimaxl==nchimax) THEN
          DO nchi=1,nchimaxl
             gz(nrho,nchi)=gzl(nchi)
          END DO
       ELSE
          gzl(nchimaxl+1)=gzl(1)
          CALL SPL1D(chinl(1:nchimaxl+1,nl),gzl,dgzl,ugzl,nchimaxl+1,4,ierr)
          DO nchi=1,nchimax
             CALL SPL1DF(chinchi(nchi),gz(nrho,nchi), &
                         chinl(1:nchimaxl+1,nl),ugzl,nchimaxl+1,ierr)
          END DO
       END IF
    END DO

    WRITE(LINE,'(A,I3,A)') '@diguv(',inum,')@'

    SELECT CASE(id)
    CASE(1)
       CALL GRD1D(NGP,rhonrho,gz,nrhomax,nrhomax,nchimax,TITLE=LINE)
    CASE(11)
       CALL GRD1D(NGP,rhonrho,ga,nrhomax,nrhomax,1,TITLE=LINE)
    END SELECT

    RETURN
  END SUBROUTINE T2_CR

END Module T2COUT
