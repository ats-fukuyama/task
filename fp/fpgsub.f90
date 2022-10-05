! fpgsub.f90

MODULE fpgsub

  USE fpcomm
  USE fpcont
  USE fpfout
  USE libgrf
  USE libmpi

  PRIVATE
  PUBLIC fp_grac

CONTAINS

  SUBROUTINE fp_grac(mode,f,string)
       
    IMPLICIT NONE
    CHARACTER(LEN=*),intent(in):: string
    REAL(rkind),INTENT(IN),DIMENSION(:,:,:,:):: f
    INTEGER,INTENT(IN):: mode
    REAL(rkind),ALLOCATABLE:: f_all(:,:,:,:)

    ALLOCATE(f_all(nthmax,npmax,nrmax,nsamax))

    CALL fp_gather_f(mode,f,f_all)
    
    IF(nrank.EQ.0) CALL fp_gracx(mode,f_all,string)
    
    RETURN
  END SUBROUTINE fp_grac
  
  SUBROUTINE fp_gather_f(mode,f,f_all)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    REAL(rkind),INTENT(IN):: f(:,:,:,:)
    REAL(rkind),INTENT(OUT):: f_all(nthmax,npmax,nrmax,nsamax)
    REAL(rkind),ALLOCATABLE:: f_temp(:,:,:,:)
    REAL(rkind),ALLOCATABLE:: dsend(:,:,:,:),drecv(:,:,:,:)
    INTEGER:: np,nth,nr,nsa,nsend

    ALLOCATE(f_temp(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend))

    SELECT CASE(mode)
    CASE(0,4) ! SPP
       CALL fp_convert_f0(f,f_temp)
    CASE(1,5) ! DPP,DPT,FPP
       CALL fp_convert_f1(f,f_temp)
    CASE(2,6) ! DTP,DTT,FTH
       CALL fp_convert_f2(f,f_temp)
    CASE(3,7) ! DRR,FRR
       CALL fp_convert_f3(f,f_temp)
    END SELECT

    nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)*(nsaend-nsastart+1)
    ALLOCATE(dsend(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend))
    ALLOCATE(drecv(nthmax,npmax,nrmax,nsamax))

    DO nsa=nsastart,nsaend
       DO nr=nrstart,nrend
          DO np=npstart,npend
             DO nth=1,nthmax
                dsend(nth,np,nr,nsa)=f_temp(nth,np,nr,nsa)
             END DO
          END DO
       END DO
    END DO
    CALL mtx_gather_real8(dsend,nsend,drecv) 
    IF(nrank.EQ.0) THEN
       DO nsa=1,nsamax
          DO nr=1,nrmax
             DO np=1,npmax
                DO nth=1,nthmax
                   f_all(nth,np,nr,nsa)=drecv(nth,np,nr,nsa)
                END DO
             END DO
          END DO
       END DO
    END IF
    CALL mtx_reset_communicator
    RETURN
  END SUBROUTINE fp_gather_f

  SUBROUTINE fp_convert_f0(f,f_temp)

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: &
         f(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend)
    REAL(rkind),INTENT(OUT):: &
         f_temp(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend)
    INTEGER:: nth,np,nr,nsa

    DO nsa=nsastart,nsaend
       DO nr=nrstart,nrend
          DO np=npstart,npend
             DO nth=1,nthmax
                f_temp(nth,np,nr,nsa)=f(nth,np,nr,nsa)
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE fp_convert_f0

  SUBROUTINE fp_convert_f1(f,f_temp)

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: &
         f(nthmax,npstart:npendwg,nrstart:nrendwm,nsastart:nsaend)
    REAL(rkind),INTENT(OUT):: &
         f_temp(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend)
    INTEGER:: nth,np,nr,nsa

    DO nsa=nsastart,nsaend
       DO nr=nrstart,nrend
          DO np=npstart,npend
             DO nth=1,nthmax
                f_temp(nth,np,nr,nsa) &
                     =0.5D0*(f(nth,np,nr,nsa)+f(nth,np+1,nr,nsa))
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE fp_convert_f1

  SUBROUTINE fp_convert_f2(f,f_temp)

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: &
         f(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax)
    REAL(rkind),INTENT(OUT):: &
         f_temp(nthmax,npstart:npend,nrstart:nrend,nsamax)
    INTEGER:: nth,np,nr,nsa

    DO nsa=nsastart,nsaend
       DO nr=nrstart,nrend
          DO np=npstart,npend
             DO nth=1,nthmax
                f_temp(nth,np,nr,nsa) &
                     =0.5D0*(f(nth,np,nr,nsa)+f(nth+1,np,nr,nsa))
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE fp_convert_f2

  SUBROUTINE fp_convert_f3(f,f_temp)

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: &
         f(nthmax,npstart:npend,nrstart:nrendwg,nsastart:nsaend)
    REAL(rkind),INTENT(OUT):: &
         f_temp(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend)
    INTEGER:: nth,np,nr,nsa

    DO nsa=nsastart,nsaend
       DO nr=nrstart,nrend
          DO np=npstart,npend
             DO nth=1,nthmax
                f_temp(nth,np,nr,nsa) &
                     =0.5D0*(f(nth,np,nr,nsa)+f(nth,np,nr+1,nsa))
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE fp_convert_f3

! --- ask and plot contour ---
  
  SUBROUTINE fp_gracx(mode,f,string)
!       
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    REAL(rkind),INTENT(IN):: f(nthmax,npmax,nrmax,nsamax)
    CHARACTER(LEN=*),INTENT(IN):: string
    REAL,ALLOCATABLE:: gf(:,:)
    CHARACTER(LEN=80):: string1
    INTEGER:: nr,nsa,np,nth

    ALLOCATE(gf(npmax,nthmax))

    IF(NGRAPH.EQ.0) THEN
       CALL fp_gtof(mode,f,string)
       RETURN
    ENDIF

    nr=1
    nsa=1
1   CONTINUE
    IF(nrmax.GT.1.OR.nsamax.GT.1) THEN
       WRITE(6,'(A,I3,A,I2,A)') &
            '# INPUT NR (1..',nrmax,'),NSA(1..',nsamax,'): 0 for end ?'
       READ(5,*,ERR=1,END=9000) nr,nsa
       IF(nr.LE.0) GOTO 9000
       IF(nr.GT.nrmax) GOTO 1
       IF(nsa.LE.0) GOTO 9000
       IF(nsa.GT.nsamax) GOTO 1
    END IF
    WRITE(string1,'(A,A,I3,A,I2)') string,' : NR=',nr,',NSA=',nsa

    SELECT CASE(mode)
    CASE(0:3)
       DO nth=1,nthmax
          DO np=1,npmax
             gf(np,nth)=gdclip(f(nth,np,nr,nsa))
          END DO
       END DO
    CASE(4:7)
       DO nth=1,nthmax
          DO np=1,npmax
             IF(f(nth,np,nr,nsa).LT.1.D-70) THEN
                gf(np,nth)=-70.0
             ELSE
                gf(np,nth)=gdclip(LOG10(ABS(f(nth,np,nr,nsa))))
             ENDIF
          END DO
       END DO
    END SELECT
!
    CALL fp_gracxx(mode,nsa,gf,string1)

    IF(nrmax.GT.1) GOTO 1
!
9000 RETURN
  END SUBROUTINE fp_gracx

  SUBROUTINE fp_gracxx(mode,nsa,gf,string)
!       
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode,nsa
    REAL,INTENT(IN):: gf(npmax,nthmax)
    CHARACTER(LEN=*),INTENT(IN):: string
    REAL:: gp(npmax),gth(nthmax)
    INTEGER,PARAMETER:: nglm=30
    REAL:: zl(nglm),rgb(3,nglm),wln(nglm)
    INTEGER:: iln(nglm)
    REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
    REAL:: gpmax,gfmin,gfmax,gfmin1,gfmax1,gfstep,gpmin1,gpmax1,gpstep
    REAL:: GLIN,GFFMAX
    INTEGER:: lmode,ngl,nglmax
    INTEGER:: ns,np,nth

    lmode=mode/4
    ns=ns_nsa(nsa)

    DO np=1,npmax
       gp(np)=gdclip(pg(np,ns))
    END DO
    DO nth=1,nthmax
       gth(nth)=gdclip(thm(nth))
    END DO

    gpmax=gdclip(pmax(ns))

    CALL PAGES
    CALL SETLNW(0.07)
    CALL SETCHS(.3,0.)
    CALL SETFNT(32)

    CALL GMNMX2(gf,npmax,1,npmax,1,1,nthmax,1,gfmin,gfmax)
    CALL GQSCAL(gfmin,gfmax,gfmin1,gfmax1,gfstep)
    CALL GQSCAL(0.0,gpmax,gpmin1,gpmax1,gpstep)

    CALL GDEFIN(3.,23.,2.,12.,-gpmax,gpmax,0.,gpmax)
    CALL GFRAME
    CALL SETLNW(0.035)
    CALL GSCALE(0.,gpstep,0.,gpstep,1.0,0)
    CALL SETLNW(0.07)
    CALL GVALUE(0.,gpstep*2,0.,gpstep*2,ngslen(2*gpstep))

    IF(LMODE.EQ.0) THEN
       IF(GFMIN*GFMAX.GE.0.0) THEN
          IF(GFMIN.GE.0.0) THEN
             NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
             DO NGL=1,NGLMAX
                ZL(NGL)=GFMIN1+0.5*GFSTEP*(NGL-1)
                RGB(1,NGL)=1.D0
                RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                ILN(NGL)=0
                WLN(NGL)=0.07
             ENDDO
             CALL CONTG4X(gf,gp,gth,npmax,npmax,nthmax,zl,rgb,iln,wln,nglmax,0)
          ELSE
             NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
             DO NGL=1,NGLMAX
                ZL(NGL)=GFMAX1-0.5*GFSTEP*(NGL-1)
                RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                RGB(3,NGL)=1.D0
                ILN(NGL)=0
                WLN(NGL)=0.07
             ENDDO
             CALL CONTG4X(gf,gp,gth,npmax,npmax,nthmax,zl,rgb,iln,wln,nglmax,0)
          ENDIF
       ELSE
          GFFMAX=MAX(ABS(GFMAX1),ABS(GFMIN1))
          NGLMAX=INT(GFFMAX/(0.5*GFSTEP))
          DO NGL=1,NGLMAX
             ZL(NGL)=-0.25*GFSTEP-0.5*GFSTEP*(NGL-1)
             RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
             RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
             RGB(3,NGL)=1.D0
             ILN(NGL)=0
             WLN(NGL)=0.07
          ENDDO
          CALL CONTG4X(gf,gp,gth,npmax,npmax,nthmax,zl,rgb,iln,wln,nglmax,0)
          DO NGL=1,NGLMAX
             ZL(NGL)= 0.25*GFSTEP+0.5*GFSTEP*(NGL-1)
             RGB(1,NGL)=1.D0
             RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
             RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
             ILN(NGL)=0
             WLN(NGL)=0.07
          ENDDO
          CALL CONTG4X(gf,gp,gth,npmax,npmax,nthmax,zl,rgb,iln,wln,nglmax,0)
       ENDIF
    ELSE
       DO NGL=1,NGLINE
          ZL(NGL)=GFMAX-0.020*(NGL-1)**2
          CALL fp_setrgb(1.0-FLOAT(NGL-1)/FLOAT(NGLINE-1),RGB(1,NGL))
          ILN(NGL)=0
          WLN(NGL)=0.07
       ENDDO
       CALL CONTG4X(gf,gp,gth,npmax,npmax,nthmax,zl,rgb,iln,wln,ngline,0)
    ENDIF

    CALL SETLIN(0,2,7)
    CALL MOVE(24.0,1.0)
    CALL TEXT('PPARA',5)
    CALL MOVE(1.0,13.5)
    CALL TEXT('PPERP',5)

    CALL MOVE(3.0,12.5)
    CALL TEXT(STRING,LEN(STRING))
    CALL MOVE(8.0,12.5)
    CALL TEXT('FMIN =',6)
    CALL NUMBR(GFMIN,'(1PE12.4)',12)
    CALL MOVE(13.0,12.5)
    CALL TEXT('FMAX =',6)
    CALL NUMBR(GFMAX,'(1PE12.4)',12)
    IF(LMODE.EQ.0) THEN
       CALL MOVE(18.0,12.5)
       CALL TEXT('STEP =',6)
       CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
    ENDIF
    CALL PAGEE

    RETURN
  END SUBROUTINE fp_gracxx

  SUBROUTINE fp_setrgb(F,RGB)
    IMPLICIT NONE
    REAL,INTENT(IN):: F
    REAL,DIMENSION(3),INTENT(OUT):: RGB
    INTEGER,PARAMETER:: NFMAX=8
    REAL,DIMENSION(3,NFMAX):: RGBC
    DATA RGBC/ 0.0,0.0,0.0, &
               0.0,0.0,1.0, &
               0.0,0.8,1.0, &
               0.0,0.8,0.0, &
               1.0,0.8,0.0, &
               1.0,0.4,0.0, &
               1.0,0.0,0.0, &
               1.0,1.0,1.0/
    REAL(rkind):: GF,DF
    INTEGER:: IM

    GF=F*DBLE(NFMAX-1)+1
    IM=MIN(INT(GF),NFMAX-1)
    DF=GF-IM
    RGB(1)=RGBC(1,IM)*(1.D0-DF)+RGBC(1,IM+1)*DF
    RGB(2)=RGBC(2,IM)*(1.D0-DF)+RGBC(2,IM+1)*DF
    RGB(3)=RGBC(3,IM)*(1.D0-DF)+RGBC(3,IM+1)*DF
    RETURN
  END SUBROUTINE fp_setrgb

  SUBROUTINE fp_gtof(mode,f,string)

    INTEGER,INTENT(IN):: mode
    REAL(rkind),INTENT(IN):: f(nthmax,npmax,nrmax,nsamax)
    CHARACTER(LEN=*),INTENT(IN):: string
    real(rkind):: fl(npmax,nthmax)
    INTEGER:: nr,nsa,nth,np

    nr=1
    nsa=1

1   CONTINUE
    IF(nrmax.GT.1.OR.nsamax.GT.1) THEN
       WRITE(6,'(A,I3,A,I2,A)') &
            '# INPUT NR (1..',nrmax,'),NSA(1..',nsamax,'): 0 for end ?'
       READ(5,*,ERR=1,END=9000) nr,nsa
       IF(nr.LE.0) GOTO 9000
       IF(nr.GT.nrmax) GOTO 1
       IF(nsa.LE.0) GOTO 9000
       IF(nsa.GT.nsamax) GOTO 1
    END IF

    DO nth=1,nthmax
       DO np=1,npmax
          fl(np,nth)=f(nth,np,nr,nsa)
       ENDDO
    ENDDO

    CALL fp_foutx4(string,npmax,npmax,nthmax,fl)
    IF(nrmax.GT.1) GOTO 1

9000 RETURN
  END SUBROUTINE fp_gtof

  SUBROUTINE fp_foutx4(string,n1m,n1max,n2max,f)

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: string
    INTEGER,INTENT(IN):: n1m,n1max,n2max
    REAL(rkind),INTENT(IN):: f(n1m,n2max)
    INTEGER:: n1,n2

    WRITE(22,'(A)') string
    WRITE(22,'(2I10)') n1max,n2max
    DO n2=1,n2max
       WRITE(22,'(1P5E15.7)') (f(n1,n2),n1=1,n1max)
    ENDDO
    RETURN
  END SUBROUTINE fp_foutx4
END MODULE fpgsub
