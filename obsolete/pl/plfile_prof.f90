MODULE plfile_prof_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC plfile_prof_read,plfile_prof_ntu,plfile_prof_q

  INTEGER(4):: nsmaxl,nrmaxl,nrmaxlq
  REAL(8),DIMENSION(:),POINTER:: prfr,prfrq
  REAL(8),DIMENSION(:,:),POINTER:: prfn,prft,prfu
  REAL(8),DIMENSION(:),POINTER:: prfq
  REAL(8),DIMENSION(:,:,:),POINTER:: uprfn,uprft,uprfu
  REAL(8),DIMENSION(:,:),POINTER:: uprfq

  CONTAINS

  SUBROUTINE plfile_prof_read(modeln,modelq,ierr)
    INTEGER(4),INTENT(IN):: modeln,modelq
    INTEGER(4),INTENT(OUT):: ierr
    CHARACTER(80):: knamfln,knamflt,knamflq
    REAL(8),DIMENSION(:,:),POINTER:: data
    REAL(8),DIMENSION(:),POINTER:: deriv
    INTEGER(4):: ndata,nrdata,ns,nr

    knamfln= 'taebm/77788_ne_fit_44985.dat'
    knamflt= 'taebm/77788_Te_fit_44985.dat'
    knamflq= 'taebm/77788_eq_data_44985.dat'

    if(modeln.eq.6) then

       nsmaxl=2
       nrdata=201

       ndata=2
       CALL plfile_prof_data(knamfln,ndata,data,nrdata,ierr)
       IF(ierr.NE.0) GOTO 900
       nrmaxl=nrdata
       IF(ASSOCIATED(prfr)) DEALLOCATE(prfr)
       ALLOCATE(prfr(nrmaxl))
       IF(ASSOCIATED(prfn)) DEALLOCATE(prfn)
       ALLOCATE(prfn(nrmaxl,nsmaxl))
       DO nr=1,nrmaxl
          prfr(nr)=data(nr,1)
          DO ns=1,nsmaxl
             prfn(nr,ns)=0.1D0*data(nr,2)
          ENDDO
       ENDDO
    
       ndata=2
       CALL plfile_prof_data(knamflt,ndata,data,nrdata,ierr)
       IF(ierr.NE.0) GOTO 900
       if(nrdata.ne.nrmaxl) goto 800
       IF(ASSOCIATED(prft)) DEALLOCATE(prft)
       ALLOCATE(prft(nrmaxl,nsmaxl))
       DO nr=1,nrmaxl
          DO ns=1,nsmaxl
             prft(nr,ns)=data(nr,2)
          ENDDO
       ENDDO
    
       IF(ASSOCIATED(prfu)) DEALLOCATE(prfu)
       ALLOCATE(prfu(nrmaxl,nsmaxl))
       DO ns=1,nsmaxl
          DO nr=1,nrmaxl
             prfu(nr,ns)=0.d0
          ENDDO
       ENDDO

       ALLOCATE(deriv(nrmaxl))
       IF(ASSOCIATED(uprfn)) DEALLOCATE(uprfn)
       IF(ASSOCIATED(uprft)) DEALLOCATE(uprft)
       IF(ASSOCIATED(uprfu)) DEALLOCATE(uprfu)
       ALLOCATE(uprfn(4,nrmaxl,nsmaxl))
       ALLOCATE(uprft(4,nrmaxl,nsmaxl))
       ALLOCATE(uprfu(4,nrmaxl,nsmaxl))

       DO ns=1,nsmaxl
          deriv(1)=0.D0
          CALL SPL1D(prfr,prfn(1:nrmaxl,ns),deriv, &
                     uprfn(1:4,1:nrmaxl,ns),nrmaxl,1,ierr)
       ENDDO
       DO ns=1,nsmaxl
          deriv(1)=0.D0
          CALL SPL1D(prfr,prft(1:nrmaxl,ns),deriv, &
                     uprft(1:4,1:nrmaxl,ns),nrmaxl,1,ierr)
       ENDDO
       DO ns=1,nsmaxl
          deriv(1)=0.D0
          CALL SPL1D(prfr,prfu(1:nrmaxl,ns),deriv, &
                     uprfu(1:4,1:nrmaxl,ns),nrmaxl,1,ierr)
       ENDDO
       DEALLOCATE(deriv)
    ENDIF

    IF(modelq.EQ.6) THEN

       ndata=9
       nrdata=257
       CALL plfile_prof_data(knamflq,ndata,data,nrdata,ierr)
       IF(ierr.NE.0) GOTO 900
       nrmaxlq=nrdata
       IF(ASSOCIATED(prfrq)) DEALLOCATE(prfrq)
       IF(ASSOCIATED(prfq)) DEALLOCATE(prfq)
       ALLOCATE(prfrq(nrmaxlq))
       ALLOCATE(prfq(nrmaxlq))
       DO nr=1,nrmaxlq
          prfrq(nr)=data(nr,1)
          prfq(nr)=data(nr,2)
       ENDDO
       ALLOCATE(deriv(nrmaxlq))
       IF(ASSOCIATED(uprfq)) DEALLOCATE(uprfq)
       ALLOCATE(uprfq(4,nrmaxlq))
       deriv(1)=0.D0
       CALL SPL1D(prfrq,prfq,deriv,uprfq,nrmaxlq,1,ierr)
       DEALLOCATE(deriv)
    ENDIF

!    DO nr=1,nrmaxl
!       WRITE(6,'(I5,1P5E12.4)') nr,prfr(nr),prfn(nr,1),prfn(nr,2), &
!            prft(nr,1),prft(nr,2)
!    ENDDO
!    DO nr=1,nrmaxlq
!       WRITE(6,'(I5,1P2E12.4)') nr,prfrq(nr),prfq(nr)
!    ENDDO


    RETURN

800 write(6,*) 'nrdata.NE.nrmaxl'
    ierr=20001
    return
900 return
  END SUBROUTINE plfile_prof_read

  SUBROUTINE plfile_prof_data(knamfl,ndata,data,nrdata,ierr)
    INTEGER(4),INTENT(IN):: ndata,nrdata
    CHARACTER(LEN=*),INTENT(IN):: knamfl
    REAL(8),DIMENSION(:,:),POINTER:: data
    INTEGER(4),INTENT(OUT):: ierr
    REAL(8):: rhol
    INTEGER(4):: nfl,modef,modep,i,ist,nr

    nfl=31
    modef=1 ! formatted
    modep=0 ! no prompt
    CALL FROPEN(nfl,knamfl,modef,modep,'plfile_prof_data',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX plfile_prof: fropen: irr=',ierr
       RETURN
    ENDIF

    IF(ASSOCIATED(data)) DEALLOCATE(data)
    ALLOCATE(data(nrdata,ndata))
    DO nr=1,nrdata
       READ(nfl,*,IOSTAT=ist,ERR=800,END=900) (data(nr,i),i=1,ndata)
    ENDDO
    CLOSE(nfl)
    RETURN
800 WRITE(6,*) 'XX plfile_prof: read error: IOSTAT=',ist
    ierr=10000+ist
    return
900 WRITE(6,*) 'XX plfile_prof: read end: IOSTAT=',ist
    ierr=10000+ist
    return
  END SUBROUTINE plfile_prof_data


!     ***** Interpolation of profile at a given point *****

  SUBROUTINE plfile_prof_ntu(rhol,ns,pnl,ptl,pul)

    REAL(8),INTENT(IN):: RHOL
    INTEGER(4),INTENT(IN):: NS
    REAL(8),INTENT(OUT)::  PNL,PTL,PUL 
    INTEGER(4):: IERR

    IF(NS.GT.NSMAXL) THEN
       WRITE(6,*) 'XX PLFILE_PROF: NS.GT.NSMAXL:',NS,NSMAXL
       STOP
    ENDIF
    IF (RHOL.LT.0.0D0) THEN
       PNL=prfn(1,NS)
       PTL=prft(1,NS)
       PUL=prfu(1,NS)
    ELSE IF(RHOL.GT.1.0D0) THEN
       PNL = 0.D0
       PTL = prft(nrmaxl,NS)
       PUL = 0.D0
    ELSE
       CALL SPL1DF(Rhol,PNL,prfr,uprfn(1:4,1:nrmaxl,ns),nrmaxl,IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX PLFILE_PROF: SPL1DF:pnl: IERR=',IERR
          STOP
       ENDIF
       CALL SPL1DF(Rhol,PTL,prfr,uprft(1:4,1:nrmaxl,ns),nrmaxl,IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX PLFILE_PROF: SPL1DF:ptl: IERR=',IERR
          STOP
       ENDIF
       CALL SPL1DF(Rhol,PUL,prfr,uprfu(1:4,1:nrmaxl,ns),nrmaxl,IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX PLFILE_PROF: SPL1DF:pul: IERR=',IERR
          STOP
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE plfile_prof_ntu

!     ***** Interpolation of profile at a given point *****

  SUBROUTINE plfile_prof_q(rhol,pql)

    REAL(8),INTENT(IN):: RHOL
    REAL(8),INTENT(OUT)::  PQL
    INTEGER(4):: IERR

    IF (rhol.LT.0.0D0) THEN
       pql=prfq(1)
    ELSE IF(rhol.GT.1.0D0) THEN
       pql=prfq(nrmaxlq)*rhol**2
    ELSE
       CALL SPL1DF(rhol,pql,prfrq,uprfq,nrmaxlq,IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX PLFILE_PROF_Q: SPL1DF:pql: IERR=',IERR
          STOP
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE plfile_prof_q
  
END MODULE plfile_prof_mod

