! wmfile.f90

MODULE wmfile

  PRIVATE
  PUBLIC wm_save,wm_load

CONTAINS

!***********************************************************************

!     Save field data

!***********************************************************************

  SUBROUTINE wm_save(ierr)

    USE wmcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: I,MD,ND,NR,NS
    COMPLEX(rkind):: CRF

    CRF=CMPLX(RF,RFI)

    CALL FWOPEN(21,KNAMWM,0,MODEFW,'WM',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX wm_save: FWOPEN: IERR=',IERR
       RETURN
    ENDIF

!      WRITE(6,*) MDSIZ,NDSIZ,NRMAX
!      WRITE(6,'(4I5,1P2E12.4)') 
!     &           ((((NR,ND,MD,I,CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),
!     &                  ND=1,NDSIZ),NR=1,NRMAX)

    REWIND(21)
    WRITE(21) MDSIZ,NDSIZ,NRMAX,NSMAX
    WRITE(21) RA,RR,BB,CRF,NPH0,NTH0
    WRITE(21) ((((CEFLD(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((CBFLD(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((PABS(MD,ND,NR,NS),        MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX),NS=1,NSMAX)
    WRITE(21)  (((PCUR(MD,ND,NR),           MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((CEN(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((CEP(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)

    WRITE(21) CPRAD,PABSTT,PCURT
    WRITE(21) (PABST(NS),NS=1,NSMAX)
    WRITE(21) ((((CEFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((CBFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    WRITE(21) ((((PABSK(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX),NS=1,NSMAX)
    WRITE(21) NPHMAX,NHC,NHHMAX

    CLOSE(21)

    WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

    RETURN
  END SUBROUTINE wm_save
  
!***********************************************************************

!     Load field data

!***********************************************************************

  SUBROUTINE wm_load(ierr)

    USE wmcomm
    USE wmsetg
    USE wmexec
    USE libfio
    USE dpparm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: I,MD,ND,NR,NS
    COMPLEX(rkind):: CRF

    IERR=0
    MODEEG=0
    MODELK=1

    CALL wm_setg(IERR)
    IF(IERR.NE.0) RETURN

    CALL dpprep_local(ierr)
    IF(IERR.NE.0) RETURN

    CALL wm_setj(IERR)
    IF(IERR.NE.0) RETURN

    CALL wm_setew

    CALL FROPEN(21,KNAMWM,0,0,'WM',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX WMSAVE: FWOPEN: IERR=',IERR
       RETURN
    ENDIF

    REWIND(21)
    READ(21) MDSIZ,NDSIZ,NRMAX,NSMAX
    READ(21) RA,RR,BB,CRF,NPH0,NTH0
    READ(21) ((((CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((CBFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((PABS(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX),NS=1,NSMAX)
    READ(21)  (((PCUR(MD,ND,NR),           MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((CEN(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((CEP(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)

    READ(21) CPRAD,PABSTT,PCURT
    READ(21) (PABST(NS),NS=1,NSMAX)
    READ(21) ((((CEFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((CBFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX)
    READ(21) ((((PABSK(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ), &
         NR=1,NRMAX),NS=1,NSMAX)
    READ(21) NPHMAX,NHC,NHHMAX

    CLOSE(21)

    WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED TO THE FILE.'

    CALL PREGOUT  

    RETURN
  END SUBROUTINE wm_load


  SUBROUTINE  PREGOUT

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: DS(:),DSS(:,:,:)
    INTEGER:: NR,NHH,NTH,NS,NHHF,NTHF
    REAL(rkind):: DTH,DHH,DPSIPDRHO,DSSS,DRHO
    COMPLEX(rkind):: CRF
    
    ALLOCATE(DS(nrmax+1),DSS(nthmax,nhhmax,nrmax+1))

    DTH=2.D0*PI/DBLE(NTHMAX)
    DHH=2.D0*PI/DBLE(NHHMAX*NHC)

    IF(NRANK.EQ.0) THEN

       DO NR=1,NRMAX
          PCURR(NR)=0.D0
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH*DHH
             ENDDO
          ENDDO
       ENDDO

       DO NS=1,NSMAX
          DO NR=1,NRMAX
             PABSR(NR,NS)=0.D0
             DO NHH=1,NHHMAX
                DO NTH=1,NTHMAX
                   PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DHH
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       PCURT=0.D0
       DO NR=1,NRMAX
          PCURT=PCURT+PCURR(NR)
       ENDDO

       DO NS=1,NSMAX
          PABST(NS)=0.D0
          DO NR=1,NRMAX
             PABST(NS)=PABST(NS)+PABSR(NR,NS)
          ENDDO
       ENDDO

       PABSTT=0.D0
       DO NS=1,NSMAX
          PABSTT=PABSTT+PABST(NS)
       ENDDO

       NR=1
          DS(NR)=0.D0
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                DSS(NTH,NHH,NR)=0.D0
             ENDDO
          ENDDO
       DO NR=2,NRMAX
          DS(NR)=0.D0
          DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                NHHF=(NHH-1)*FACTOR_NHH +1
                NTHF=(NTH-1)*FACTOR_NTH +1
                IF(MODELG.EQ.3) THEN
                   DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
                ELSE
                   DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
                ENDIF
                DSSS=DPSIPDRHO*DRHO*RJ(NTHF,NHHF,NR)
                DSS(NTH,NHH,NR)=1.D0/DSSS
                DS(NR)=DS(NR)+DSSS*DTH*DHH*NHC
             ENDDO
          ENDDO
          DS(NR)=1.D0/DS(NR)
       ENDDO
    ENDIF

  END SUBROUTINE PREGOUT
END MODULE wmfile
