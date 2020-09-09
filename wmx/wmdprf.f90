! wmdprf.f90

MODULE wmdprf

  PRIVATE
  PUBLIC wm_dprf

CONTAINS

!     ***** READ DIII-D PROFILE DATA *****

  SUBROUTINE wm_dprf(IERR)

!     DATA FORMAT (4E15.8)
!       RHO, NE, TE, TI

!       number of data is number of lines

    USE wmcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind),ALLOCATABLE:: &
         RHOFL(:),RNEFL(:),RTEFL(:),RTIFL(:),DERIV(:), &
         URNEFL(:,:),URTEFL(:,:),URTIFL(:,:)
    INTEGER:: NR,IST,NPFDSK,NS,NRFLMAX
    REAL(rkind):: D1,D2,D3,D4,PNI,XRHOL,RNEL,RTEL,RTIL

    ALLOCATE(RHOFL(nrmax+1),RNEFL(nrmax+1),RTEFL(nrmax+1),RTIFL(nrmax+1))
    ALLOCATE(URNEFL(4,nrmax+1),URTEFL(4,nrmax+1),URTIFL(4,nrmax+1))
    ALLOCATE(DERIV(nrmax+1))

    NPFDSK=22
    CALL FROPEN(NPFDSK,KNAMPF,1,0,'PF',IERR)
    IF(IERR.NE.0) RETURN

    REWIND(NPFDSK)
    NR=0
40  READ(NPFDSK,'(4E15.8)',IOSTAT=IST,ERR=80,END=90) D1,D2,D3,D4
    NR=NR+1
    RHOFL(NR)=D1
    RNEFL(NR)=D2
    RTEFL(NR)=D3
    RTIFL(NR)=D4
    GOTO 40

80  WRITE(6,*) 'XX WMDPRF: FILE READ ERROR: IOSTAT=',IST
    IERR=8003
    GOTO 9000

90  REWIND(NPFDSK)
    CLOSE(NPFDSK)
    NRFLMAX=NR

!     ----- Convert density in 10^19 to 10^20 -----

    DO NR=1,NRFLMAX
       RNEFL(NR)=RNEFL(NR)*0.1D0
    ENDDO

!     ----  Set coefficient for spline

    DERIV(1)=0
    CALL SPL1D(RHOFL,RNEFL,DERIV,URNEFL,NRFLMAX,1,IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX WMDPRF: SPL1D: RNEFL: IERR=',IERR
       GO TO 9000
    END IF
    DERIV(1)=0
    CALL SPL1D(RHOFL,RTEFL,DERIV,URTEFL,NRFLMAX,1,IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX WMDPRF: SPL1D: RTEFL: IERR=',IERR
       GO TO 9000
    END IF
    DERIV(1)=0
    CALL SPL1D(RHOFL,RTIFL,DERIV,URTIFL,NRFLMAX,1,IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX WMDPRF: SPL1D: RTIFL: IERR=',IERR
       GO TO 9000
    END IF

!     ----  Set profile data at the point calculated in wm-code.

    PNI=0.D0
    DO NS=2,NSMAX
       PNI=PNI+PZ(NS)*PN(NS)
    ENDDO
    DO NR=1,NRMAX+1
       XRHOL=XRHO(NR)
       IF (XRHOL.GT.1.0D0) THEN
          RNPRF(NR,1)=RNEFL(NRFLMAX)
          RTPRF(NR,1)=RTEFL(NRFLMAX)
          DO NS=2,NSMAX
             RNPRF(NR,NS)=RNEFL(NRFLMAX)*PN(NS)/PNI
             RTPRF(NR,NS)=RTIFL(NRFLMAX)
          ENDDO
       ELSE
          CALL SPL1DF(XRHOL,RNEL,RHOFL,URNEFL,NRFLMAX,IERR)
          IF(IERR.NE.0) THEN
             WRITE(6,*) 'XX WMDPRF: SPL1DF: RNEFL: IERR=',IERR
             GOTO 9000
          ENDIF
          RNPRF(NR,1)=RNEL
          DO NS=2,NSMAX
             RNPRF(NR,NS)=RNEL*PN(NS)/PNI
          ENDDO
          CALL SPL1DF(XRHOL,RTEL,RHOFL,URTEFL,NRFLMAX,IERR)
          IF(IERR.NE.0) THEN
             WRITE(6,*) 'XX WMDPRF: SPL1DF: RTEFL: IERR=',IERR
             GO TO 9000
          ENDIF
          RTPRF(NR,1)=RTEL
          CALL SPL1DF(XRHOL,RTIL,RHOFL,URTIFL,NRFLMAX,IERR)
          IF(IERR.NE.0) THEN
             WRITE(6,*) 'XX WMDPRF: SPL1DF: RTIFL: IERR=',IERR
             GO TO 9000
          ENDIF
          DO NS=2,NSMAX
             RTPRF(NR,NS)=RTIL
          ENDDO
          !            WRITE(6,'(I5,1P4E12.4)') NR,XRHOL,RNEL,RTEL,RTIL
       ENDIF
    ENDDO

!     ----  Change the value at center and surface

    DO NS=1,NSMAX
       PN(NS)  =RNPRF(1,NS)
       PTPR(NS)=RTPRF(1,NS)
       PTPP(NS)=RTPRF(1,NS)
       PNS(NS) =RNPRF(NRMAX+1,NS)
       PTS(NS) =RTPRF(NRMAX+1,NS)
    ENDDO

9000 RETURN
  END SUBROUTINE wm_dprf
END MODULE wmdprf
