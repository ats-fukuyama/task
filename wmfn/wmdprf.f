C     $Id$
C
C     ***** READ DIII-D PROFILE DATA *****
C
      SUBROUTINE WMDPRF(IERR)
C
C     DATA FORMAT (4E15.8)
C       RHO, NE, TE, TI
C
C       number of data is number of lines
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RHOFL(NRM),RNEFL(NRM),RTEFL(NRM),RTIFL(NRM)
      DIMENSION DERIV(NRM),URNEFL(4,NRM),URTEFL(4,NRM),URTIFL(4,NRM)
C
      NPFDSK=22
      CALL FROPEN(NPFDSK,KNAMPF,1,0,'PF',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(NPFDSK)
      NR=0
   40 READ(NPFDSK,'(4E15.8)',IOSTAT=IST,ERR=80,END=90) D1,D2,D3,D4
      NR=NR+1
      RHOFL(NR)=D1
      RNEFL(NR)=D2
      RTEFL(NR)=D3
      RTIFL(NR)=D4
      GOTO 40
C
   80 WRITE(6,*) 'XX WMDPRF: FILE READ ERROR: IOSTAT=',IST
      IERR=8003
      GOTO 9000
C
   90 REWIND(NPFDSK)
      CLOSE(NPFDSK)
      NRFLMAX=NR
C
C     ----- Convert density in 10^19 to 10^20 -----
C
      DO NR=1,NRFLMAX
         RNEFL(NR)=RNEFL(NR)*0.1D0
      ENDDO
C
C     ----  Set coefficient for spline
C
      DERIV(1)=0
      CALL SPL1D(RHOFL,RNEFL,DERIV,URNEFL,NRFLMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMDPRF: SPL1D: RNEFL: IERR=',IERR
      IF(IERR.NE.0) GO TO 9000
      DERIV(1)=0
      CALL SPL1D(RHOFL,RTEFL,DERIV,URTEFL,NRFLMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMDPRF: SPL1D: RTEFL: IERR=',IERR
      IF(IERR.NE.0) GO TO 9000
      DERIV(1)=0
      CALL SPL1D(RHOFL,RTIFL,DERIV,URTIFL,NRFLMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMDPRF: SPL1D: RTIFL: IERR=',IERR
      IF(IERR.NE.0) GO TO 9000
C
C     ----  Set profile data at the point calculated in wm-code.
C
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
C            WRITE(6,'(I5,1P4E12.4)') NR,XRHOL,RNEL,RTEL,RTIL
         ENDIF
      ENDDO
C
C     ----  Change the value at center and surface
C
      DO NS=1,NSMAX
         PN(NS)  =RNPRF(1,NS)
         PTPR(NS)=RTPRF(1,NS)
         PTPP(NS)=RTPRF(1,NS)
         PNS(NS) =RNPRF(NRMAX+1,NS)
         PTS(NS) =RTPRF(NRMAX+1,NS)
      ENDDO
C
 9000 RETURN
      END
