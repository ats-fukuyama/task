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
      INCLUDE 'wmcomm.h'
C
      DIMENSION RHOFL(NRM),RNEFL(NRM),RTEFL(NRM),RTIFL(NRM)
      DIMENSION DERIV(NRM),URNEFL(4,NRM),URTEFL(4,NRM),URTIFL(4,NRM)
      LOGICAL LEX
C
      IERR=0
      NPFDSK=22
      IF(KNAMPF(1:2).EQ.'  ') GOTO 9000
      INQUIRE(FILE=KNAMPF,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         OPEN(NPFDSK,FILE=KNAMPF,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='FORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMPF,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         IERR=8001
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX FILE (',KNAMPF,') NOT FOUND'
         IERR=8002
         GOTO 9000
      ENDIF
c
   30 NR=0
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
   90 NRFLMAX=NR
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
      DO NR=1,NRMAX+1
         XRHOL=XRHO(NR)
         IF (XRHOL.GT.1.0D0) THEN
            DO NS=1,NSMAX
               RNPRF(NR,NS)=0.0D0
               RTPRF(NR,NS)=0.0D0
            ENDDO
         ELSE
            CALL SPL1DF(XRHOL,RNEL,RHOFL,URNEFL,NRFLMAX,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX WMDPRF: SPL1DF: RNEFL: IERR=',IERR
               GOTO 9000
            ENDIF
            RNPRF(NR,1)=RNEL
            RNPRF(NR,2)=RNEL/PZ(2)
            DO NS=3,NSMAX
               RNPRF(NR,NS)=RNEL
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
            WRITE(6,'(I5,1P4E12.4)') NR,XRHOL,RNEL,RTEL,RTIL
         ENDIF
      ENDDO
C
C     ----  Change the value at center and surface
C
      PN(1)  =RNEFL(1)
      PTPR(1)=RTEFL(1)
      PTPP(1)=RTEFL(1)
      PNS(1) =RNEFL(NRFLMAX)
      PTS(1) =RTEFL(NRFLMAX)
C
      PN(2)  =RNEFL(1)/PZ(2)
      PTPR(2)=RTIFL(1)
      PTPP(2)=RTIFL(1)
      PNS(2) =RNEFL(NRFLMAX)/PZ(2)
      PTS(2) =RTIFL(NRFLMAX)
      DO NS=3,NSMAX
         PN(NS)  =0.D0
         PTPR(NS)=RTIFL(1)
         PTPP(NS)=RTIFL(1)
         PNS(NS) =0.D0
         PTS(NS) =RTIFL(NRFLMAX)
      ENDDO
C
 9000 RETURN
      END
