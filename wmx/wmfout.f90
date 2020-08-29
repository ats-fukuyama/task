C     $Id: wmfout.f,v 1.6 2009/03/23 04:56:08 honda Exp $
C
C***********************************************************************
C
C     Save field data
C
C***********************************************************************
C
      SUBROUTINE WMSAVE
C
      INCLUDE './wmcomm.inc'
C
      CALL FWOPEN(21,KNAMWM,0,MODEFW,'WM',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WMSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF
C
C      WRITE(6,*) MDSIZ,NDSIZ,NRMAX
C      WRITE(6,'(4I5,1P2E12.4)') 
C     &           ((((NR,ND,MD,I,CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),
C     &                  ND=1,NDSIZ),NR=1,NRMAX)
C
      REWIND(21)
      WRITE(21) MDSIZ,NDSIZ,NRMAX,NSMAX
      WRITE(21) RA,RR,BB,CRF,NPH0,NTH0
      WRITE(21) ((((CEFLD(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((CBFLD(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((PABS(MD,ND,NR,NS),        MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX),NS=1,NSMAX)
      WRITE(21)  (((PCUR(MD,ND,NR),           MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((CEN(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((CEP(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)

      WRITE(21) CRADTT,PABSTT,PCURT
      WRITE(21) (PABST(NS),NS=1,NSMAX)
      WRITE(21) ((((CEFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((CBFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((PABSK(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX),NS=1,NSMAX)
      WRITE(21) NPHMAX,NHC,NHHMAX

      IF (NPHMAX > 1)THEN
      ENDIF
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
C***********************************************************************
C
C     Load field data
C
C***********************************************************************
C
      SUBROUTINE WMLOAD
C
      INCLUDE './wmcomm.inc'
C
      IERR=0
      MODEEG=0
      MODELK=1
C      MODELK=0
C
      CALL WMSETG(IERR)
      IF(IERR.NE.0) RETURN
      CALL DPCHEK(NTHMAX,NRMAX+1,XRHO(1),XRHO(NRMAX+1),RR,IERR)
      IF(IERR.NE.0) RETURN
      CALL WMSETJ(IERR)
      IF(IERR.NE.0) RETURN
      CALL WMSETEW

      CALL FROPEN(21,KNAMWM,0,0,'WM',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WMSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF
C
C      WRITE(6,*) MDSIZ,NDSIZ,NRMAX
C      WRITE(6,'(4I5,1P2E12.4)')
C     &           ((((NR,ND,MD,I,CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),
C     &                  ND=1,NDSIZ),NR=1,NRMAX)
C
      REWIND(21)
      READ(21) MDSIZ,NDSIZ,NRMAX,NSMAX
      READ(21) RA,RR,BB,CRF,NPH0,NTH0
      READ(21) ((((CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((CBFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((PABS(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX),NS=1,NSMAX)
      READ(21)  (((PCUR(MD,ND,NR),           MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((CEN(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((CEP(I,MD,ND,NR),I=1,3), MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)

      READ(21) CRADTT,PABSTT,PCURT
      READ(21) (PABST(NS),NS=1,NSMAX)
      READ(21) ((((CEFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((CBFLDK(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      READ(21) ((((PABSK(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX),NS=1,NSMAX)
      READ(21) NPHMAX,NHC,NHHMAX

      IF (NPHMAX > 1)THEN
      ENDIF
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED TO THE FILE.'
C
      CALL PREGOUT  
C
      RETURN
      END
C
C
C
C
      SUBROUTINE  PREGOUT
C
      INCLUDE './wmcomm.inc'
      DIMENSION DS(NRM),DSS(NTHM,NHHM,NRM)
C

      CW=2.D0*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NHHMAX)

      IF(NRANK.EQ.0) THEN
C
        DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
C            PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH*DPH*NHC
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
!         print *,NS,NR
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
C            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DPH
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DPH*NHC
!            PRINT *,PABS(NTH,NHH,NR,NS)*DTH*DPH
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C      PCURT=0.D0
C      DO NR=1,NRMAX
C         PCURT=PCURT+PCURR(NR)
C      ENDDO
C
C      DO NS=1,NSMAX
C         PABST(NS)=0.D0
C         DO NR=1,NRMAX
C            PABST(NS)=PABST(NS)+PABSR(NR,NS)
C         ENDDO
C      ENDDO
C
C      PABSTT=0.D0
C      DO NS=1,NSMAX
C         PABSTT=PABSTT+PABST(NS)
C      ENDDO
C
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
C            DS(NR)=DS(NR)+DSSS*DTH*DPH
            DS(NR)=DS(NR)+DSSS*DTH*DPH*NHC
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
      ENDIF

      END
