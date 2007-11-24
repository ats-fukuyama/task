C     $Id$
C
C     ***** READ AND INTERPOLATE wmdata *****
C
      SUBROUTINE FPWMREAD(IERR)
C
      INCLUDE 'fpcomm.inc'
      COMPLEX*16 CEWL(NWTHM,NWRM),CEWX(NWTHM,NWRM),CEWY(NWTHM,NWRM)
      COMPLEX*16 CEWXY(NWTHM,NWRM)
C
      CALL FROPEN(21,KNAMWM,0,MODEFR,'WM',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWMREAD: FROPEN: IERR=',IERR
         RETURN
      ENDIF
C
      REWIND(21)
      READ(21) MDSIZ,NDSIZ,NRMAX1,NSMAX1
      WRITE(6,*) 'MDSIZ,NDSIZ,NRMAX1,NSMAX1=',MDSIZ,NDSIZ,NRMAX1,NSMAX1
      IF(MDSIZ.GT.NWTHM) THEN
         IERR=101
         RETURN
      ENDIF
      IF(NDSIZ.NE.1) THEN
         IERR=102
         RETURN
      ENDIF
      IF(NRMAX1.GT.NWRM) THEN
         IERR=103
         RETURN
      ENDIF
      READ(21) RAW,RRW,BBW,RFWR,RFWI,NPH0W,NTH0W
C
      READ(21) ((((CEWV(I,MD,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX1)
      CLOSE(21)
C
      NWTHMAX=MDSIZ
      NWRMAX=NRMAX1
      DO NWR=1,NWRMAX
         DO I=1,3
            DO NWTH=1,NWTHMAX
               CEWL(NWTH,1)=FACTWM*CEWV(I,NWTH,NWR)
            ENDDO
            CALL FPFFT(CEWL,NWTHMAX,1)
            DO NWTH=1,NWTHMAX
               CEWV(I,NWTH,NWR)=CEWL(NWTH,1)
            ENDDO
         ENDDO
      ENDDO
C
      DO NWR=1,NWRMAX
         RWSPL(NWR)=DBLE(NWR-1)/DBLE(NWRMAX-1)
      ENDDO
      DO NWTH=1,NWTHMAX+1
         THWSPL(NWTH)=DBLE(NWTH-1)*2.D0*PI/DBLE(NWTHMAX-1)
      ENDDO
      DO I=1,3
         DO NWR=1,NWRMAX
            DO NWTH=1,NWTHMAX
               CEWL(NWTH,NWR)=CEWV(I,NWTH,NWR)
            ENDDO
         ENDDO
         CALL CSPL2D(THWSPL,RWSPL,CEWL,CEWX,CEWY,CEWXY,UCEW(1,1,1,1,I),
     &               NWTHM,NWTHMAX+1,NWRMAX,4,0,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX FPWMGET: CSPL2DD: IERR=',IERR
            IERR=201
            RETURN
         ENDIF
      ENDDO
      RETURN
      END
C
C     ***** check INTERPOLATE wmdata *****
C
      SUBROUTINE FPWMCHEK
C
      INCLUDE 'fpcomm.inc'
      COMPLEX*16 CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
      DATA RL,THL/0.D0,0.D0/
C
 1010 CONTINUE
         WRITE(6,'(A,1P2E12.4)') 'Input RL,THL:',RL,THL
         READ(5,*,END=9000,ERR=1010) RL,THL
         IF(RL.EQ.0.D0) GOTO 9000
         CALL FPWMGET(RL,THL,CEWR1,CEWTH1,CEWPH1,
     &                       CKWR1,CKWTH1,CKWPH1,IERR)
         WRITE(6,'(2I10)') NPH0W,NTH0W
         WRITE(6,'(1P2E12.4)') RRW,RAW
         WRITE(6,'(1P6E12.4)') CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
         GOTO 1010
C
 9000    RETURN
         END
C
C     ***** INTERPOLATE wmdata *****
C
      SUBROUTINE FPWMGET(RL,THL,CEWR1,CEWTH1,CEWPH1,
     &                          CKWR1,CKWTH1,CKWPH1,IERR)
C
      INCLUDE 'fpcomm.inc'
      COMPLEX*16 CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
      COMPLEX*16 CEWDTH,CEWDR
C
      IERR=0
      CALL CSPL2DD(THL,RL,CEWR1,CEWDTH,CEWDR,THWSPL,RWSPL,
     &             UCEW(1,1,1,1,1),NWTHM,NWTHMAX+1,NWRMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWMGET: CSPL2DD: IERR=',IERR
         IERR=1
         RETURN
      ENDIF
      CKWR1=-CI*CEWDR/(CEWR1*RAW)
C
      CALL CSPL2DD(THL,RL,CEWTH1,CEWDTH,CEWDR,THWSPL,RWSPL,
     &             UCEW(1,1,1,1,2),NWTHM,NWTHMAX+1,NWRMAX,IERR)
      IF(IERR.NE.0) THEN
         IERR=2
         RETURN
      ENDIF
      IF(RL.LE.0.D0) THEN
         CKWTH1=0.D0
      ELSE
         CKWTH1=-CI*CEWDTH/(CEWTH1*RL*RAW)
      ENDIF
C
      CALL CSPL2DD(THL,RL,CEWPH1,CEWDTH,CEWDR,THWSPL,RWSPL,
     &             UCEW(1,1,1,1,3),NWTHM,NWTHMAX+1,NWRMAX,IERR)
      IF(IERR.NE.0) THEN
         IERR=3
         RETURN
      ENDIF
      CKWPH1=-NPH0W/RRW
C
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE FPFFT(CA,N,KEY)
C
      INCLUDE 'fpcomm.inc'
C
      COMPLEX*16 CA(N)
      DATA NS/0/
C
      IF(N.NE.1) THEN
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IND=1
            NS=N
         ENDIF
         IF(KEY.EQ.0) THEN
            CALL FFT2L(CA,CFFT,RFFT,LFFT,N,IND,KEY)
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CA(IX)=CFFT(I)
            ENDDO
         ELSE
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CFFT(I)=CA(IX)
            ENDDO
            CALL FFT2L(CFFT,CA,RFFT,LFFT,N,IND,KEY)
         ENDIF
      ENDIF
C
      RETURN
      END
