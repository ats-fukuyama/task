C     $Id$
C
C     ***************************************************************
C
C        Save transport data
C
C     ***************************************************************
C
      SUBROUTINE TXSAVE
C
      INCLUDE 'txcomm.h'
C
      CHARACTER TXFNAM*100, STR*1, RCSId*100
      LOGICAL LEX
C
      RCSId = ' '
C
   10 CONTINUE
      WRITE(6,*) '# INPUT : SAVE FILE NAME'
      CALL GUFLSH
      READ(*,'(A100)',END=50) TXFNAM
      INQUIRE(FILE=TXFNAM,EXIST=LEX)
      IF (LEX) THEN
         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ', 
     &              'ARE YOU SURE {Y/N} ?'
         CALL GUFLSH
         READ(*,'(A1)') STR
         IF (STR.NE.'Y' .AND. STR.NE.'y') GOTO 10
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',ERR=20, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR OUTPUT.'
         GOTO 40
   20    CONTINUE
         WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ELSE
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',ERR=30, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (', TXFNAM, ') IS CREATED FOR OUTPUT.'
         GOTO 40
   30    CONTINUE
         WRITE(6,*) 'XX  NEW FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ENDIF
C
   40 CONTINUE
      WRITE(21) SLID
      WRITE(21) RCSId
C
      WRITE(21) RA,RB,RR,BB
      WRITE(21) PA,PZ,Zeff
      WRITE(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
      WRITE(21) De0,Di0,rMue0,rMui0,WPM0
      WRITE(21) Chie0,Chii0
      WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
      WRITE(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
      WRITE(21) rLn,rLT
      WRITE(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
      WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
      WRITE(21) DLT,DT,EPS
      WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
      WRITE(21) rG1
      WRITE(21) rIPs,rIPe
C
      WRITE(21) DR,TIME,TMAX,NT,NQMAX,IERR
      WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)
C
      WRITE(21) NGT,NGYTM,NGYVM
      WRITE(21) (GTX(I), I=0, NGT)
      WRITE(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
      WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
      CLOSE(21)
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'
C
   50 CONTINUE
C
      RETURN
      END
C
C     ***************************************************************
C
C        Load transport data
C
C     ***************************************************************
C
      SUBROUTINE TXLOAD
C
      INCLUDE 'txcomm.h'
C
      CHARACTER TXFNAM*100, LOADSLID*8, RCSId*100
      LOGICAL LEX
C
C tmp : NGYT
C
   10 CONTINUE
      WRITE(6,*) '# INPUT : LOAD FILE NAME'
      CALL GUFLSH
      READ(*,'(A100)',END=50) TXFNAM
      INQUIRE(FILE=TXFNAM,EXIST=LEX)
      IF (LEX) THEN
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',ERR=20, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR INPUT.'
         GOTO 30
   20    CONTINUE
         WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ELSE
         WRITE(6,*) 'XX  FILE (', TXFNAM, ') DOES NOT EXIST !'
         GOTO 10
      ENDIF
C
   30 CONTINUE
      READ(21,ERR=80) LOADSLID
C      IF(LOADSLID(1:5).EQ.'tx210') THEN
         READ(21) RCSId
C
         READ(21) RA,RB,RR,BB
         READ(21) PA,PZ,Zeff
         READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
         READ(21) De0,Di0,rMue0,rMui0,WPM0
         READ(21) Chie0,Chii0
         READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
         READ(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
         READ(21) rLn,rLT
         READ(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
         READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
         READ(21) DLT,DT,EPS
         READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
         READ(21) rG1
         READ(21) rIPs,rIPe
C
         READ(21) DR,TIME,TMAX,NT,NQMAX,IERR
         READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)
C
         READ(21) NGT,NGYT,NGYV
         READ(21) (GTX(I), I=0, NGT)
         READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
         READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
C      ENDIF
      CLOSE(21)
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
      IF (NQMAX .EQ. 16) THEN
         DO NR = 0, NRMAX - 1
            X(LQm5,NR) = BB
            X(LQn2,NR) = 0
         ENDDO
         DO NR = 0, NRMAX
            X(LQm2,NR) = 0
            X(LQb3,NR) = 0
         ENDDO
      ENDIF
C
      NGT=-1
      NGR=-1
      NGVV=-1
      rIP=rIPs
C
      CALL TXCALM
      CALL TXPRFG
      CALL TXCALV(X)
      CALL TXCALC
      CALL TXCALA
      CALL TXGLOB
      CALL TXSTGT(SNGL(TIME))
      CALL TXSTGV(SNGL(TIME))
      CALL TXSTGR
      CALL TXWDAT
      CALL TXWDAT2
C
   50 CONTINUE
      RETURN
C
   80 WRITE(6,*) 'XX READ ERROR in TXLOAD !'
      CLOSE(21)
      RETURN
      END
C
C     ***************************************************************
C
C        Save graphic data
C
C     ***************************************************************
C
      SUBROUTINE TXGSAV
C
      INCLUDE 'txcomm.h'
C
      CHARACTER TXFNAM*100, STR*1, RCSId*100
      LOGICAL LEX
C
      RCSId = ' '
C
   10 CONTINUE
      WRITE(6,*) '# INPUT : SAVE FILE NAME'
      CALL GUFLSH
      READ(*,'(A100)',END=50) TXFNAM
      INQUIRE(FILE=TXFNAM,EXIST=LEX)
      IF (LEX) THEN
         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ', 
     &              'ARE YOU SURE {Y/N} ?'
         CALL GUFLSH
         READ(*,'(A1)') STR
         IF (STR.NE.'Y' .AND. STR.NE.'y') GOTO 10
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',ERR=20, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR OUTPUT.'
         GOTO 40
   20    CONTINUE
         WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ELSE
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',ERR=30, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (', TXFNAM, ') IS CREATED FOR OUTPUT.'
         GOTO 40
   30    CONTINUE
         WRITE(6,*) 'XX  NEW FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ENDIF
C
   40 CONTINUE
      WRITE(21) SLID
      WRITE(21) RCSId
C
      WRITE(21) RA,RB,RR,BB
      WRITE(21) PA,PZ,Zeff
      WRITE(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
      WRITE(21) De0,Di0,rMue0,rMui0,WPM0
      WRITE(21) Chie0,Chii0
      WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
      WRITE(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
      WRITE(21) rLn,rLT
      WRITE(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
      WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
      WRITE(21) DLT,DT,EPS
      WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
      WRITE(21) rG1
      WRITE(21) rIPs,rIPe
C
      WRITE(21) DR,TIME,TMAX,NT,NQMAX,IERR
      WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)
C
      WRITE(21) NGR,NGYRM
      WRITE(21) NGT,NGYTM
      WRITE(21) NGVV,NGYVM
      WRITE(21) (GT(IGR), IGR=0, NGR)
      WRITE(21)(((GY(I,IGR,IGYR), I=0,NRMAX), IGR=0,NGR), IGYR=1,NGYRM)
      WRITE(21) (GTX(I), I=0, NGT)
      WRITE(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
      WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
      CLOSE(21)
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'
C
   50 CONTINUE
C
      RETURN
      END
C
C     ***************************************************************
C
C        Load graphic data
C
C     ***************************************************************
C
      SUBROUTINE TXGLOD
C
      INCLUDE 'txcomm.h'
C
      CHARACTER TXFNAM*100, LOADSLID*8, RCSId*100
      LOGICAL LEX
C
C tmp : NGYT
C
   10 CONTINUE
      WRITE(6,*) '# INPUT : LOAD FILE NAME'
      CALL GUFLSH
      READ(*,'(A100)',END=50) TXFNAM
      INQUIRE(FILE=TXFNAM,EXIST=LEX)
      IF (LEX) THEN
         OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',ERR=20, 
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR INPUT.'
         GOTO 30
   20    CONTINUE
         WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
         GOTO 10
      ELSE
         WRITE(6,*) 'XX  FILE (', TXFNAM, ') DOES NOT EXIST !'
         GOTO 10
      ENDIF
C
   30 CONTINUE
      READ(21,ERR=80) LOADSLID
C      IF(LOADSLID(1:5).EQ.'tx210') THEN
         READ(21) RCSId
C
         READ(21) RA,RB,RR,BB
         READ(21) PA,PZ,Zeff
         READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
         READ(21) De0,Di0,rMue0,rMui0,WPM0
         READ(21) Chie0,Chii0
         READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
         READ(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
         READ(21) rLn,rLT
         READ(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
         READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
         READ(21) DLT,DT,EPS
         READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
         READ(21) rG1
         READ(21) rIPs,rIPe
C
         READ(21) DR,TIME,TMAX,NT,NQMAX,IERR
         READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)
C
         READ(21) NGR,NGYR
         READ(21) NGT,NGYT
         READ(21) NGVV,NGYV
         READ(21) (GT(IGR), IGR=0, NGR)
         READ(21) (((GY(I,IGR,IGYR), I=0, NRMAX), IGR=0, NGR),
     &                               IGYR=1, NGYR)
         READ(21) (GTX(I), I=0, NGT)
         READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
         READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
C      ENDIF
      CLOSE(21)
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
      IF (NQMAX .EQ. 16) THEN
         DO NR = 0, NRMAX - 1
            X(LQm5,NR) = BB
            X(LQn2,NR) = 0
         ENDDO
         DO NR = 0, NRMAX
            X(LQm2,NR) = 0
            X(LQb3,NR) = 0
         ENDDO
      ENDIF
C
      CALL TXCALM
      CALL TXPRFG
      CALL TXCALV(X)
      CALL TXCALC
      CALL TXGLOB
      CALL TXWDAT
      CALL TXWDAT2
C
   50 CONTINUE
      RETURN
C
   80 WRITE(6,*) 'XX READ ERROR in TXGLOD !'
      CLOSE(21)
      RETURN
      END
