C
C     ***********************************************************
C
C           MODELLING TEST
C
C     ***********************************************************
C
      SUBROUTINE TRMDLT
C
      INCLUDE 'trcomm.inc'
      CHARACTER KFILE*80,KSTAT*3
      LOGICAL LEX
C
      IF(MDLUF.EQ.0) RETURN
C
C     *** RATIO AND MEAN SQUARE DEVIATION ***
C
      MODE=1
      NTL=1
      WSIM=0.D0
      WEXP=0.D0
      WSIME=0.D0
      WEXPE=0.D0
      WSIMI=0.D0
      WEXPI=0.D0
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(MODE.EQ.0) THEN
         NRHO09=INT(0.9*NRMAX)
         DO NR=1,NRHO09
            TESIML=RT(NR,1)
            TISIML=RT(NR,2)
            TEEXPL=RTU(NTL,NR,1)
            TIEXPL=RTU(NTL,NR,2)
            WSIM=WSIM+1.5D0*( RN(NR,1)*TESIML
     &                       +RN(NR,2)*TISIML)*DVRHO(NR)*DR
            WEXP=WEXP+1.5D0*( RNU_ORG(NTL,NR,1)*TEEXPL
     &                       +RNU_ORG(NTL,NR,2)*TIEXPL)*DVRHO(NR)*DR
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         NRHO09=INT(0.9*NRMAX)
         DO NR=1,NRHO09
            TESIML=RT(NR,1)-RT(NRHO09,1)
            TISIML=RT(NR,2)-RT(NRHO09,2)
            TEEXPL=RTU(NTL,NR,1)-RTU(NTL,NRHO09,1)
            TIEXPL=RTU(NTL,NR,2)-RTU(NTL,NRHO09,2)
            WSIM=WSIM+1.5D0*( RN(NR,1)*TESIML
     &                       +RN(NR,2)*TISIML)*DVRHO(NR)*DR
            WEXP=WEXP+1.5D0*( RNU_ORG(NTL,NR,1)*TEEXPL
     &                       +RNU_ORG(NTL,NR,2)*TIEXPL)*DVRHO(NR)*DR
CCC            WSIME=WSIME+1.5D0*RN(NR,1)*TESIML*DVRHO(NR)*DR
CCC            WSIMI=WSIMI+1.5D0*RN(NR,2)*TISIML*DVRHO(NR)*DR
CCC            WEXPE=WEXPE+1.5D0*RNU_ORG(NTL,NR,1)*TEEXPL*DVRHO(NR)*DR
CCC            WEXPI=WEXPI+1.5D0*RNU_ORG(NTL,NR,2)*TIEXPL*DVRHO(NR)*DR
         ENDDO
      ELSE
         DO NR=1,NRMAX
            TESIML=RT(NR,1)
            TISIML=RT(NR,2)
            TEEXPL=RTU(NTL,NR,1)
            TIEXPL=RTU(NTL,NR,2)
            WSIM=WSIM+1.5D0*( RN(NR,1)*TESIML
     &                       +RN(NR,2)*TISIML)*DVRHO(NR)*DR
            WEXP=WEXP+1.5D0*( RNU_ORG(NTL,NR,1)*TEEXPL
     &                       +RNU_ORG(NTL,NR,2)*TIEXPL)*DVRHO(NR)*DR
         ENDDO
      ENDIF
      WTO1=WSIM/WEXP-1.D0
      WTO2=(WSIM/WEXP-1.D0)**2
CCC      WTO2E=(WSIME/WEXPE-1.D0)**2
CCC      WTO2I=(WSIMI/WEXPI-1.D0)**2
C
C     *** STANDARD DEVIATION AND OFFSET ***
C
      RNUME=0.D0
      RNUMI=0.D0
      RDESE=0.D0
      RDESI=0.D0
      RDEOE=0.D0
      RDEOI=0.D0
      DO NR=1,NRMAX
         RP=DBLE(NR-1)*DR
         IF(RP.GE.0.2D0.AND.RP.LE.0.9D0) THEN
            RNUME=RNUME+RTU(NTL,NR,1)**2
            RNUMI=RNUMI+RTU(NTL,NR,2)**2
            RDESE=RDESE+(RT(NR,1)-RTU(NTL,NR,1))**2
            RDESI=RDESI+(RT(NR,2)-RTU(NTL,NR,2))**2
            RDEOE=RDEOE+(RT(NR,1)-RTU(NTL,NR,1))
            RDEOI=RDEOI+(RT(NR,2)-RTU(NTL,NR,2))
         ENDIF
      ENDDO
C
      NRHO02=INT(0.2*NRMAX)
      NRDMSH=NRHO09-NRHO02
      STDE=SQRT(RDESE)/SQRT(RNUME)
      STDI=SQRT(RDESI)/SQRT(RNUMI)
      OFFE=RDEOE      /SQRT(RNUME*DBLE(NRDMSH))
      OFFI=RDEOI      /SQRT(RNUMI*DBLE(NRDMSH))
C
C     *** WRITE OUT ***
C
      CALL KTRIM(KUFDEV,IUFDEV)
      CALL KTRIM(KUFDCG,IUFDCG)
      KFILE=KUFDEV(1:IUFDEV)//KUFDCG(1:IUFDCG)//'.dat'
      INQUIRE(FILE=KFILE,EXIST=LEX)
      IF(LEX) THEN
         KSTAT='OLD'
      ELSE
         KSTAT='NEW'
      ENDIF
      OPEN(21,FILE=KFILE,IOSTAT=IST,STATUS=KSTAT,
     &     ERR=10,FORM='FORMATTED')
      WRITE(6,600) '# NEW FILE (',KFILE,') IS CREATED'
      GOTO 20
 10   WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
      CLOSE(21)
      RETURN
C
 20   WRITE(21,610) 'DEV   = ',KUFDEV
      WRITE(21,610) 'SHOT  = ',KUFDCG
      WRITE(21,620) 'AVERW =' ,WTO1
      WRITE(21,620) 'RW    =' ,WTO2
CCC      WRITE(21,620) 'AVERW =' ,WTO2E
CCC      WRITE(21,620) 'RW    =' ,WTO2I
      WRITE(21,620) 'STDE  =' ,STDE
      WRITE(21,620) 'STDI  =' ,STDI
      WRITE(21,620) 'OFFE  =' ,OFFE
      WRITE(21,620) 'OFFI  =' ,OFFI
      CLOSE(21)
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
 600  FORMAT(' ',A12,A20,A12)
 610  FORMAT(' ',A8,A20)
 620  FORMAT(' ',A7,E15.7)
      RETURN
      END
