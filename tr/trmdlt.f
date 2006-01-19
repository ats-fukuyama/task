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
CCC      IF(MDLUF.EQ.2) THEN
         NTL=1
CCC      ENDIF
      MODE=0
      WSIM=0.D0
      WEXP=0.D0
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      RTE09 =0.5D0*(RTU(NTL,NRHO09,1)+RTU(NTL,NRHO09+1,1))
      RTI09 =0.5D0*(RTU(NTL,NRHO09,2)+RTU(NTL,NRHO09+1,2))
C      NRHO09=NRAMAX
C      RTE09 =PNSA(1)
C      RTI09 =PNSA(2)
      IF(MODE.EQ.0) THEN
         DO NR=1,NRMAX
            IF(RG(NR).GE.0.2D0.AND.RG(NR).LE.0.9D0) THEN
               RNESIM= 0.5D0*(RN(NR,1)+RN(NR+1,1))
               RNISIM= 0.5D0*(RN(NR,2)+RN(NR+1,2))
               RNEEXP= 0.5D0*(RN(NR,1)+RN(NR+1,1))
               RNIEXP= 0.5D0*(RN(NR,2)+RN(NR+1,2))
               TESIM = 0.5D0*(RT(NR,1)+RT(NR+1,1))-RTE09
               TISIM = 0.5D0*(RT(NR,2)+RT(NR+1,2))-RTI09
               TEEXP = 0.5D0*(RTU(NTL,NR,1)+RTU(NTL,NR+1,1))-RTE09
               TIEXP = 0.5D0*(RTU(NTL,NR,2)+RTU(NTL,NR+1,2))-RTI09
               WSIM=WSIM+1.5D0*( RNESIM*TESIM
     &                          +RNISIM*TISIM)*DVRHO(NR)*DR
               WEXP=WEXP+1.5D0*( RNEEXP*TEEXP
     &                          +RNIEXP*TIEXP)*DVRHO(NR)*DR
            ENDIF
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
            IF(RG(NR).LE.0.9D0) THEN
               RNESIM= 0.5D0*(RN(NR,1)+RN(NR+1,1))
               RNISIM= 0.5D0*(RN(NR,2)+RN(NR+1,2))
               RNEEXP= 0.5D0*(RN(NR,1)+RN(NR+1,1))
               RNIEXP= 0.5D0*(RN(NR,2)+RN(NR+1,2))
               TESIM = 0.5D0*(RT(NR,1)+RT(NR+1,1))-RTE09
               TISIM = 0.5D0*(RT(NR,2)+RT(NR+1,2))-RTI09
               TEEXP = 0.5D0*(RTU(NTL,NR,1)+RTU(NTL,NR+1,1))-RTE09
               TIEXP = 0.5D0*(RTU(NTL,NR,2)+RTU(NTL,NR+1,2))-RTI09
               WSIM=WSIM+1.5D0*( RNESIM*TESIM
     &                          +RNISIM*TISIM)*DVRHO(NR)*DR
               WEXP=WEXP+1.5D0*( RNEEXP*TEEXP
     &                          +RNIEXP*TIEXP)*DVRHO(NR)*DR
            ENDIF
         ENDDO
      ENDIF
      WTO1=WSIM/WEXP-1.D0
      WTO2=(WSIM/WEXP-1.D0)**2
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
         IF(RG(NR).GE.0.2D0.AND.RG(NR).LE.0.9D0) THEN
            TESIM=0.5D0*(RT(NR,1)+RT(NR+1,1))
            TISIM=0.5D0*(RT(NR,2)+RT(NR+1,2))
            TEEXP=0.5D0*(RTU(NTL,NR,1)+RTU(NTL,NR+1,1))
            TIEXP=0.5D0*(RTU(NTL,NR,2)+RTU(NTL,NR+1,2))
            RNUME=RNUME+TEEXP**2
            RNUMI=RNUMI+TIEXP**2
            RDESE=RDESE+(TESIM-TEEXP)**2
            RDESI=RDESI+(TISIM-TIEXP)**2
            RDEOE=RDEOE+(TESIM-TEEXP)
            RDEOI=RDEOI+(TISIM-TIEXP)
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
CCC      WRITE(21,620) 'AVERW =' ,WEXP*RKEV*1.D20*1.D-6
CCC      WRITE(21,620) 'RW    =' ,WSIM*RKEV*1.D20*1.D-6
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
