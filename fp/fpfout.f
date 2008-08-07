C     $Id$
C
C ***********************************************************
C
C                FILE OUTPUT of GRAPH DATA
C
C ***********************************************************
C
C ***********************************************************
C
C                        RADIAL PROFILE
C
C ***********************************************************
C
      SUBROUTINE FPFOTR(STRING,FR)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FR(NRM,NTG1M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
      CALL FPFOUX(STRING,NRMAX,NTG1,FR,NRM)
C
      RETURN
      END
C
C ***********************************************************
C
C                        TIME EVOLUTION
C
C ***********************************************************
C
      SUBROUTINE FPFOTT(STRING,FT)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION  FT(NTG2M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
      CALL FPFOUX(STRING,NTG2,1,FT,NTG2M)
      RETURN
      END
C
C ***********************************************************
C
C                        Momentum Dependence
C
C ***********************************************************
C
      SUBROUTINE FPFOTP(STRING,FG)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(NTHM,NPM,NRM),FL(NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
    1 WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
C
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            FL(NP,NR)=FG(NTH,NP,NR)
         ENDDO
      ENDDO
C
      CALL FPFOUX(STRING,NPMAX,NRMAX,FL,NPM)
      GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPFOTC(STRING,FG,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*)
      DIMENSION FL(NPM,NTHM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
C
      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ENDIF         
C
      CALL FPFOUX(STRING,NPG,NTHG,FL,NPM)
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPFOTC2(STRING,FG,FH,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*),FH(*)
      DIMENSION FL(NPM,NTHM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
C
      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NMP1=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               NMP2=NPM*NTHM*(NR-1)+NTHM* NP   +NTH
               FGA=0.5D0*(FG(NMP1)+FG(NMP2))
               NMT1=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               NMT2=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH+1
               FHA=0.5D0*(FH(NMT1)+FH(NMT2))
               FL(NP,NTH)=SQRT(FGA**2+FHA**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=SQRT(FG(NM)**2+FH(NM)**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=SQRT(FG(NM)**2+FH(NM)**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=LOG10(ABS(SQRT(FG(NM)**2+FH(NM)**2)))
            ENDDO
         ENDDO
      ENDIF         
C
      CALL FPFOUX(STRING,NPG,NTHG,FL,NPM)
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        FILE OUTPUT EXEC
C
C ***********************************************************
C
      SUBROUTINE FPFOUX(STRING,N1MAX,N2MAX,FL,N1M)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FL(N1M,*)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
      WRITE(22,'(A)') STRING
      WRITE(22,'(2I10)') N1MAX,N2MAX
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1,N2),N1=1,N1MAX)
      ENDDO
      RETURN
      END
