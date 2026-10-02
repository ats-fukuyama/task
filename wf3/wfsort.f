C     $Id$
C
C     ***** SLAVE ROUTINE FOR WFSORT: FOR SUBSTITUTION *****
C
      SUBROUTINE WFSRTS(I,J)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION NDELMX(4)
      SAVE IVELMX,IDELMX,KAELMX,NDELMX
C
      IF(I.EQ.0) THEN
         IVELM(J)=IVELMX
         IDELM(J)=IDELMX
         KAELM(J)=KAELMX
         NDELM(1,J)=NDELMX(1)
         NDELM(2,J)=NDELMX(2)
         NDELM(3,J)=NDELMX(3)
         NDELM(4,J)=NDELMX(4)
      ELSEIF(J.EQ.0) THEN
         IVELMX=IVELM(I)
         IDELMX=IDELM(I)
         KAELMX=KAELM(I)
         NDELMX(1)=NDELM(1,I)
         NDELMX(2)=NDELM(2,I)
         NDELMX(3)=NDELM(3,I)
         NDELMX(4)=NDELM(4,I)
      ELSE
         IVELM(J)=IVELM(I)
         IDELM(J)=IDELM(I)
         KAELM(J)=KAELM(I)
         NDELM(1,J)=NDELM(1,I)
         NDELM(2,J)=NDELM(2,I)
         NDELM(3,J)=NDELM(3,I)
         NDELM(4,J)=NDELM(4,I)
      ENDIF
      RETURN
      END
C
C     ***** SLAVE ROUTINE FOR WFSORT: FOR EXCANGE *****
C
      SUBROUTINE WFSRTX(I,J)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION NDELMX(4)
C
      IVELMX=IVELM(I)
      IDELMX=IDELM(I)
      KAELMX=KAELM(I)
      NDELMX(1)=NDELM(1,I)
      NDELMX(2)=NDELM(2,I)
      NDELMX(3)=NDELM(3,I)
      NDELMX(4)=NDELM(4,I)
C
      IVELM(I)=IVELM(J)
      IDELM(I)=IDELM(J)
      KAELM(I)=KAELM(J)
      NDELM(1,I)=NDELM(1,J)
      NDELM(2,I)=NDELM(2,J)
      NDELM(3,I)=NDELM(3,J)
      NDELM(4,I)=NDELM(4,J)
C
      IVELM(J)=IVELMX
      IDELM(J)=IDELMX
      KAELM(J)=KAELMX
      NDELM(1,J)=NDELMX(1)
      NDELM(2,J)=NDELMX(2)
      NDELM(3,J)=NDELMX(3)
      NDELM(4,J)=NDELMX(4)
C
      RETURN
      END
C
C     ***** QUICK SORT BASED ON NUM RECIPE *****
C
      SUBROUTINE WFSORT(N,ARR,SUB,SUBX)
C
      INTEGER N,M,NSTACK
      REAL*8 ARR(N)
      EXTERNAL SUB,SUBX
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,IR,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A,TEMP
C
      JSTACK=0
      L=1
      IR=N
1     IF(IR-L.LT.M)THEN
        DO J=L+1,IR
          A=ARR(J)
          CALL SUB(J,0)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A)GOTO 2
            ARR(I+1)=ARR(I)
            CALL SUB(I,I+1)
          ENDDO
          I=0
2         ARR(I+1)=A
          CALL SUB(0,I+1)
        ENDDO
        IF(JSTACK.EQ.0)RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=TEMP
        CALL SUBX(K,L+1)
        IF(ARR(L+1).GT.ARR(IR))THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=TEMP
          CALL SUBX(L+1,IR)
        ENDIF
        IF(ARR(L).GT.ARR(IR))THEN
          TEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=TEMP
          CALL SUBX(L,IR)
        ENDIF
        IF(ARR(L+1).GT.ARR(L))THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=TEMP
          CALl SUBX(L+1,L)
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        CALL SUB(L,0)
3       CONTINUE
          I=I+1
        IF(ARR(I).LT.A)GOTO 3
4       CONTINUE
          J=J-1
        IF(ARR(J).GT.A)GOTO 4
        IF(J.LT.I)GOTO 5
        TEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=TEMP
        CALL SUBX(I,J)
        GOTO 3
5       ARR(L)=ARR(J)
        ARR(J)=A
        CALL SUB(J,L)
        CALL SUB(0,J)
        JSTACK=JSTACK+2
        IF(JSTACK.GT.NSTACK) THEN
           WRITE(6,*) 'XX WFSORT NSTACK TOO SMALL IN SORT2'
        ENDIF
        IF(IR-I+1.GE.J-L)THEN
          ISTACK(JSTACK)=IR
          ISTACK(JSTACK-1)=I
          IR=J-1
        ELSE
          ISTACK(JSTACK)=J-1
          ISTACK(JSTACK-1)=L
          L=I
        ENDIF
      ENDIF
      GOTO 1
      END
C
C     ***** LOCATE INDEX *****
C
      SUBROUTINE WFLCAT(XX,N,X,J)
C
      INTEGER J,N
      REAL*8 X,XX(N)
      INTEGER JL,JM,JU
C
      JL=0
      JU=N+1
   10 IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      J=JL
      RETURN
      END
C
C     ***** LOCATE INDEX MAXIMUM LOWER THAN VALUE *****
C
      SUBROUTINE WFLCAT_MAX(XX,N,X,J)
C
      INTEGER J,N
      REAL*8 X,XX(N)
      INTEGER JL,JM,JU
C
      JL=0
      JU=N+1
   10 IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      J=JL
      IF(J.EQ.0) THEN
         J=1
         RETURN
      ENDIF
   20 IF(J.LT.N) THEN
         IF(XX(J+1).EQ.XX(J)) THEN
            J=J+1
            GOTO 20
         ENDIF
      ENDIF
      RETURN
      END
C
C     ***** LOCATE INDEX MINIMUM HIGHER THAN VALUE *****
C
      SUBROUTINE WFLCAT_MIN(XX,N,X,J)
C
      INTEGER J,N
      REAL*8 X,XX(N)
      INTEGER JL,JM,JU
C
      JL=0
      JU=N+1
   10 IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.LT.XX(JM)))THEN
          JU=JM
        ELSE
          JL=JM
        ENDIF
      GOTO 10
      ENDIF
      J=JU
      IF(J.EQ.N+1) THEN
         J=N
         RETURN
      ENDIF
   20 IF(J.GT.1) THEN
         IF(XX(J-1).EQ.XX(J)) THEN
            J=J-1
            GOTO 20
         ENDIF
      ENDIF
      RETURN
      END
