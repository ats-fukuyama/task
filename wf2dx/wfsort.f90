!     $Id: wfsort.f90,v 1.2 2011/07/19 17:54:28 maruyama Exp $

!     ***** SLAVE ROUTINE FOR WFSORT: FOR SUBSTITUTION *****

SUBROUTINE WFSRTS(I,J)
  
  use wfcomm
  implicit none
  integer :: NDELMX(3),I,J,IVELMX,IDELMX,KAELMX
  SAVE IVELMX,IDELMX,KAELMX,NDELMX
  
  IF(I.EQ.0) THEN
     IVELM(J)=IVELMX
     IDELM(J)=IDELMX
     KAELM(J)=KAELMX
     NDELM(1,J)=NDELMX(1)
     NDELM(2,J)=NDELMX(2)
     NDELM(3,J)=NDELMX(3)
  ELSEIF(J.EQ.0) THEN
     IVELMX=IVELM(I)
     IDELMX=IDELM(I)
     KAELMX=KAELM(I)
     NDELMX(1)=NDELM(1,I)
     NDELMX(2)=NDELM(2,I)
     NDELMX(3)=NDELM(3,I)
  ELSE
     IVELM(J)=IVELM(I)
     IDELM(J)=IDELM(I)
     KAELM(J)=KAELM(I)
     NDELM(1,J)=NDELM(1,I)
     NDELM(2,J)=NDELM(2,I)
     NDELM(3,J)=NDELM(3,I)
  ENDIF
  RETURN
END SUBROUTINE WFSRTS

!     ***** SLAVE ROUTINE FOR WFSORT: FOR EXCANGE *****

SUBROUTINE WFSRTX(I,J)
        
  use wfcomm
  implicit none
  integer :: NDELMX(3),IVELMX,I,IDELMX,KAELMX,J
      
  IVELMX=IVELM(I)
  IDELMX=IDELM(I)
  KAELMX=KAELM(I)
  NDELMX(1)=NDELM(1,I)
  NDELMX(2)=NDELM(2,I)
  NDELMX(3)=NDELM(3,I)

  IVELM(I)=IVELM(J)
  IDELM(I)=IDELM(J)
  KAELM(I)=KAELM(J)
  NDELM(1,I)=NDELM(1,J)
  NDELM(2,I)=NDELM(2,J)
  NDELM(3,I)=NDELM(3,J)

  IVELM(J)=IVELMX
  IDELM(J)=IDELMX
  KAELM(J)=KAELMX
  NDELM(1,J)=NDELMX(1)
  NDELM(2,J)=NDELMX(2)
  NDELM(3,J)=NDELMX(3)
  
  RETURN
END SUBROUTINE WFSRTX

!     ***** QUICK SORT BASED ON NUM RECIPE *****

SUBROUTINE WFSORT(N,ARR,SUB,SUBX)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer,parameter :: M=7
  integer,parameter :: NSTACK=50
  integer :: N,I,IR,J,JSTACK,K,L,ISTACK(NSTACK)
  real(rkind) :: ARR(N),A,TEMP
  EXTERNAL SUB,SUBX

  JSTACK=0
  L=1
  IR=N
1 IF(IR-L.LT.M)THEN
     DO J=L+1,IR
        A=ARR(J)
        CALL SUB(J,0)
        DO I=J-1,1,-1
           IF(ARR(I).LE.A)GOTO 2
           ARR(I+1)=ARR(I)
           CALL SUB(I,I+1)
        ENDDO
        I=0
2       ARR(I+1)=A
        CALL SUB(0,I+1)
     ENDDO
     if(JSTACK.eq.0) return
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
3    CONTINUE
     I=I+1
     IF(ARR(I).LT.A)GOTO 3
4    CONTINUE
     J=J-1
     IF(ARR(J).GT.A)GOTO 4
     IF(J.LT.I)GOTO 5
     TEMP=ARR(I)
     ARR(I)=ARR(J)
     ARR(J)=TEMP
     CALL SUBX(I,J)
     GOTO 3
5    ARR(L)=ARR(J)
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
END SUBROUTINE WFSORT

!     ***** LOCATE INDEX *****

SUBROUTINE WFLCAT(XX,N,X,J)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer :: J,N,JL,JM,JU
  real(rkind) :: X,XX(N)
  
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
END SUBROUTINE WFLCAT

!     ***** LOCATE INDEX MAXIMUM LOWER THAN VALUE *****

SUBROUTINE WFLCAT_MAX(XX,N,X,J)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer :: J,N,JL,JM,JU
  real(rkind) :: X,XX(N)

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
END SUBROUTINE WFLCAT_MAX

!     ***** LOCATE INDEX MINIMUM HIGHER THAN VALUE *****

SUBROUTINE WFLCAT_MIN(XX,N,X,J)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer :: J,N,JL,JM,JU
  real(rkind) :: X,XX(N)
  
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
END SUBROUTINE WFLCAT_MIN
