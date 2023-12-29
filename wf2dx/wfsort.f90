! wfsort.f90

MODULE wfsort

  PRIVATE
  PUBLIC wf_sort
  PUBLIC wf_sort_subst
  PUBLIC wf_sort_exchange
  PUBLIC wf_eval_index
  PUBLIC wf_eval_index_max
  PUBLIC wf_eval_index_min

CONTAINS

!     ***** SLAVE ROUTINE FOR WFSORT: FOR SUBSTITUTION *****

SUBROUTINE wf_sort_subst(I,J)
  
  use wfcomm
  implicit none
  integer :: node_nside_nelmX(3),I,J,IVELMX,IDELMX,KAELMX
  SAVE IVELMX,IDELMX,KAELMX,node_nside_nelmX
  
  IF(I.EQ.0) THEN
     IVELM(J)=IVELMX
     IDELM(J)=IDELMX
     KAELM(J)=KAELMX
     node_nside_nelm(1,J)=node_nside_nelmX(1)
     node_nside_nelm(2,J)=node_nside_nelmX(2)
     node_nside_nelm(3,J)=node_nside_nelmX(3)
  ELSEIF(J.EQ.0) THEN
     IVELMX=IVELM(I)
     IDELMX=IDELM(I)
     KAELMX=KAELM(I)
     node_nside_nelmX(1)=node_nside_nelm(1,I)
     node_nside_nelmX(2)=node_nside_nelm(2,I)
     node_nside_nelmX(3)=node_nside_nelm(3,I)
  ELSE
     IVELM(J)=IVELM(I)
     IDELM(J)=IDELM(I)
     KAELM(J)=KAELM(I)
     node_nside_nelm(1,J)=node_nside_nelm(1,I)
     node_nside_nelm(2,J)=node_nside_nelm(2,I)
     node_nside_nelm(3,J)=node_nside_nelm(3,I)
  ENDIF
  RETURN
END SUBROUTINE wf_sort_subst

!     ***** SLAVE ROUTINE FOR WFSORT: FOR EXCANGE *****

SUBROUTINE wf_sort_exchange(I,J)
        
  use wfcomm
  implicit none
  integer :: node_nside_nelmX(3),IVELMX,I,IDELMX,KAELMX,J
      
  IVELMX=IVELM(I)
  IDELMX=IDELM(I)
  KAELMX=KAELM(I)
  node_nside_nelmX(1)=node_nside_nelm(1,I)
  node_nside_nelmX(2)=node_nside_nelm(2,I)
  node_nside_nelmX(3)=node_nside_nelm(3,I)

  IVELM(I)=IVELM(J)
  IDELM(I)=IDELM(J)
  KAELM(I)=KAELM(J)
  node_nside_nelm(1,I)=node_nside_nelm(1,J)
  node_nside_nelm(2,I)=node_nside_nelm(2,J)
  node_nside_nelm(3,I)=node_nside_nelm(3,J)

  IVELM(J)=IVELMX
  IDELM(J)=IDELMX
  KAELM(J)=KAELMX
  node_nside_nelm(1,J)=node_nside_nelmX(1)
  node_nside_nelm(2,J)=node_nside_nelmX(2)
  node_nside_nelm(3,J)=node_nside_nelmX(3)
  
  RETURN
END SUBROUTINE wf_sort_exchange

!     ***** QUICK SORT BASED ON NUM RECIPE *****

SUBROUTINE wf_sort(N,ARR,SUB,SUBX)

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
END SUBROUTINE wf_sort

!     ***** evaluate INDEX *****

SUBROUTINE wf_eval_index(XX,N,X,J)

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
END SUBROUTINE wf_eval_index

!     ***** evaluate index MAXIMUM LOWER THAN VALUE *****

SUBROUTINE wf_eval_index_max(XX,N,X,J)

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
END SUBROUTINE wf_eval_index_max

!     ***** evaluate index MINIMUM HIGHER THAN VALUE *****

SUBROUTINE wf_eval_index_min(XX,N,X,J)

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
END SUBROUTINE wf_eval_index_min
END MODULE wfsort
