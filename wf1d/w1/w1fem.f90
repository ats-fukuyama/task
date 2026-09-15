MODULE w1fem

CONTAINS

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDA(IERR)
    USE w1comm
    USE libbnd
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    REAL(rkind):: DT0(2,2,2),DT1(2,2,2),DU0(2,2,2)
    INTEGER:: I,J,K,N,KML,L,N1,N2,M,NX
    REAL(rkind):: RKV,DX,DS01,DS02
    REAL(rkind):: DS11,DS12,DS11A,DS12A,DS21,DS22
    REAL(rkind):: DT11,DT12,DT11A,DT12A,DT21,DT22

    MWID=11
    MLEN=3*NXPMAX+4
    ALLOCATE(CF(MWID,MLEN))

    RKV=2.D6*PI*RF/VC

    DO I=1,2
       DO J=1,2
          DO K=1,2
             IF(I.EQ.J.AND.I.EQ.K) THEN
                DS0(I,J,K)=0.25D0
             ELSE
                DS0(I,J,K)=1.D0/12.D0
             ENDIF
             IF(J.EQ.K) THEN
                DS1(I,J,K)=1.D0/3.D0
             ELSE
                DS1(I,J,K)=1.D0/6.D0
             ENDIF
             IF(I.EQ.1) THEN
                DS1(I,J,K)=-DS1(I,J,K)
             ENDIF
             IF(I.EQ.J) THEN
                DS2(I,J,K)= 0.5D0
             ELSE
                DS2(I,J,K)=-0.5D0
             ENDIF
             IF(I.EQ.1) THEN
                IF(J.EQ.K) THEN
                   DT0(I,J,K)=1.D0/3.D0
                ELSE
                   DT0(I,J,K)=1.D0/6.D0
                ENDIF
             ELSE
                DT0(I,J,K)=0.D0
             ENDIF
             IF(I.EQ.1) THEN
                IF(J.EQ.1) THEN
                   DT1(I,J,K)=-0.5D0
                ELSE
                   DT1(I,J,K)= 0.5D0
                ENDIF
             ELSE
                DT1(I,J,K)=0.D0
             ENDIF
             IF(I.EQ.1.AND.J.EQ.1) THEN
                DU0(I,J,K)=0.5D0
             ELSE
                DU0(I,J,K)=0.D0
             ENDIF
          END DO
       END DO
    END DO

    DO I=1,MWID
       DO N=1,MLEN
          CF(I,N)=(0.D0,0.D0)
       END DO
    END DO
    DO N=1,MLEN
       CA(N)=(0.D0,0.0D0)
    END DO

    KML=0
    CF(KML+6,1)=CGIN(1,1)
    CF(KML+7,1)=CGIN(2,1)
    CF(KML+9,1)=(-1.D0,0.D0)
    CF(KML+5,2)=CGIN(1,3)
    CF(KML+6,2)=CGIN(2,3)
    CF(KML+9,2)=(-1.D0,0.D0)
    CF(KML+3,4)=CGIN(1,2)
    CF(KML+4,4)=CGIN(2,2)
    CF(KML+2,5)=CGIN(1,4)
    CF(KML+3,5)=CGIN(2,4)

    CF(KML+6,MLEN-1)=CGOT(1,1)
    CF(KML+7,MLEN-1)=CGOT(2,1)
    CF(KML+4,MLEN-1)=(-1.D0,0.D0)
    CF(KML+5,MLEN  )=CGOT(1,3)
    CF(KML+6,MLEN  )=CGOT(2,3)
    CF(KML+4,MLEN  )=(-1.D0,0.D0)
    CF(KML+8,MLEN-3)=CGOT(1,2)
    CF(KML+9,MLEN-3)=CGOT(2,2)
    CF(KML+7,MLEN-2)=CGOT(1,4)
    CF(KML+8,MLEN-2)=CGOT(2,4)
    CF(KML+6,MLEN-4)=1.D0

    CA(1    )=-CGIN(3,1)
    CA(2    )=-CGIN(3,3)
    CA(4    )=-CGIN(3,2)
    CA(5    )=-CGIN(3,4)
    CA(MLEN-1)=-CGOT(3,1)
    CA(MLEN  )=-CGOT(3,3)
    CA(MLEN-3)=-CGOT(3,2)
    CA(MLEN-2)=-CGOT(3,4)

    DO I=1,2
       DO J=1,2
          L=3*(J-I+2)-1
          DS01=DS0(I,J,1)
          DS02=DS0(I,J,2)
          DS11=DS1(I,J,1)+0.5D0*DS1(1,I,J)
          DS12=DS1(I,J,2)+0.5D0*DS1(2,I,J)
          DS11A=DS1(J,I,1)+0.5D0*DS1(1,J,I)
          DS12A=DS1(J,I,2)+0.5D0*DS1(2,J,I)
          DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
          DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
          DT11=DT1(I,J,1)+0.5D0*DT1(I,1,J)
          DT12=DT1(I,J,2)+0.5D0*DT1(I,2,J)
          DT11A=DT1(J,I,1)+0.5D0*DT1(J,1,I)
          DT12A=DT1(J,I,2)+0.5D0*DT1(J,2,I)

          DO NX=1,NXPMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             M=3*(NX+I-1)-1
             CF(L+1,M+1)=CF(L+1,M+1) &
                        +CD0(1,N1)*DU0(I,J,1)*DX &
                        +CD0(1,N2)*DU0(I,J,2)*DX
             CF(L+2,M+1)=CF(L+2,M+1) &
                        +CD0(2,N1)*DT0(I,J,1)*DX &
                        +CD0(2,N2)*DT0(I,J,2)*DX
             CF(L,  M+2)=CF(L,  M+2) &
                       -CD0(2,N1)*DT0(J,I,1)*DX &
                       -CD0(2,N2)*DT0(J,I,2)*DX
             CF(L+1,M+2)=CF(L+1,M+2) &
                        +CD0(3,N1)*DS01*DX &
                        +CD0(3,N2)*DS02*DX &
                        +CD2(3,N1)*DS21/DX &
                        +CD2(3,N2)*DS22/DX
             CF(L+1,M+3)=CF(L+1,M+3) &
                        +CD0(4,N1)*DS01*DX &
                        +CD0(4,N2)*DS02*DX &
                        +CD2(4,N1)*DS21/DX &
                        +CD2(4,N2)*DS22/DX
             CF(L+3,M+1)=CF(L+3,M+1) &
                        -CI*CD1(1,N1)*DT11 &
                        -CI*CD1(1,N2)*DT12
             CF(L+2,M+2)=CF(L+2,M+2) &
                        +CI*CD1(2,N1)*DS11 &
                        +CI*CD1(2,N2)*DS12
             CF(L-1,M+3)=CF(L-1,M+3) &
                        +CI*CD1(1,N1)*DT11A &
                        +CI*CD1(1,N2)*DT12A
             CF(L,  M+3)=CF(L,  M+3) &
                        +CI*CD1(2,N1)*DS11A &
                        +CI*CD1(2,N2)*DS12A
          END DO
       END DO
    END DO

!     WRITE(6,999) ((L,M,CF(L,M),L=30,41),M=195,197)
! 999 FORMAT((3(3X,I3,I4,1P2E14.6)))

    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

601 FORMAT('!! ERROR IN BANDCD : IND = ',I5)
  END SUBROUTINE W1BNDA

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWA(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    REAL(rkind):: DT0(2,2,2),DT1(2,2,2),DU0(2,2,2)
    INTEGER:: NS,NX,I,J,K,N1,N2,M,L
    REAL(rkind):: RKV,RCE,DS01,DS02,DX,PABSL,PFLXL
    REAL(rkind):: DS11,DS12,DS11A,DS12A,DS21,DS22
    REAL(rkind):: DT11,DT12,DT11A,DT12A
    COMPLEX(rkind):: CABSL,CFLXL
    COMPLEX(rkind):: CM11,CM12,CM21,CM22,CM33,CM13,CM23,CM31,CM32
    COMPLEX(rkind):: CD11,CD12,CD21,CD22,CD33,CD13,CD23,CD31,CD32

    RKV=2.D6*PI*RF/VC
    RCE=VC*EPS0

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABS(NX,NS)=0.D0
       END DO
    END DO
    DO NX=1,NXPMAX
       FLUX(NX)=0.D0
    END DO

    DO NX=1,NXPMAX
       CE2DA(NZ,NX,1)=CA(3*NX)
       CE2DA(NZ,NX,2)=CA(3*NX+1)
       CE2DA(NZ,NX,3)=CA(3*NX+2)
    END DO

    DO I=1,2
       DO J=1,2
          DO K=1,2
             IF(I.EQ.J.AND.I.EQ.K) THEN
                DS0(I,J,K)=0.25D0
             ELSE
                DS0(I,J,K)=1.D0/12.D0
             ENDIF
             IF(J.EQ.K) THEN
                DS1(I,J,K)=1.D0/3.D0
             ELSE
                DS1(I,J,K)=1.D0/6.D0
             ENDIF
             IF(I.EQ.1) THEN
                DS1(I,J,K)=-DS1(I,J,K)
             ENDIF
             IF(I.EQ.J) THEN
                DS2(I,J,K)= 0.5D0
             ELSE
                DS2(I,J,K)=-0.5D0
             ENDIF
             IF(I.EQ.1) THEN
                IF(J.EQ.K) THEN
                   DT0(I,J,K)=1.D0/3.D0
                ELSE
                   DT0(I,J,K)=1.D0/6.D0
                ENDIF
             ELSE
                DT0(I,J,K)=0.D0
             ENDIF
             IF(I.EQ.1) THEN
                IF(J.EQ.1) THEN
                   DT1(I,J,K)=-0.5D0
                ELSE
                   DT1(I,J,K)= 0.5D0
                ENDIF
             ELSE
                DT1(I,J,K)=0.D0
             ENDIF
             IF(I.EQ.1.AND.J.EQ.1) THEN
                DU0(I,J,K)=0.5D0
             ELSE
                DU0(I,J,K)=0.D0
             ENDIF
          END DO
       END DO
    END DO

    DO I=1,2
       DO J=1,2
          DS01=DS0(I,J,1)
          DS02=DS0(I,J,2)
          DS11=DS1(I,J,1)+0.5D0*DS1(1,I,J)
          DS12=DS1(I,J,2)+0.5D0*DS1(2,I,J)
          DS11A=DS1(J,I,1)+0.5D0*DS1(1,J,I)
          DS12A=DS1(J,I,2)+0.5D0*DS1(2,J,I)
          DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
          DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
          DT11=DT1(I,J,1)+0.5D0*DT1(I,1,J)
          DT12=DT1(I,J,2)+0.5D0*DT1(I,2,J)
          DT11A=DT1(J,I,1)+0.5D0*DT1(J,1,I)
          DT12A=DT1(J,I,2)+0.5D0*DT1(J,2,I)

          DO NS=1,NSMAX
             DO NX=1,NXPMAX-1
                N1=NX
                N2=NX+1
                DX=RKV*(XA(N2)-XA(N1))
                M=3*(NX+I-1)-1
                L=3*(NX+J-1)-1
                CM11       =+CM0(1,N1,NS)*DU0(I,J,1)*DX &
                            +CM0(1,N2,NS)*DU0(I,J,2)*DX
                CM12       =+CM0(2,N1,NS)*DT0(I,J,1)*DX &
                            +CM0(2,N2,NS)*DT0(I,J,2)*DX
                CM21       =-CM0(2,N1,NS)*DT0(J,I,1)*DX &
                            -CM0(2,N2,NS)*DT0(J,I,2)*DX
                CM22       =+CM0(3,N1,NS)*DS01*DX &
                            +CM0(3,N2,NS)*DS02*DX &
                            +CM2(3,N1,NS)*DS21/DX &
                            +CM2(3,N2,NS)*DS22/DX
                CM33       =+CM0(4,N1,NS)*DS01*DX &
                            +CM0(4,N2,NS)*DS02*DX &
                            +CM2(4,N1,NS)*DS21/DX &
                            +CM2(4,N2,NS)*DS22/DX
                CM13       =-CI*CM1(1,N1,NS)*DT11 &
                            -CI*CM1(1,N2,NS)*DT12
                CM23       =+CI*CM1(2,N1,NS)*DS11 &
                            +CI*CM1(2,N2,NS)*DS12
                CM31       =+CI*CM1(1,N1,NS)*DT11A &
                            +CI*CM1(1,N2,NS)*DT12A
                CM32       =+CI*CM1(2,N1,NS)*DS11A &
                            +CI*CM1(2,N2,NS)*DS12A

                CABSL &
                  =DCONJG(CA(M+1))*(CM11*CA(L+1)+CM12*CA(L+2)+CM13*CA(L+3)) &
                  +DCONJG(CA(M+2))*(CM21*CA(L+1)+CM22*CA(L+2)+CM23*CA(L+3)) &
                  +DCONJG(CA(M+3))*(CM31*CA(L+1)+CM32*CA(L+2)+CM33*CA(L+3))
                PABSL=-CI*RCE*CABSL
                IF(I.EQ.1.AND.J.EQ.1) THEN
                   PABS(NX  ,NS)=PABS(NX  ,NS)+PABSL
                ELSEIF(I.EQ.2.AND.J.EQ.2) THEN
                   PABS(NX+1,NS)=PABS(NX+1,NS)+PABSL
                ELSE
                   PABS(NX  ,NS)=PABS(NX  ,NS)+0.5D0*PABSL
                   PABS(NX+1,NS)=PABS(NX+1,NS)+0.5D0*PABSL
                ENDIF
             END DO
          END DO
       END DO
    END DO

    DO I=1,2
       DO J=1,2
          DS21=DS1(J,I,1)+0.25D0*DS1(1,I,J)
          DS22=DS1(J,I,2)+0.25D0*DS1(2,I,J)
          DS11=DS0(J,I,1)
          DS12=DS0(J,I,2)

          DO NX=1,NXPMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             L=3*(NX+I-2)+2
             M=3*(NX+J-2)+2
             CD11=+CD2(1,N1)*DS21+CD2(1,N2)*DS22
             CD12=+CD2(2,N1)*DS21+CD2(2,N2)*DS22
             CD21=-CD2(2,N1)*DS21-CD2(2,N2)*DS22
             CD22=+CD2(3,N1)*DS21+CD2(3,N2)*DS22
             CD33=+CD2(4,N1)*DS21+CD2(4,N2)*DS22
             CD13=0.D0
             CD23=+CI*CD1(2,N1)*DS11*DX+CI*CD1(2,N2)*DS12*DX
             CD31=+CI*CD1(1,N1)*DS11*DX+CI*CD1(1,N2)*DS12*DX
             CD32=0.D0

             CFLXL=DCONJG(CA(L+1))*(CD11*CA(M+1)+CD12*CA(M+2)+CD13*CA(M+3)) &
                  +DCONJG(CA(L+2))*(CD21*CA(M+1)+CD22*CA(M+2)+CD23*CA(M+3)) &
                  +DCONJG(CA(L+3))*(CD31*CA(M+1)+CD32*CA(M+2)+CD33*CA(M+3))
             PFLXL=CI*RCE*CFLXL/DX
             FLUX(NX)     =FLUX(NX)     +PFLXL
          END DO
       END DO
    END DO
    FLUX(NXPMAX)=FLUX(NXPMAX-1)
    RETURN
  END SUBROUTINE W1EPWA

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDC(IERR)
    USE w1comm
    USE libbnd
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    INTEGER:: I,J,K,N,KML,N1,N2,M,L,NX
    REAL(rkind):: RKV,DS01,DS02,DS21,DS22,DS11,DS12,DT11,DT12,DX

    RKV=2.D6*PI*RF/VC
    MWID=11
    MLEN=3*NXPMAX+4
    ALLOCATE(CF(MWID,MLEN))

    DO I=1,2
       DO J=1,2
          DO K=1,2
             IF(I.EQ.J.AND.I.EQ.K) THEN
                DS0(I,J,K)=0.25D0
             ELSE
                DS0(I,J,K)=1.D0/12.D0
             ENDIF
             IF(J.EQ.K) THEN
                DS1(I,J,K)=1.D0/3.D0
             ELSE
                DS1(I,J,K)=1.D0/6.D0
             ENDIF
             IF(I.EQ.1) THEN
                DS1(I,J,K)=-DS1(I,J,K)
             ENDIF
             IF(I.EQ.J) THEN
                DS2(I,J,K)= 0.5D0
             ELSE
                DS2(I,J,K)=-0.5D0
             ENDIF
          END DO
       END DO
    END DO

    DO I=1,11
       DO N=1,MLEN
          CF(I,N)=(0.D0,0.D0)
       END DO
    END DO
    DO  N=1,MLEN
       CA(N)=(0.D0,0.0D0)
    END DO

    KML=0
    CF(KML+6,1)=CGIN(1,1)
    CF(KML+7,1)=CGIN(2,1)
    CF(KML+9,1)=(-1.D0,0.D0)
    CF(KML+5,2)=CGIN(1,3)
    CF(KML+6,2)=CGIN(2,3)
    CF(KML+9,2)=(-1.D0,0.D0)
    CF(KML+3,4)=CGIN(1,2)
    CF(KML+4,4)=CGIN(2,2)
    CF(KML+2,5)=CGIN(1,4)
    CF(KML+3,5)=CGIN(2,4)

    CF(KML+6,MLEN-1)=CGOT(1,1)
    CF(KML+7,MLEN-1)=CGOT(2,1)
    CF(KML+4,MLEN-1)=(-1.D0,0.D0)
    CF(KML+5,MLEN  )=CGOT(1,3)
    CF(KML+6,MLEN  )=CGOT(2,3)
    CF(KML+4,MLEN  )=(-1.D0,0.D0)
    CF(KML+8,MLEN-3)=CGOT(1,2)
    CF(KML+9,MLEN-3)=CGOT(2,2)
    CF(KML+7,MLEN-2)=CGOT(1,4)
    CF(KML+8,MLEN-2)=CGOT(2,4)

    CA(1    )=-CGIN(3,1)
    CA(2    )=-CGIN(3,3)
    CA(4    )=-CGIN(3,2)
    CA(5    )=-CGIN(3,4)
    CA(MLEN-1)=-CGOT(3,1)
    CA(MLEN  )=-CGOT(3,3)
    CA(MLEN-3)=-CGOT(3,2)
    CA(MLEN-2)=-CGOT(3,4)

    DO I=1,2
       DO J=1,2
          L=3*(J-I+2)-1
          DS01=DS0(I,J,1)
          DS02=DS0(I,J,2)
          DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
          DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
          DS11=(DS1(J,I,1)+0.5D0*DS1(1,I,J))
          DS12=(DS1(J,I,2)+0.5D0*DS1(2,I,J))
          DT11=(DS1(I,J,1)+0.5D0*DS1(1,I,J))
          DT12=(DS1(I,J,2)+0.5D0*DS1(2,I,J))

          DO NX=1,NXPMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             M=3*(NX+I-1)-1
             CF(L+1,M+1)=CF(L+1,M+1) &
                        +CD0(1,N1)*DS01*DX &
                        +CD0(1,N2)*DS02*DX &
                        +CD2(1,N1)*DS21/DX &
                        +CD2(1,N2)*DS22/DX
             CF(L+2,M+1)=CF(L+2,M+1) &
                        +CD0(2,N1)*DS01*DX &
                        +CD0(2,N2)*DS02*DX &
                        +CD2(2,N1)*DS21/DX &
                        +CD2(2,N2)*DS22/DX
             CF(L,  M+2)=CF(L,  M+2) &
                        -CD0(2,N1)*DS01*DX &
                        -CD0(2,N2)*DS02*DX &
                        -CD2(2,N1)*DS21/DX &
                        -CD2(2,N2)*DS22/DX
             CF(L+1,M+2)=CF(L+1,M+2) &
                        +CD0(3,N1)*DS01*DX &
                        +CD0(3,N2)*DS02*DX &
                        +CD2(3,N1)*DS21/DX &
                        +CD2(3,N2)*DS22/DX
             CF(L+1,M+3)=CF(L+1,M+3) &
                        +CD0(4,N1)*DS01*DX &
                        +CD0(4,N2)*DS02*DX &
                        +CD2(4,N1)*DS21/DX &
                        +CD2(4,N2)*DS22/DX
             CF(L+3,M+1)=CF(L+3,M+1) &
                        -CI*CD1(1,N1)*DS11 &
                        -CI*CD1(1,N2)*DS12
             CF(L+2,M+2)=CF(L+2,M+2) &
                        +CI*CD1(2,N1)*DT11 &
                        +CI*CD1(2,N2)*DT12
             CF(L-1,M+3)=CF(L-1,M+3) &
                        +CI*CD1(1,N1)*DT11 &
                        +CI*CD1(1,N2)*DT12
             CF(L,  M+3)=CF(L,  M+3) &
                        +CI*CD1(2,N1)*DS11 &
                        +CI*CD1(2,N2)*DS12
          END DO
       END DO
    END DO

!     WRITE(6,999) ((L,M,CF(L,M),L=30,41),M=195,197)
! 999 FORMAT((3(3X,I3,I4,1P2E14.6)))

    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

  601 FORMAT('!! ERROR IN BANDCD : IND = ',I5)
  END SUBROUTINE W1BNDC

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWC(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    INTEGER:: NS,NX,I,J,K,N1,N2,M,L
    REAL(rkind):: RKV,RCE,DS01,DS02,DS21,DS22,DS11,DS12,DT11,DT12,DX
    REAL(rkind):: PABSL,PFLXL
    COMPLEX(rkind):: CABSL,CFLXL
    COMPLEX(rkind):: CM11,CM12,CM13,CM21,CM22,CM23,CM31,CM32,CM33
    COMPLEX(rkind):: CD11,CD12,CD13,CD21,CD22,CD23,CD31,CD32,CD33

    RKV=2.D6*PI*RF/VC
    RCE=VC*EPS0

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABS(NX,NS)=0.D0
       END DO
    END DO
    DO NX=1,NXPMAX
       FLUX(NX)=0.D0
    END DO

    DO NX=1,NXPMAX
       CE2DA(NZ,NX,1)=CA(3*NX)
       CE2DA(NZ,NX,2)=CA(3*NX+1)
       CE2DA(NZ,NX,3)=CA(3*NX+2)
    END DO

    DO I=1,2
       DO J=1,2
          DO K=1,2
             IF(I.EQ.J.AND.I.EQ.K) THEN
                DS0(I,J,K)=0.25D0
             ELSE
                DS0(I,J,K)=1.D0/12.D0
             ENDIF
             IF(J.EQ.K) THEN
                DS1(I,J,K)=1.D0/3.D0
             ELSE
                DS1(I,J,K)=1.D0/6.D0
             ENDIF
             IF(I.EQ.1) THEN
                DS1(I,J,K)=-DS1(I,J,K)
             ENDIF
             IF(I.EQ.J) THEN
                DS2(I,J,K)= 0.5D0
             ELSE
                DS2(I,J,K)=-0.5D0
             ENDIF
          END DO
       END DO
    END DO

    DO  I=1,2
       DO J=1,2
          DS01=DS0(I,J,1)
          DS02=DS0(I,J,2)
          DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
          DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
          DS11=(DS1(J,I,1)+0.5D0*DS1(1,I,J))
          DS12=(DS1(J,I,2)+0.5D0*DS1(2,I,J))
          DT11=(DS1(I,J,1)+0.5D0*DS1(1,I,J))
          DT12=(DS1(I,J,2)+0.5D0*DS1(2,I,J))

          DO NS=1,NSMAX
             DO NX=1,NXPMAX-1
                N1=NX
                N2=NX+1
                DX=RKV*(XA(N2)-XA(N1))
                M=3*(NX+I-1)-1
                L=3*(NX+J-1)-1
                CM11=+CM0(1,N1,NS)*DS01*DX &
                     +CM0(1,N2,NS)*DS02*DX &
                     +CM2(1,N1,NS)*DS21/DX &
                     +CM2(1,N2,NS)*DS22/DX
                CM12=+CM0(2,N1,NS)*DS01*DX &
                     +CM0(2,N2,NS)*DS02*DX &
                     +CM2(2,N1,NS)*DS21/DX &
                     +CM2(2,N2,NS)*DS22/DX
                CM21=-CM0(2,N1,NS)*DS01*DX &
                     -CM0(2,N2,NS)*DS02*DX &
                     -CM2(2,N1,NS)*DS21/DX &
                     -CM2(2,N2,NS)*DS22/DX
                CM22=+CM0(3,N1,NS)*DS01*DX &
                     +CM0(3,N2,NS)*DS02*DX &
                     +CM2(3,N1,NS)*DS21/DX &
                     +CM2(3,N2,NS)*DS22/DX
                CM33=+CM0(4,N1,NS)*DS01*DX &
                     +CM0(4,N2,NS)*DS02*DX &
                     +CM2(4,N1,NS)*DS21/DX &
                     +CM2(4,N2,NS)*DS22/DX
                CM13=-CI*CM1(1,N1,NS)*DS11 &
                     -CI*CM1(1,N2,NS)*DS12
                CM23=+CI*CM1(2,N1,NS)*DT11 &
                     +CI*CM1(2,N2,NS)*DT12
                CM31=+CI*CM1(1,N1,NS)*DT11 &
                     +CI*CM1(1,N2,NS)*DT12
                CM32=+CI*CM1(2,N1,NS)*DS11 &
                     +CI*CM1(2,N2,NS)*DS12

                CABSL &
                  =DCONJG(CA(M+1))*(CM11*CA(L+1)+CM12*CA(L+2)+CM13*CA(L+3)) &
                  +DCONJG(CA(M+2))*(CM21*CA(L+1)+CM22*CA(L+2)+CM23*CA(L+3)) &
                  +DCONJG(CA(M+3))*(CM31*CA(L+1)+CM32*CA(L+2)+CM33*CA(L+3))
                PABSL=-CI*RCE*CABSL
                IF(I.EQ.1.AND.J.EQ.1) THEN
                   PABS(NX  ,NS)=PABS(NX  ,NS)+PABSL
                ELSEIF(I.EQ.2.AND.J.EQ.2) THEN
                   PABS(NX+1,NS)=PABS(NX+1,NS)+PABSL
                ELSE
                   PABS(NX  ,NS)=PABS(NX  ,NS)+0.5D0*PABSL
                   PABS(NX+1,NS)=PABS(NX+1,NS)+0.5D0*PABSL
                ENDIF
             END DO
          END DO
       END DO
    END DO

    DO I=1,2
       DO J=1,2
          DS21=DS1(J,I,1)+0.25D0*DS1(1,I,J)
          DS22=DS1(J,I,2)+0.25D0*DS1(2,I,J)
          DS11=DS0(J,I,1)
          DS12=DS0(J,I,2)

          DO NX=1,NXPMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             L=3*(NX+I-2)+2
             M=3*(NX+J-2)+2
             CD11=+CD2(1,N1)*DS21+CD2(1,N2)*DS22
             CD12=+CD2(2,N1)*DS21+CD2(2,N2)*DS22
             CD21=-CD2(2,N1)*DS21-CD2(2,N2)*DS22
             CD22=+CD2(3,N1)*DS21+CD2(3,N2)*DS22
             CD33=+CD2(4,N1)*DS21+CD2(4,N2)*DS22
             CD13=0.D0
             CD23=+CI*CD1(2,N1)*DS11*DX+CI*CD1(2,N2)*DS12*DX
             CD31=+CI*CD1(1,N1)*DS11*DX+CI*CD1(1,N2)*DS12*DX
             CD32=0.D0

             CFLXL=DCONJG(CA(L+1))*(CD11*CA(M+1)+CD12*CA(M+2)+CD13*CA(M+3)) &
                  +DCONJG(CA(L+2))*(CD21*CA(M+1)+CD22*CA(M+2)+CD23*CA(M+3)) &
                  +DCONJG(CA(L+3))*(CD31*CA(M+1)+CD32*CA(M+2)+CD33*CA(M+3))
             PFLXL=CI*RCE*CFLXL/DX
             FLUX(NX)     =FLUX(NX)     +PFLXL
          END DO
       END DO
    END DO
    FLUX(NXPMAX)=FLUX(NXPMAX-1)
    RETURN
  END SUBROUTINE W1EPWC
END MODULE w1fem
