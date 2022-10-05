! w1exec7.f90

MODULE w1exec7

  USE w1comm,ONLY: rkind

  PRIVATE
  PUBLIC w1_exec7

  REAL(rkind):: DXD
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: CL
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NCLA
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XM,YK
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CS0,CS1,CS2,CS3
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: YX
  INTEGER:: NXDMIN,NXDMAX,NXWMAX,NXLMAX,NCLMAX

CONTAINS

  SUBROUTINE w1_exec7(NZP,IERR)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZP
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    CALL W1BCND
    CALL W1DSP
    CALL W1BND(IERR)
       IF(IERR.NE.0) RETURN
    CALL W1EPW(NZP)
    CALL W1CLCD(NZP)
    CALL W1CLP(NZP)
    RETURN
  END SUBROUTINE w1_exec7

!     ****** SET BOUNDARY CONDITIONS AT R=RA ******

  SUBROUTINE W1BCND
    USE w1comm
    IMPLICIT NONE
    INTEGER:: I,J
    REAL(rkind):: RW,RKV,RKPR,RNPR,RCE,RKKV2
    COMPLEX(rkind):: CKKV2,CKKV

    RW=2.D6*PI*RF
    RKV=RW/VC
    RKPR=RKZ
    RNPR=VC*RKPR/RW
    RCE=VC*EPS0
    RKKV2=RNPR*RNPR-1.D0
    IF(RKKV2.GE.0.D0) THEN
       CKKV2=RKKV2
       CKKV =SQRT(RKKV2)
    ELSE
       CKKV2=RKKV2
       CKKV=CI*SQRT(-RKKV2)
    ENDIF
!    IF(ABS(WALLR).GT.1.D-12) THEN
!       CKKV2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
!       CKKV =SQRT(CKKV2)
!    ELSE 
!       CKKV2=RNPR*RNPR-1.D0
!       CKKV =SQRT(CKKV2)
!    ENDIF

    DO I=1,4    ! I=1:DEY I=2:DEZ I=3:EY I=4:EZ
       DO J=1,4 ! J=1:N-1, J=2:N, J=3:N+1, J=4:RHS
          CGIN(J,I)=(0.D0,0.D0)
       END DO
    END DO
    DO I=1,4    ! I=1:EY I=2:EZ I=3:DEY I=4:DEZ
       DO J=1,4 ! J=1:N-1, J=2:N, J=3:N+1, J=4:RHS
          CGOT(J,I)=(0.D0,0.D0)
       END DO
    END DO

    IF(MDLWG.EQ.4) THEN    ! HFS Xmode
       CGIN(2,1)=-CI*CKKV             ! y component derivative
       CGIN(3,1)=+CI*CKKV
    END IF

    IF(MDLWG.EQ.3) THEN    ! HFS Omode
       CGIN(2,2)=-CI*CKKV             ! z component derivative
       CGIN(3,2)=+CI*CKKV
    END IF
 
    IF(MDLWG.EQ.2) THEN    ! LFS Xmode
       CGOT(1,3)=-CI*CKKV             ! y component derivative
       CGOT(2,3)=+CI*CKKV
    END IF

    IF(MDLWG.EQ.1) THEN      ! LFS Omode
       CGOT(1,4)=-CI*CKKV             ! z component derivative
       CGOT(2,4)=+CI*CKKV
    END IF

    RETURN
  END SUBROUTINE W1BCND

!     ******* LOCAL PLASMA PARAMETERS for WARM *******

  SUBROUTINE W1DSP
    USE w1comm
    USE w1alfa
    USE libdsp
    USE w1lib
    IMPLICIT NONE

!      CM0 , CD0         CM1 , CD1          CM2 , CD2
!      1   2   0         0   0   1          1   2   0
!    -(2)  3   0         0   0   2        -(2)  3   0
!      0   0   4         1 -(2)  0          0   0   4

    INTEGER:: NS,NX
    REAL(rkind):: RW,RKPR,RNPR
    COMPLEX(rkind):: CW,CWP0,CWC0,CWP,CWC
    
    RW  = 2.D6*PI*RF
    RKPR=RKZ
    RNPR=VC*RKPR/RW

    DO NS=1,NSMAX
       CW  = (1.D0 + CI*PZCL(NS))*RW
       CWP0 = 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0*CW*RW)
       CWC0 = AEE*PZ(NS)*BB/(AMP*PA(NS)*CW)
       DO NX=1,NXMAX
          CWP  = CWP0*PROFPN(NX,NS)
          CWC  = CWC0*PROFB(NX)
          CM0(1,NX,NS) = -CWP/(1.D0-CWC*CWC)
          CM0(2,NX,NS) = -CI*CWP*CWC/(1.D0-CW*CWC)
          CM0(3,NX,NS) = -CWP/(1.D0-CWC*CWC)
          CM0(4,NX,NS) = -CWP
          CM1(1,NX,NS) = 0.D0
          CM1(2,NX,NS) = 0.D0
          CM2(1,NX,NS) = 0.D0
          CM2(2,NX,NS) = 0.D0
          CM2(3,NX,NS) = 0.D0
          CM2(4,NX,NS) = 0.D0
       END DO
    END DO

    DO NX = 1 , NXMAX
       CD0( 1 , NX ) =  1.D0 - RNPR*RNPR
       CD0( 2 , NX ) =  0.D0
       CD0( 3 , NX ) =  1.D0 - RNPR*RNPR
       CD0( 4 , NX ) =  1.D0
       CD1( 1 , NX ) =  RNPR
       CD1( 2 , NX ) =  0.D0
       CD2( 1 , NX ) =  0.D0
       CD2( 2 , NX ) =  0.D0
       CD2( 3 , NX ) = -1.D0
       CD2( 4 , NX ) = -1.D0
    END DO

    DO NS = 1 , NSMAX
       DO NX = 1 , NXMAX
          CD0( 1 , NX ) = CD0( 1 , NX ) + CM0( 1 , NX , NS )
          CD0( 2 , NX ) = CD0( 2 , NX ) + CM0( 2 , NX , NS )
          CD0( 3 , NX ) = CD0( 3 , NX ) + CM0( 3 , NX , NS )
          CD0( 4 , NX ) = CD0( 4 , NX ) + CM0( 4 , NX , NS )
          CD1( 1 , NX ) = CD1( 1 , NX ) + CM1( 1 , NX , NS )
          CD1( 2 , NX ) = CD1( 2 , NX ) + CM1( 2 , NX , NS )
          CD2( 1 , NX ) = CD2( 1 , NX ) + CM2( 1 , NX , NS )
          CD2( 2 , NX ) = CD2( 2 , NX ) + CM2( 2 , NX , NS )
          CD2( 3 , NX ) = CD2( 3 , NX ) + CM2( 3 , NX , NS )
          CD2( 4 , NX ) = CD2( 4 , NX ) + CM2( 4 , NX , NS )
       END DO
    END DO
    RETURN
  END SUBROUTINE W1DSP

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BND(IERR)
    USE w1comm
    USE libbnd
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    REAL(rkind):: DT0(2,2,2),DT1(2,2,2),DU0(2,2,2)
    INTEGER:: I,J,K,N,L,N1,N2,M,NX
    REAL(rkind):: RW,RKV,DX,DS01,DS02
    REAL(rkind):: DS11,DS12,DS11A,DS12A,DS21,DS22
    REAL(rkind):: DT11,DT12,DT11A,DT12A,DT21,DT22

    MWID=11
    MLEN=3*NXMAX+4
    ALLOCATE(CF(MWID,MLEN))
    MCEN=6

    RW=2.D6*PI*RF
    RKV=RW/VC

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

    DO I=1,MLEN
       DO J=1,MWID
          CF(J,I)=(0.D0,0.D0)
       END DO
    END DO
    DO I=1,MLEN
       CA(I)=(0.D0,0.0D0)
    END DO

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

          DO NX=1,NXMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
!             WRITE(6,'(A,2I5,1P4E12.4)') 'N1,N2,XA,DX,RKV=', &
!                  N1,N2,XA(N1),XA(N2),DX,RKV
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

    IF(MDLWG.EQ.4) THEN
       DO I=1,MWID
          CF(I,1)=0.D0
       END DO
       CF(MCEN,1)=1.D0
       CA(4)=-CF(MCEN  ,4)*CFWG4
       CF(MCEN  ,4)=0.D0
    ELSEIF(MDLWG.EQ.-4) THEN
       CF(MCEN  ,1)=CGIN(2,1)
       CF(MCEN+3,1)=CGIN(3,1)
       CF(MCEN-3,4)=CF(MCEN  ,4)
       CA(4)=-CF(MCEN  ,4)*CFWG4
       CF(MCEN  ,4)=0.D0
    ELSE
       DO I=1,MWID
          CF(I,1)=0.D0
          CF(I,4)=0.D0
       END DO
       CF(MCEN,1)=1.D0
       CF(MCEN,4)=1.D0
    END IF

    IF(MDLWG.EQ.3) THEN
       DO I=1,MWID
          CF(I,2)=0.D0
       END DO
       CF(MCEN,2)=1.D0
       CA(5)=-CF(MCEN  ,5)*CFWG3
       CF(MCEN  ,5)=0.D0
    ELSEIF(MDLWG.EQ.-3) THEN
       CF(MCEN  ,2)=CGIN(2,2)
       CF(MCEN+3,2)=CGIN(3,2)
       CF(MCEN-3,5)=CF(MCEN  ,5)
       CA(5)=-CF(MCEN  ,5)*CFWG3
       CF(MCEN  ,5)=0.D0
    ELSE
       DO I=1,MWID
          CF(I,2)=0.D0
          CF(I,5)=0.D0
       END DO
       CF(MCEN,2)=1.D0
       CF(MCEN,5)=1.D0
    END IF

    CF(MCEN,  MLEN-4)=1.D0    ! ER

    IF(MDLWG.EQ.2) THEN
       CA(MLEN-3)=-CF(MCEN,MLEN-3)*CFWG2
       CF(MCEN,  MLEN-3)=0.D0
       DO I=1,MWID
          CF(I,MLEN-1)=0.D0
!          CF(I,MLEN-2)=0.D0
!          CF(I,MLEN  )=0.D0
!          CF(I,1)=0.D0
!          CF(I,2)=0.D0
!          CF(I,4)=0.D0
!          CF(I,5)=0.D0
       END DO
       CF(MCEN,MLEN-1)=1.D0
!       CF(MCEN,MLEN-2)=1.D0
!       CF(MCEN,MLEN  )=1.D0
!       CF(MCEN,1)=1.D0
!       CF(MCEN,2)=1.D0
!       CF(MCEN,4)=1.D0
!       CF(MCEN,5)=1.D0
    ELSEIF(MDLWG.EQ.-2) THEN
       CF(MCEN+2,MLEN-3)=CF(MCEN,  MLEN-3)
       CA(MLEN-3)=-CF(MCEN,MLEN-3)*CFWG2
       CF(MCEN,  MLEN-3)=0.D0
       CF(MCEN-2,MLEN-1)=CGOT(1,3)
       CF(MCEN,  MLEN-1)=CGOT(2,3)
    ELSE
       DO I=1,MWID
          CF(I,MLEN-3)=0.D0
          CF(I,MLEN-1)=0.D0
       END DO
       CF(MCEN,MLEN-3)=1.D0
       CF(MCEN,MLEN-1)=1.D0
    END IF

    IF(MDLWG.EQ.1) THEN
       CA(MLEN-2)=-CF(MCEN,  MLEN-2)*CFWG1
       CF(MCEN,  MLEN-2)=0.D0
       DO I=1,MWID
          CF(I,MLEN  )=0.D0
       END DO
       CF(MCEN,MLEN  )=1.D0
    ELSEIF(MDLWG.EQ.-1) THEN
       CF(MCEN-2,MLEN  )=CGOT(1,4)
       CF(MCEN,  MLEN  )=CGOT(2,4)
       CF(MCEN+2,MLEN-2)=CF(MCEN,  MLEN-2)
       CA(MLEN-2)=-CF(MCEN,  MLEN-2)*CFWG1
       CF(MCEN,  MLEN-2)=0.D0
    ELSE
       DO I=1,MWID
          CF(I,MLEN-2)=0.D0
          CF(I,MLEN  )=0.D0
       END DO
       CF(MCEN,MLEN-2)=1.D0
       CF(MCEN,MLEN  )=1.D0
    END IF

    CA(3*(NXANT1-1)+2+2)=CI*CFJY1/(RW*EPS0)
    CA(3*(NXANT1-1)+2+3)=CI*CFJZ1/(RW*EPS0)
    CA(3*(NXANT2-1)+2+2)=CI*CFJY2/(RW*EPS0)
    CA(3*(NXANT2-1)+2+3)=CI*CFJZ2/(RW*EPS0)

    IF(NZMAX.EQ.1) WRITE(6,'(A,2I5)') 'MLEN,MWID=',MLEN,MWID
    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

601 FORMAT('!! ERROR IN BANDCD : IND = ',I5)
  END SUBROUTINE W1BND

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPW(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    REAL(rkind):: DT0(2,2,2),DT1(2,2,2),DU0(2,2,2)
    INTEGER:: NS,NX,I,J,K,N1,N2,M,L
    REAL(rkind):: RW,RKV,RCE,DS01,DS02,DX,PABSL,PFLXL
    REAL(rkind):: DS11,DS12,DS11A,DS12A,DS21,DS22
    REAL(rkind):: DT11,DT12,DT11A,DT12A
    COMPLEX(rkind):: CABSL,CFLXL
    COMPLEX(rkind):: CM11,CM12,CM21,CM22,CM33,CM13,CM23,CM31,CM32
    COMPLEX(rkind):: CD11,CD12,CD21,CD22,CD33,CD13,CD23,CD31,CD32

    RW=2.D6*PI*RF
    RKV=RW/VC
    RCE=VC*EPS0

    DO NS=1,NSMAX
       DO NX=1,NXMAX
          PABS(NX,NS)=0.D0
       END DO
    END DO
    DO NX=1,NXMAX
       FLUX(NX)=0.D0
    END DO

    DO NX=1,NXMAX
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
             DO NX=1,NXMAX-1
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

          DO NX=1,NXMAX-1
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
    FLUX(NXMAX)=FLUX(NXMAX-1)
    RETURN
  END SUBROUTINE W1EPW

!     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******

  SUBROUTINE W1CLP(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS,I
    REAL(rkind):: FACT

    IF(NSYM.NE.0.AND.NZ.NE.1.AND.NZ.NE.(NZMAX/2+1)) THEN
       FACT = 2.0*RZ
    ELSEIF(NSYM.EQ.-1.AND.NZ.EQ.(NZMAX/2+1)) THEN
       FACT = 0.0
    ELSE
       FACT =     RZ
    ENDIF

    DO NS=1,NSMAX
       DO NX=1,NXMAX
          PABSX(NX,NS)=PABSX(NX,NS)+PABS(NX,NS)*FACT
       END DO
    END DO

    DO NX=1,NXMAX
       FLUXX(NX)=FLUXX(NX)+FLUX(NX)*FACT
    END DO

    DO NS=1,NSMAX
       PAK(NZ,NS)=0.D0
       DO NX=1,NXMAX
          PAK(NZ,NS)=PAK(NZ,NS)+PABS(NX,NS)*RZ
       END DO
    END DO

    CPANTK(NZ)=-RZ*DCONJG(CE2DA(NZ,NXANT1,2))*CFJY1 &
               -RZ*DCONJG(CE2DA(NZ,NXANT2,2))*CFJY2 &
               -RZ*DCONJG(CE2DA(NZ,NXANT1,3))*CFJZ1 &
               -RZ*DCONJG(CE2DA(NZ,NXANT2,3))*CFJZ2

    PAKT(NZ,3)=0.D0
    DO NS=1,NSMAX
       PAKT(NZ,3)=PAKT(NZ,3)+PAK(NZ,NS)
    END DO

    PAKT(NZ,2)=0.
    DO I=1,3
       PAKT(NZ,2)=PAKT(NZ,2)+( FLUXX(NXMAX) &
                              -FLUXX(1          ))
    END DO
    PAKT(NZ,1)=PAKT(NZ,2)+PAKT(NZ,3)

    RETURN
  END SUBROUTINE W1CLP

!     ****** CALCULATE DRIVEN CURRENT ******

  SUBROUTINE W1CLCD(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS
    REAL(rkind):: XR,VTE,VPH,W,EFCD,RKPR
    REAL(rkind):: Y0,Y1,Y2,Y3,E0,E1,E2,E3,RLNLMD,AJCD

    AJCDK(NZ)=0.D0

    DO NX=1,NXMAX
       XR=XAM(NX)/RR
       DO NS=1,NSMAX
          IF(IELEC(NS).EQ.1) THEN
             RKPR=RKZ
             IF(ABS(RKPR).LT.1.D-6) RKPR=1.D-6
             VTE=SQRT(PROFTR(NX,NS)*AEE*1.D3/AME)
             VPH=2.D0*PI*RF*1.D6/RKPR
             W=ABS(VPH/VTE)
             IF(ABS(VPH).GT.VC) THEN
                EFCD=0.D0
             ELSE
                IF(WVYSIZ.LE.0.D0) THEN
                   EFCD=W1CDEF(W,ZEFF,XR,0.D0,NCDTYP)
                ELSE
                   Y0=0.00D0*WVYSIZ/RR
                   Y1=0.25D0*WVYSIZ/RR
                   Y2=0.50D0*WVYSIZ/RR
                   Y3=0.75D0*WVYSIZ/RR
                   E0=W1CDEF(W,ZEFF,XR,Y0,NCDTYP)
                   E1=W1CDEF(W,ZEFF,XR,Y1,NCDTYP)
                   E2=W1CDEF(W,ZEFF,XR,Y2,NCDTYP)
                   E3=W1CDEF(W,ZEFF,XR,Y3,NCDTYP)
!     *** WEIGHTING WITH PARABOLIC ELECTRIC FIELD PROFILE ***
                   EFCD=(256.D0*E0+450.D0*E1+288.D0*E2+98.D0*E3)/1092.D0
!     *** WEIGHTING WITH PARABOLIC POWER DENSITY PROFILE ***
!                   EFCD=(16.D0*E0+30.D0*E1+24.D0*E2+14.D0*E3)/84.D0
                ENDIF
             ENDIF
             IF(PROFPN(NX,1).LE.0.D0) THEN
                AJCD=0.D0
             ELSE
                RLNLMD=16.1D0 - 1.15D0*LOG10(PROFPN(NX,1)) &
                              + 2.30D0*LOG10(PROFTR(NX,1))
                AJCD=0.384D0*PROFTR(NX,NS)*EFCD/(PROFPN(NX,1)*RLNLMD) &
                     *PABS(NX,NS)
             END IF
             IF(VPH.LT.0.D0) AJCD=-AJCD
             AJCDX(NX)=AJCDX(NX)+AJCD
             AJCDK(NZ)=AJCDK(NZ)+AJCD
          ENDIF
       END DO
    END DO
    RETURN
  END SUBROUTINE W1CLCD

!     ****** CURRENT DRIVE EFFICIENCY ******

!      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
!      Z  = ZEFF   : EFFECTIVE Z
!      XR = X / RR : NORMALIZED X
!      YR = Y / RR : NORMALIZED Y
!      ID : 0 : LANDAU DAMPING
!           1 : TTMP

  FUNCTION W1CDEF(WT,Z,XR,YR,ID)
    USE w1comm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: WT,Z,XR,YR
    INTEGER,INTENT(IN):: ID
    REAL(rkind):: W1CDEF,R,D,C,A,RM,RC,W,EFF0,EFF1,Y2,Y1,EFF2,YT,ARG,EFF3

    R=SQRT(XR*XR+YR*YR)
    IF(ID.EQ.0) THEN
       D=3.D0/Z
       C=3.83D0
       A=0.D0
       RM=1.38D0
       RC=0.389D0
    ELSE
       D=11.91D0/(0.678D0+Z)
       C=4.13D0
       A=12.3D0
       RM=2.48D0
       RC=0.0987D0
    ENDIF
    IF(WT.LE.1.D-20) THEN
       W=1.D-20
    ELSE
       W=WT
    ENDIF
    EFF0=D/W+C/Z**0.707D0+4.D0*W*W/(5.D0+Z)
    EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)

    Y2=(R+XR)/(1.D0+R)
    IF(Y2.LT.0.D0) Y2=0.D0
    Y1=SQRT(Y2)
    EFF2=1.D0+A*(Y1/W)**3

    IF(Y2.LE.1.D-20) THEN
       YT=(1.D0-Y2)*WT*WT/1.D-60
    ELSE
       YT=(1.D0-Y2)*WT*WT/Y2
    ENDIF
    IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
       ARG=(RC*YT)**RM
       IF(ARG.LE.100.D0) THEN
          EFF3=1.D0-MIN(EXP(-ARG),1.D0)
       ELSE
          EFF3=1.D0
       ENDIF
    ELSE
       EFF3=1.D0
    ENDIF

    W1CDEF=EFF0*EFF1*EFF2*EFF3
    RETURN
  END FUNCTION W1CDEF

END MODULE w1exec7
