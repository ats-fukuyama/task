! w1exec9.f90

MODULE w1exec9

  USE w1comm,ONLY: rkind
  USE w1fflr,ONLY: w1fnmn

  PRIVATE
  PUBLIC w1_exec9

  REAL(rkind):: DXD
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: CL
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NCLA
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XM,YK
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CS0,CS1,CS2,CS3
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: YX
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: SF,SG,AF,AG
  INTEGER:: NXDMIN,NXDMAX,NXWMAX,NXLMAX,NCLMAX

CONTAINS

  SUBROUTINE w1_exec9(NZ,IERR)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    CALL W1DSPC
    CALL W1BNDC(NZ,IERR)
       IF(IERR.NE.0) RETURN
    CALL W1EPWC(NZ)
    CALL W1CLCD(NZ)
    CALL W1CLPW(NZ)
    RETURN
  END SUBROUTINE w1_exec9

!     ******* LOCAL PLASMA PARAMETERS *******

  SUBROUTINE W1DSPC
    USE w1comm
    USE w1lib
    IMPLICIT NONE

!      CM0 , CD0         CM1 , CD1          CM2 , CD2
!      1   2   0         0   0   1          1   2   0
!    -(2)  3   0         0   0   2        -(2)  3   0
!      0   0   4         1 -(2)  0          0   0   4

    INTEGER:: NS,NX
    REAL(rkind)::RT2,RW,WPF,WCF,WP,WC
    COMPLEX(rkind):: CW

    RT2= SQRT(2.D0)
    RW=2.D6*PI*RF

    DO NS=1,NSMAX
       CW=RW*DCMPLX(1.D0,PZCL(NS))
       WPF= 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0)
       WCF= AEE*PZ(NS)*BB/(AMP*PA(NS))

       DO NX=1,NXPMAX
          WP  = WPF*PROFPN(NX,NS)
          WC  = WCF*PROFB(NX)
          CM0(1,NX,NS)=-(WP/RW)*CW/(CW*CW-WC*WC)
          CM0(2,NX,NS)=-CI*(WP/RW)*WC/(CW*CW-WC*WC)
          CM0(3,NX,NS)=-(WP/RW)*CW/(CW*CW-WC*WC)
          CM0(4,NX,NS)=-(WP/RW)/CW
       END DO
    END DO

    DO NX = 1 , NXPMAX
       CD0( 1 , NX ) =  0.D0
       CD0( 2 , NX ) =  0.D0
       CD0( 3 , NX ) =  0.D0
       CD0( 4 , NX ) =  0.D0
    END DO

    DO NS = 1 , NSMAX
       DO NX = 1 , NXPMAX
          CD0( 1 , NX ) = CD0( 1 , NX ) + CM0( 1 , NX , NS )
          CD0( 2 , NX ) = CD0( 2 , NX ) + CM0( 2 , NX , NS )
          CD0( 3 , NX ) = CD0( 3 , NX ) + CM0( 3 , NX , NS )
          CD0( 4 , NX ) = CD0( 4 , NX ) + CM0( 4 , NX , NS )
       END DO
    END DO
    RETURN
  END SUBROUTINE W1DSPC

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDC(NZ,IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2,2),DS1(2,2,2),DS2(2,2,2)
    REAL(rkind):: DT0(2,2,2),DT1(2,2,2),DU0(2,2,2)
    REAL(rkind):: DV0(2,2),DV1(2,2),DV2(2,2),DV3(2,2)

    INTEGER:: I,J,K,NX,L,N1,N2,M,NS
    REAL(rkind):: RW,RKV,DTT0,DSS0,DTS1,DTS2,DTTW,RKPR,RNPR,DX,RCE
    COMPLEX(rkind):: CKKV,CKKV2

    MWID=11
    MCEN=6
    MLEN=3*NXMAX+4
    ALLOCATE(CF(MWID,MLEN))
    ALLOCATE(CFS(MWID,MLEN,NSMAX))

    RW=2.D6*PI*RF
    RKV=RW/VC
    RKPR=RKZ
    IF(ABS(RKPR).LE.1.D-6) RKPR=1.D-6
    RNPR=VC*RKPR/RW
    RCE=VC*EPS0

    IF(ABS(WALLR).GT.1.D-12) THEN
       IF(1.D0-RNPR*RNPR.GT.0.D0) THEN
          CKKV=CI*SQRT(1.D0-RNPR*RNPR+DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR)))
       ELSE
          CKKV=-SQRT(RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR)))
       END IF
    ELSE
       IF(1.D0-RNPR*RNPR.GT.0.D0) THEN
          CKKV=CI*SQRT(1.D0-RNPR*RNPR)
       ELSE
          CKKV=-SQRT(RNPR*RNPR-1.D0)
       END IF
    END IF
    CKKV2=CKKV*CKKV

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

    DV0(1,1)=1.D0/3.D0
    DV0(2,1)=1.D0/6.D0
    DV0(1,2)=1.D0/6.D0
    DV0(2,2)=1.D0/3.D0
    DV1(1,1)=1.D0
    DV1(2,1)=0.D0
    DV1(1,2)=0.D0
    DV1(2,2)=0.D0
    DV2(1,1)=-1.D0
    DV2(2,1)= 1.D0
    DV2(1,2)=0.D0
    DV2(2,2)=0.D0
    DV3(1,1)= 1.D0
    DV3(2,1)=-1.D0
    DV3(1,2)=-1.D0
    DV3(2,2)= 1.D0

    DO I=1,MLEN
       DO J=1,MWID
          CF(J,I)=(0.D0,0.D0)
          DO NS=1,NSMAX
             CFS(J,I,NS)=(0.D0,0.D0)
          END DO
       END DO
    END DO
    DO I=1,MLEN
       CA(I)=(0.D0,0.0D0)
    END DO

! --- Vacuum contribution ---

    DO I=1,2
       DO J=1,2
          L=MCEN+3*(J-I)-1
          DTT0=DV0(I,J)
          DSS0=DV1(I,J)
          DTS1=DV2(I,J)
          DTS2=DV2(J,I)
          DTTW=DV3(I,J)

          DO NX=1,NXMAX-1
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             M=3*(NX-1)+3*(I-1)+2
             CF(L+1,M+1)=CF(L+1,M+1) &
                        +(1.D0-RNPR*RNPR)*DSS0*DX
             CF(L+1,M+2)=CF(L+1,M+2) &
                        +(1.D0-RNPR*RNPR)*DTT0*DX &
                        -DTTW/DX
             CF(L+1,M+3)=CF(L+1,M+3) &
                        +DTT0*DX &
                        -DTTW/DX
             CF(L+3,M+1)=CF(L+3,M+1) &
                        -CI*RNPR*DTS2
             CF(L-1,M+3)=CF(L-1,M+3) &
                        +CI*RNPR*DTS1
          END DO
       END DO
    END DO

! --- plasma contributioin ---

    DO I=1,2
       DO J=1,2
          L=3*(J-I+2)-1
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
                        +CD0(3,N1)*DS0(I,J,1)*DX &
                        +CD0(3,N2)*DS0(I,J,2)*DX
             CF(L+1,M+3)=CF(L+1,M+3) &
                        +CD0(4,N1)*DS0(I,J,1)*DX &
                        +CD0(4,N2)*DS0(I,J,2)*DX
             DO NS=1,NSMAX
                CFS(L+1,M+1,NS)=CFS(L+1,M+1,NS) &
                               +CM0(1,N1,NS)*DU0(I,J,1)*DX &
                               +CM0(1,N2,NS)*DU0(I,J,2)*DX
                CFS(L+2,M+1,NS)=CFS(L+2,M+1,NS) &
                               +CM0(2,N1,NS)*DT0(I,J,1)*DX &
                               +CM0(2,N2,NS)*DT0(I,J,2)*DX
                CFS(L,  M+2,NS)=CFS(L,  M+2,NS) &
                               -CM0(2,N1,NS)*DT0(J,I,1)*DX &
                               -CM0(2,N2,NS)*DT0(J,I,2)*DX
                CFS(L+1,M+2,NS)=CFS(L+1,M+2,NS) &
                               +CM0(3,N1,NS)*DS0(I,J,1)*DX &
                               +CM0(3,N2,NS)*DS0(I,J,2)*DX
                CFS(L+1,M+3,NS)=CFS(L+1,M+3,NS) &
                               +CM0(4,N1,NS)*DS0(I,J,1)*DX &
                               +CM0(4,N2,NS)*DS0(I,J,2)*DX
             END DO
          END DO
       END DO
    END DO

! --- wave guide B.C. ---

    CF(MCEN+3,1)=1.D0
    CF(MCEN  ,1)=-1.D0
    CA(1)=CFWG4

    CF(MCEN+3,2)=1.D0
    CF(MCEN  ,2)=-1.D0
    CA(1)=CFWG3

    CF(MCEN-3,4)=-CKKV
    CA(4)=-CKKV*CFWG4

    CF(MCEN-2,5)=CF(MCEN-2,5)+CI*RNPR
    CF(MCEN-3,5)=-CKKV
    CA(5)=-CKKV*CFWG3

    CF(MCEN,  MLEN-4)=1.D0    ! ER over RHS wall

    CF(MCEN+2,MLEN-3)=CKKV
    CA(MLEN-3)=CKKV*CFWG2

    CF(MCEN-2,MLEN-2)=CF(MCEN-2,MLEN-2)+CI*RNPR
    CF(MCEN+2,MLEN-2)=CKKV
    CA(MLEN-2)=CKKV*CFWG1

    CF(MCEN-2,MLEN-1)=1.D0
    CF(MCEN,MLEN-1)=-1.D0
    CA(MLEN-1)=CFWG2

    CF(MCEN-2,MLEN)=1.D0
    CF(MCEN,MLEN)=-1.D0
    CA(MLEN)=CFWG1

    CA(3*(NXANT1-1)+2+2)=CI*CFJY1/(RW*EPS0)
    CA(3*(NXANT1-1)+2+3)=CI*CFJZ1/(RW*EPS0)
    CA(3*(NXANT2-1)+2+2)=CI*CFJY2/(RW*EPS0)
    CA(3*(NXANT2-1)+2+3)=CI*CFJZ2/(RW*EPS0)

!    WRITE(6,'(A,3I5)') 'NZ,MLEN,MWID=',NZ,MLEN,MWID
    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR

!    DO NQ=1,8
!       WRITE(6,'(I5,1P2E12.3)') &
!            NQ,CA(NQ)
!    END DO
       
    DEALLOCATE(CF)
    RETURN

601 FORMAT('!! ERROR IN BANDCD : IND = ',I5)
  END SUBROUTINE W1BNDC

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWC(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NS,NX,NX1,NXD,NCL,IL,NN,MM,IA,IB,L
    REAL(rkind):: RW,RKV,RCE,DX,PABSL,RNZ
    REAL(rkind):: PIN1,PIN2,PIN3,PIN4,PIN,POUT1,POUT2,POUT3,POUT4,POUT,PCONV
    COMPLEX(rkind):: CABSL,CDEY,CDEZ,CBY,CBZ

    RW=2.D6*PI*RF
    RKV=2.D6*PI*RF/VC
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

    WRITE(6,'(A,I5,1PE12.4)') &
         '== Wave electric field == NZ,RKZ',NZ,AKZ(NZ)
    WRITE(6,'(A,1P5E12.4)') &
         'HFS-O:',CFWG4,CA(1),ABS(CA(1))**2, &
         'HFS-X:',CFWG3,CA(2),ABS(CA(2))**2, &
         'LFS-O:',CFWG2,CA(MLEN-1),ABS(CA(MLEN-1))**2, &
         'LFS-X:',CFWG1,CA(MLEN),ABS(CA(MLEN))**2
    RNZ=AKZ(NZ)*VC/RW
    IF(ABS(RNZ).LT.1.D0) THEN
       PIN1=ABS(CFWG1)**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       PIN2=ABS(CFWG2)**2*SQRT(1.D0-RNZ**2)
       PIN3=ABS(CFWG3)**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       PIN4=ABS(CFWG4)**2*SQRT(1.D0-RNZ**2)
       PIN=PIN1+PIN2+PIN3+PIN4
       POUT1=ABS(CA(MLEN))**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       POUT2=ABS(CA(MLEN-1))**2*SQRT(1.D0-RNZ**2)
       POUT3=ABS(CA(2))**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       POUT4=ABS(CA(1))**2*SQRT(1.D0-RNZ**2)
       POUT=POUT1+POUT2+POUT3+POUT4
       PCONV=PIN-POUT
    ELSE
       PIN=1.D0
       POUT=0.D0
       PCONV=PIN-POUT
    END IF
    WRITE(6,'(A,F7.2,F7.3,1P5E12.4)') &
         'REFL:',AKZ(NZ),RNZ,  &
                 POUT1/PIN,POUT2/PIN,POUT3/PIN,POUT4/PIN,PCONV/PIN

    DO NS=1,NSMAX
       DO NX=1,NXMAX-1
          DX=RKV*(XA(NX+1)-XA(NX))
          DO NX1=MAX(1,NX-1),MIN(NX+1,NXMAX-1)
             CABSL=0.D0
             NN=3*NX-1
             MM=3*NX1-1
             DO IA=1,3
                DO IB=1,3
                   L=NN+IA-(MM+IB)+MCEN
                   CABSL=CABSL &
                        +0.5D0*CONJG(CA(NN+IA))*CFS(L,MM+IB,NS)*CA(MM+IB) &
                        -0.5D0*CA(NN+IA)*CONJG(CFS(L,MM+IB,NS)*CA(MM+IB))
                END DO
             END DO
             PABSL=-CI*RCE*CABSL*RKV
!                   WRITE(6,'(A,2I5,1P5E12.4)') &
!                        'PABS:',NS,NX,RCE,CABSL,RKV,PABSL
             PABS(NX,   NS)=PABS(NX,   NS)+0.5D0*PABSL
             PABS(NX1,  NS)=PABS(NX1,  NS)+0.5D0*PABSL
          END DO
       END DO
    END DO

    DO NX=1,NXMAX
       IF(NX.EQ.1) THEN
          CDEY=(CE2DA(NZ,NX+1,2)-CE2DA(NZ,NX,2))/(XA(NX+1)-XA(NX))
          CDEZ=(CE2DA(NZ,NX+1,3)-CE2DA(NZ,NX,3))/(XA(NX+1)-XA(NX))
       ELSEIF(NX.EQ.NXMAX+1) THEN
          CDEY=(CE2DA(NZ,NX,2)-CE2DA(NZ,NX-1,2))/(XA(NX)-XA(NX-1))
          CDEZ=(CE2DA(NZ,NX,3)-CE2DA(NZ,NX-1,3))/(XA(NX)-XA(NX-1))
       ELSE
          CDEY=(CE2DA(NZ,NX+1,2)-CE2DA(NZ,NX-1,2))/(XA(NX+1)-XA(NX-1))
          CDEZ=(CE2DA(NZ,NX+1,3)-CE2DA(NZ,NX-1,3))/(XA(NX+1)-XA(NX-1))
       END IF
       CBY=-CDEZ/(-CI*RW)
       CBZ= CDEY/(-CI*RW)
       FLUX(NX)=FLUX(NX) &
            +CONJG(CE2DA(NZ,NX,2))*CBZ-CONJG(CE2DA(NZ,NX,3))*CBY
    END DO

    DEALLOCATE(CFS)
    RETURN
  END SUBROUTINE W1EPWC

!     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******

  SUBROUTINE W1CLPW(NZ)
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
  END SUBROUTINE W1CLPW

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

END MODULE w1exec9
