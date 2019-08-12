MODULE w1intg

  USE w1comm,ONLY: rkind
  USE w1fflr,ONLY: w1fnmn
  REAL(rkind):: DXD
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: CL
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NCLA
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XM,YK
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CS0,CS1,CS2,CS3
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: YX
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: SF,SG,AF,AG
  INTEGER:: NXDMIN,NXDMAX,NXWMAX,NXLMAX,NCLMAX

CONTAINS

!     ******* LOCAL PLASMA PARAMETERS *******

  SUBROUTINE W1DSPQ
    USE w1comm
    USE libdsp
    USE libgrf
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CGZ,CZ,CDZ
    REAL(rkind):: X1(4)
    COMPLEX(rkind):: CSB(3,3,4)
    INTEGER:: NLW,NCL,IL,IA,IB,NX,NX1,NS,NCMAXS,NC,NN,NXD
    REAL(rkind):: RT2,RW,FWP,FWC,FVT,RKPR,WC,UD,AKPR,ARG,RT,XD
    REAL(rkind):: VXA,VKA,VXB,VKB,VXC,VKC,VXD,VKD,VX,VK,RLI,DXA,DXB,DELTAX
    COMPLEX(rkind):: CT0A,CT1A,CT2A,CT3A,CT0B,CT1B,CT2B,CT3B
    COMPLEX(rkind):: CT0C,CT1C,CT2C,CT3C,CT0D,CT1D,CT2D,CT3D
    COMPLEX(rkind):: CT0,CT1,CT2,CT3
    
! Allocation of local data array

    IF(ALLOCATED(XM)) DEALLOCATE(XM,YX,YK,CS0,CS1,CS2,CS3)
    ALLOCATE(XM(NXPMAX),YX(NXPMAX,NSMAX),YK(NXPMAX))
    ALLOCATE(CS0(NXPMAX),CS1(NXPMAX),CS2(NXPMAX),CS3(NXPMAX))
    ALLOCATE(CGZ(NXPMAX,NCMAX),CZ(NXPMAX,NCMAX),CDZ(NXPMAX,NCMAX))

! Evaluation of NXWMAX (Maximum number of mesh to be integrated)

    NXDMIN=-1
    NXDMAX= 1
    DO NS=1,NSMAX
       FWC = AEE*PZ(NS)*BB/(AMP*PA(NS))
       FVT = AEE*1.D3/(AMP*PA(NS))
       DO NX=1,NXPMAX
          WC = FWC*PROFB(NX)
          YX(NX,NS)=WC/SQRT(FVT*PROFTP(NX,NS))   ! inverse of Larmor radius
          DO NX1=1,NXPMAX
             XD=ABS(XA(NX1)-XA(NX))*ABS(YX(NX,NS))
             IF(XD.LE.XDMAX) THEN
                NXD=NX1-NX
                IF(NXD.LT.NXDMIN) NXDMIN=NXD
                IF(NXD.GT.NXDMAX) NXDMAX=NXD
             END IF
          END DO
       END DO
    END DO
    NXDMAX=MAX(NXDMAX,-NXDMIN)
    NXDMIN=-NXDMAX
    NXWMAX=NXDMAX-NXDMIN+1
    IF(NZMAX.EQ.1) &
         WRITE(6,'(A,3I5)') 'NXDMIN,NXDMAX,NXWMAX=',NXDMIN,NXDMAX,NXWMAX

! Allocation of data index array

    IF(ALLOCATED(NCLA)) DEALLOCATE(NCLA)
    ALLOCATE(NCLA(NXWMAX,NXPMAX,NSMAX))
    NCLA(1:NXWMAX,1:NXPMAX,1:NSMAX)=0

! Evaluation of NCLMAX (Size of data array)

    NCL=0
    DO NS=1,NSMAX
       DO NX=1,NXPMAX-1
          DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXPMAX-1)
             XD=ABS(XA(NX1)-XA(NX))*ABS(YX(NX,NS))
             IF(ABS(NX1-NX).LE.1.OR.XD.LE.XDMAX) THEN
                NXD=NX1-NX+NXDMAX+1
                IF(NCLA(NXD,NX,NS).EQ.0) THEN
                   NCL=NCL+1
                   NCLA(NXD,NX,NS)=NCL
!                   WRITE(22,'(A,4I5)') 'nxd,nx,ns,ncl=',NXD,NX,NS,NCL
                END IF
             END IF
          END DO
       END DO
    END DO
    NCLMAX=NCL

    IF(NZMAX.EQ.1) WRITE(6,'(A,I5)') 'NCLMAX=',NCLMAX

! Allocation of data array

    IF(ALLOCATED(CL)) DEALLOCATE(CL)
    ALLOCATE(CL(3,3,4,NCLMAX))
    CL(1:3,1:3,1:4,1:NCLMAX)=(0.D0,0.D0)

! Calculation of CL

    RT2  = SQRT ( 2.D0 )
    RW  = 2.D6*PI*RF

    DO NS=1,NSMAX
       FWP = 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0*RW*RW)
       FWC = AEE*PZ(NS)*BB/(AMP*PA(NS))
       FVT = AEE*1.D3/(AMP*PA(NS))
       DO NX = 1 , NXPMAX
          RKPR=RKZ
          WC = FWC*PROFB(NX)
          UD = SQRT(FVT*PROFPU(NX,NS))
          AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))

          NCMAXS=2*ABS(IHARM(NS))+1
          DO NC=1,NCMAXS
             NN=NC-ABS(IHARM(NS))-1
             ARG=(RW-NN*WC)/AKPR
             CGZ(NX,NC)= ARG
          END DO
       END DO

       DO NC=1,NCMAXS
!          DO NX=1,NXPMAX
!             CALL DSPFN(CGZ(NX,NC),CZ(NX,NC),CDZ(NX,NC))
!          END DO
          CALL DSPFNA(NXPMAX,CGZ(1:NXPMAX,NC),CZ(1:NXPMAX,NC),CDZ(1:NXPMAX,NC))
       END DO

!       NX=NXPMAX/2
!      NC=3
!      WRITE(6,'(A,3I5,1P4E12.4)') 'DSP:',NS,NX,NC,CGZ(NX,NC),CZ(NX,NC)

       DO NC=1,NCMAXS
          NN=NC-ABS(IHARM(NS))-1
          DO NX=1,NXPMAX
             RKPR = RKZ
             WC=FWC*PROFB(NX)
             UD = SQRT(FVT*PROFPU(NX,NS))
             AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
             RT = PROFTP(NX,NS)/PROFTR(NX,NS)
             XM(NX)=XAM(NX)
             YK(NX)=FWP*PROFPN(NX,NS)*ABS(YX(NX,NS))*RW/AKPR
             CS0(NX)=CGZ(NX,NC)*CDZ(NX,NC)
             CS1(NX)=CZ(NX,NC)+0.5D0*(1.D0-RT)*AKPR*CDZ(NX,NC)/RW
             CS2(NX)=(RT+(1.D0-RT)*NN*WC/RW)*CDZ(NX,NC) &
                    /SQRT(2.D0*RT)
             CS3(NX)=(1.D0-1.D0/RT)*CS0(NX)*WC/RW
          END DO

          DO NX=1,NXPMAX-1
             DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXPMAX-1)
                NXD=NX1-NX+NXDMAX+1 ! positive
                NCL=NCLA(NXD,NX,NS)
                IF(NCL.NE.0) THEN
                   DXA=XA(NX+1)-XA(NX)
                   DXB=XA(NX1+1)-XA(NX1)
                   X1(1)=XA(NX)-XA(NX1)
                   X1(2)=XA(NX+1)-XA(NX1)
                   X1(3)=XA(NX)-XA(NX1+1)
                   X1(4)=XA(NX+1)-XA(NX1+1)
                   IF(NX.EQ.NX1) THEN
                      X1(1)= 0.D0
                      X1(4)= 0.D0
                   ELSEIF(NX+1.EQ.NX1) THEN
                      X1(2)=0.D0
                   ELSEIF(NX.EQ.NX1+1) THEN
                      X1(3)=0.D0
                   ENDIF

                   CALL W1QLNI(0.5D0*(XA(NX  )+XA(NX1  )), &
                               VXA,VKA,CT0A,CT1A,CT2A,CT3A,NS)
                   CALL W1QLNI(0.5D0*(XA(NX+1)+XA(NX1  )), &
                               VXB,VKB,CT0B,CT1B,CT2B,CT3B,NS)
                   CALL W1QLNI(0.5D0*(XA(NX  )+XA(NX1+1)), &
                               VXC,VKC,CT0C,CT1C,CT2C,CT3C,NS)
                   CALL W1QLNI(0.5D0*(XA(NX+1)+XA(NX1+1)), &
                               VXD,VKD,CT0D,CT1D,CT2D,CT3D,NS)
                   VX =0.25D0*(VXA +VXB +VXC +VXD )
                   VK =0.25D0*(VKA +VKB +VKC +VKD )
                   CT0=0.25D0*(CT0A+CT0B+CT0C+CT0D)
                   CT1=0.25D0*(CT1A+CT1B+CT1C+CT1D)
                   CT2=0.25D0*(CT2A+CT2B+CT2C+CT2D)
                   CT3=0.25D0*(CT3A+CT3B+CT3C+CT3D)
                   CALL W1QCAL(X1,DXA,DXB,VX,CT0,CT1,CT2,CT3,CSB,NN)

                   DO  IL=1,4
                      DO IA=1,3
                         DO IB=1,3
                            CL(IA,IB,IL,NCL)=CL(IA,IB,IL,NCL) &
                                            +VK*CSB(IA,IB,IL)*DXA*DXB
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
    END DO

!    DO NCL=600,602
!       DO IL=1,4
!          DO IB=1,3
!             WRITE(6,'(A,I4,2I2,1P6E11.3)') &
!                  'CL:',NCL,IL,IB,(CL(IA,IB,IL,NCL),IA=1,3)
!          END DO
!       END DO
!    END DO
    RETURN
  END SUBROUTINE W1DSPQ

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDQ(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2),DS1(2,2),DS2(2,2),DS3(2,2)
!  DATA DS0/0.33333 33333 33333D0,0.16666 66666 66667D0,
!     &         0.16666 66666 66667D0,0.33333 33333 33333D0/
!      DATA DS1/1.D0,0.D0,0.D0,0.D0/
!      DATA DS2/-1.D0,1.D0,0.D0,0.D0/
!      DATA DS3/1.D0,-1.D0,-1.D0,1.D0/
    INTEGER:: I,J,L,N,KML,N1,N2,M,NX,NX1,NS,IA,IB,NCL,NXD
    REAL(rkind):: RW,RKV,DTT0,DSS0,DTS1,DTS2,DTTW,DX
    REAL(rkind):: RKPR,RNPR

    DS0(1,1)=1.D0/3.D0
    DS0(2,1)=1.D0/6.D0
    DS0(1,2)=1.D0/6.D0
    DS0(2,2)=1.D0/3.D0
    DS1(1,1)=1.D0
    DS1(2,1)=0.D0
    DS1(1,2)=0.D0
    DS1(2,2)=0.D0
    DS2(1,1)=-1.D0
    DS2(2,1)= 1.D0
    DS2(1,2)=0.D0
    DS2(2,2)=0.D0
    DS3(1,1)= 1.D0
    DS3(2,1)=-1.D0
    DS3(1,2)=-1.D0
    DS3(2,2)= 1.D0

    MWID=3*(2*NXDMAX+4)-1
    MLEN=3*NXPMAX+4
    ALLOCATE(CF(MWID,MLEN))

    RW=2.D6*PI*RF
    RKV=RW/VC

    DO I=1,MLEN
       DO J=1,MWID
          CF(J,I)=(0.D0,0.D0)
       END DO
    END DO
    DO I=1,MLEN
       CA(I)=(0.D0,0.0D0)
    END DO

    KML=3*NXDMAX
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

    CF(KML+3,MLEN-4)=-1.D0
    CF(KML+6,MLEN-4)= 1.D0

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
          L=3*(J-I+NXDMAX+2)-1
          DTT0=DS0(I,J)
          DSS0=DS1(I,J)
          DTS1=DS2(I,J)
          DTS2=DS2(J,I)
          DTTW=DS3(I,J)

          DO NX=1,NXPMAX-1
             RKPR=RKZ
             RNPR=VC*RKPR/RW
             N1=NX
             N2=NX+1
             DX=RKV*(XA(N2)-XA(N1))
             M=3*(NX+I-1)-1
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

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          M=3*(NX-1)+2
          DO NXD=1,NXWMAX
             L=3*NXD+2
             NCL=NCLA(NXD,NX,NS)
             IF(NCL.NE.0) THEN
                DO IA=1,3
                   DO IB=1,3
                      CF(L+IA-IB+1,M+IB  )=CF(L+IA-IB+1,M+IB  ) &
                                          +CL(IB,IA,1,NCL)*RKV
                      CF(L+IA-IB+4,M+IB  )=CF(L+IA-IB+4,M+IB  ) &
                                          +CL(IB,IA,3,NCL)*RKV
                      CF(L+IA-IB-2,M+IB+3)=CF(L+IA-IB-2,M+IB+3) &
                                          +CL(IB,IA,2,NCL)*RKV
                      CF(L+IA-IB+1,M+IB+3)=CF(L+IA-IB+1,M+IB+3) &
                                          +CL(IB,IA,4,NCL)*RKV
                   END DO
                END DO
             END IF
          END DO
       END DO
    END DO

!    DO I=MLEN/2-3,MLEN/2+3
!       DO J=1,15,3
!          WRITE(6,'(A,2I4,1P6E11.3)') &
!               'CF:',I,J,CF(J,I),CF(J+1,I),CF(J+2,I)
!       END DO
!       J=16
!          WRITE(6,'(A,2I4,1P4E11.3)') &
!               'CF:',I,J,CF(J,I),CF(J+1,I)
!    END DO

    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

601 FORMAT('!! ERROR IN BANDCD : IND = ',I5)
  END SUBROUTINE W1BNDQ

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWQ(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NS,NX,NX1,NXD,NCL,IL,NN,MM,IA,IB
    REAL(rkind):: RKV,RCE,DX,PABSL
    COMPLEX(rkind):: CABSL

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

    DO NS=1,NSMAX
       DO NX=1,NXPMAX-1
          DX=RKV*(XA(NX+1)-XA(NX))
          DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXPMAX-1)
             NXD=NX1-NX+NXDMAX+1 ! positive
             NCL=NCLA(NXD,NX,NS)
             IF(NCL.NE.0) THEN
                DO IL=1,4
                   CABSL=0.D0
                   NN=3*NX-1
                   MM=3*NX1-1
                   IF(IL.EQ.2) NN=NN+3
                   IF(IL.EQ.3) MM=MM+3
                   IF(IL.EQ.4) NN=NN+3
                   IF(IL.EQ.4) MM=MM+3
                   DO IA=1,3
                      DO IB=1,3
                         CABSL=CABSL &
                         -0.5D0*CONJG(CA(NN+IA))*CL(IA,IB,IL,NCL)*CA(MM+IB) &
                         -0.5D0*CA(NN+IA)*CONJG(CL(IB,IA,IL,NCL)*CA(MM+IB))
                      END DO
                   END DO
                   PABSL=-CI*RCE*CABSL*RKV
!                   WRITE(6,'(A,2I5,1P5E12.4)') &
!                        'PABS:',NS,NX,RCE,CABSL,RKV,PABSL
                   IF(IL.EQ.1) THEN
                      PABS(NX,   NS)=PABS(NX,   NS)+0.5D0*PABSL
                      PABS(NX1,  NS)=PABS(NX1,  NS)+0.5D0*PABSL
                   ELSEIF(IL.EQ.2) THEN
                      PABS(NX+1, NS)=PABS(NX+1, NS)+0.5D0*PABSL
                      PABS(NX1,  NS)=PABS(NX1,  NS)+0.5D0*PABSL
                   ELSEIF(IL.EQ.3) THEN
                      PABS(NX,   NS)=PABS(NX,   NS)+0.5D0*PABSL
                      PABS(NX1+1,NS)=PABS(NX1+1,NS)+0.5D0*PABSL
                   ELSEIF(IL.EQ.4) THEN
                      PABS(NX+1, NS)=PABS(NX+1, NS)+0.5D0*PABSL
                      PABS(NX1+1,NS)=PABS(NX1+1,NS)+0.5D0*PABSL
                   END IF
                END DO
             END IF
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE W1EPWQ

!      ****** MAKE TABLE OF F,G,H ******

  SUBROUTINE W1QTBL
    USE w1comm
    IMPLICIT NONE
    INTEGER,SAVE:: NDMAX_save=0
    INTEGER,SAVE:: NCMAX_save=0
    REAL(rkind),SAVE:: XDMAX_save=0.D0
    INTEGER:: NC,NN,ND
    
    IF(NDMAX.EQ.NDMAX_save.AND. &
       NCMAX.EQ.NCMAX_save.AND. &
       ABS(XDMAX-XDMAX_save).LE.1.D-32) RETURN

    IF(ALLOCATED(SF)) DEALLOCATE(SF,SG,AF,AG)
    ALLOCATE(SF(NDMAX+1,NCMAX+1,3),SG(NDMAX+1,NCMAX+1,3))
    ALLOCATE(AF(NDMAX+1,NCMAX+1,2),AG(NDMAX+1,NCMAX+1,2))

    DXD=XDMAX/NDMAX
    DO NC=1,NCMAX+1
       NN=NC-1
       DO ND=1,NDMAX+1
          SF(ND,NC,1)= W1FNMN((ND-1)*DXD,1,NN,1)
          SF(ND,NC,2)= W1FNMN((ND-1)*DXD,1,NN,3)
          SF(ND,NC,3)= W1FNMN((ND-1)*DXD,1,NN,5)
          SG(ND,NC,1)= W1FNMN((ND-1)*DXD,2,NN,0)
          SG(ND,NC,2)= W1FNMN((ND-1)*DXD,2,NN,2)
          SG(ND,NC,3)= W1FNMN((ND-1)*DXD,2,NN,4)
          AF(ND,NC,1)=-W1FNMN((ND-1)*DXD,3,NN,1)
          AF(ND,NC,2)=-W1FNMN((ND-1)*DXD,3,NN,3)
          AG(ND,NC,1)=-W1FNMN((ND-1)*DXD,4,NN,0)
          AG(ND,NC,2)=-W1FNMN((ND-1)*DXD,4,NN,2)
       END DO
    END DO
    NDMAX_save=NDMAX
    NCMAX_save=NCMAX
    XDMAX_save=XDMAX
    RETURN
  END SUBROUTINE W1QTBL

!     ******** LINEAR INTERPOLATION  ******************************

  SUBROUTINE W1QCAL(GX,DXA,DXB,YX,CT0,CT1,CT2,CT3,CSB,NN)
    USE w1comm,ONLY: rkind,XDMAX,NDMAX
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NN
    REAL(rkind),INTENT(IN):: GX(4),DXA,DXB,YX
    COMPLEX(rkind),INTENT(IN):: CT0,CT1,CT2,CT3
    COMPLEX(rkind),INTENT(OUT):: CSB(3,3,4)
    REAL(rkind):: AF0(4,3),AG0(4,3), &
                  AF1(4,2),AF2(4,2),AF3(4,1),AF4(4,1), &
                  AG1(4,2),AG2(4,2),AG3(4,1),AG4(4,1), &
                  AAF(4),AAG(4),ADF(4),ADG(4),AQF(4),AQG(4),AQH(4), &
                  ABF(4),ACF(4),ABG(4),BDG(4),GV(4)
    INTEGER:: NNA,NX,IV,IR,NPV,NPV1,NPV2,IU
    REAL(rkind):: ESGN,VSGN,YSGN,V,PV,PVS,PVT,EI,EJ,E0,E1,E2,E3,E4
    REAL(rkind):: FA0,FA1,FA2,FA3,FA4,GA0,GA1,GA2,GA3,GA4,HA0,HA1,HA2
    REAL(rkind):: FF,FFS,FFT,FFST,FG,FGS,FGT,FGST
    REAL(rkind):: DF,DFS,DFT,DFST,DG,DGS,DGT,DGST
    REAL(rkind):: QF,QFS,QFT,QFST,QG,QGS,QGT,QGST,QH,QHS,QHT,QHST

    NNA=ABS(NN)+1
    ESGN=SIGN(1.D0,DBLE(NN))

    DO IV=1,4
       GV(IV)=YX*GX(IV)
       VSGN=SIGN(1.D0,GV(IV))
       YSGN=SIGN(1.D0,YX)
       V=ABS(GV(IV))
       IF(V.GT.XDMAX) THEN
          DO IR=1,3
             AF0(IV,IR)=0.D0
             AG0(IV,IR)=0.D0
          END DO
          AF1(IV,1)=0.D0
          AG1(IV,1)=0.D0
          AF1(IV,2)=0.D0
          AG1(IV,2)=0.D0
       ELSEIF(ABS(GX(IV)).LT.1.D-70) THEN
          AF0(IV,1)= SF(1,NNA,1)
          AF0(IV,2)= SF(1,NNA,2)
          AF0(IV,3)= SF(1,NNA,3)
          AG0(IV,1)= SG(1,NNA,1)*ESGN
          AG0(IV,2)= SG(1,NNA,2)*ESGN
          AG0(IV,3)= SG(1,NNA,3)*ESGN
          IF(IV.EQ.1.OR.IV.EQ.3) THEN
             AF1(IV,1)= AF(1,NNA,1)*YSGN
             AF1(IV,2)= AF(1,NNA,2)*YSGN
             AG1(IV,1)= AG(1,NNA,1)*YSGN*ESGN
             AG1(IV,2)= AG(1,NNA,2)*YSGN*ESGN
          ELSE
             AF1(IV,1)=-AF(1,NNA,1)*YSGN
             AF1(IV,2)=-AF(1,NNA,2)*YSGN
             AG1(IV,1)=-AG(1,NNA,1)*YSGN*ESGN
             AG1(IV,2)=-AG(1,NNA,2)*YSGN*ESGN
          ENDIF
       ELSE
          PV=V/DXD
          NPV=INT(PV)
          IF(NPV.GE.NDMAX) THEN
             NPV1=NDMAX+1
             NPV2=NDMAX+1
          ELSE
             NPV1=NPV+1
             NPV2=NPV+2
          ENDIF
          PVS=PV-DBLE(NPV)
          PVT=1.D0-PVS
          AF0(IV,1)= PVT*SF(NPV1,NNA,1)+PVS*SF(NPV2,NNA,1)
          AF0(IV,2)= PVT*SF(NPV1,NNA,2)+PVS*SF(NPV2,NNA,2)
          AF0(IV,3)= PVT*SF(NPV1,NNA,3)+PVS*SF(NPV2,NNA,3)
          AG0(IV,1)=(PVT*SG(NPV1,NNA,1)+PVS*SG(NPV2,NNA,1))*ESGN
          AG0(IV,2)=(PVT*SG(NPV1,NNA,2)+PVS*SG(NPV2,NNA,2))*ESGN
          AG0(IV,3)=(PVT*SG(NPV1,NNA,3)+PVS*SG(NPV2,NNA,3))*ESGN
          AF1(IV,1)=(PVT*AF(NPV1,NNA,1)+PVS*AF(NPV2,NNA,1))*VSGN
          AF1(IV,2)=(PVT*AF(NPV1,NNA,2)+PVS*AF(NPV2,NNA,2))*VSGN
          AG1(IV,1)=(PVT*AG(NPV1,NNA,1)+PVS*AG(NPV2,NNA,1))*VSGN*ESGN
          AG1(IV,2)=(PVT*AG(NPV1,NNA,2)+PVS*AG(NPV2,NNA,2))*VSGN*ESGN
       ENDIF
       AF2(IV,1)= 4.D0*AF0(IV,2)+GV(IV)*AF1(IV,1)
       AF2(IV,2)= 4.D0*AF0(IV,3)+GV(IV)*AF1(IV,2)
       AF3(IV,1)=(4.D0*AF1(IV,2)+GV(IV)*AF2(IV,1))/2.D0
       AF4(IV,1)=(4.D0*AF2(IV,2)+GV(IV)*AF3(IV,1))/3.D0
       AG2(IV,1)= 4.D0*AG0(IV,2)+GV(IV)*AG1(IV,1)
       AG2(IV,2)= 4.D0*AG0(IV,3)+GV(IV)*AG1(IV,2)
       AG3(IV,1)=(4.D0*AG1(IV,2)+GV(IV)*AG2(IV,1))/2.D0
       AG4(IV,1)=(4.D0*AG2(IV,2)+GV(IV)*AG3(IV,1))/3.D0
    END DO

    EI=YX*DXA
    EJ=YX*DXB
    E0=EI
    E1=EI*EJ
    E2=EI*EI*EJ
    E3=EI*EJ*EJ
    E4=EI*EI*EJ*EJ

    FA0=AF0(1,1)-AF0(2,1)-AF0(3,1)+AF0(4,1)
    FA1=AF1(1,1)-AF1(2,1)-AF1(3,1)+AF1(4,1)
    FA2=AF2(1,1)-AF2(2,1)-AF2(3,1)+AF2(4,1)
    FA3=AF3(1,1)-AF3(2,1)-AF3(3,1)+AF3(4,1)
    FA4=AF4(1,1)-AF4(2,1)-AF4(3,1)+AF4(4,1)
    GA0=AG0(1,1)-AG0(2,1)-AG0(3,1)+AG0(4,1)
    GA1=AG1(1,1)-AG1(2,1)-AG1(3,1)+AG1(4,1)
    GA2=AG2(1,1)-AG2(2,1)-AG2(3,1)+AG2(4,1)
    GA3=AG3(1,1)-AG3(2,1)-AG3(3,1)+AG3(4,1)
    GA4=AG4(1,1)-AG4(2,1)-AG4(3,1)+AG4(4,1)
    HA0=AF0(1,2)-AF0(2,2)-AF0(3,2)+AF0(4,2)
    HA1=AF1(1,2)-AF1(2,2)-AF1(3,2)+AF1(4,2)
    HA2=AF2(1,2)-AF2(2,2)-AF2(3,2)+AF2(4,2)

    FF  =-FA2/E1
    FFS = FA3/E3-(AF2(4,1)-AF2(2,1))/E1
    FFT =-FA3/E2-(AF2(4,1)-AF2(3,1))/E1
    FFST= FA4/E4-(AF3(4,1)-AF3(2,1))/E3 &
                +(AF3(4,1)-AF3(3,1))/E2-AF2(4,1)/E1
    FG  =-GA2/E1
    FGS = GA3/E3-(AG2(4,1)-AG2(2,1))/E1
    FGT =-GA3/E2-(AG2(4,1)-AG2(3,1))/E1
    FGST= GA4/E4-(AG3(4,1)-AG3(2,1))/E3 &
                +(AG3(4,1)-AG3(3,1))/E2-AG2(4,1)/E1
    DF  =-FA1/E1
    DFS = FA2/E3-(AF1(4,1)-AF1(2,1))/E1
    DFT =-FA2/E2-(AF1(4,1)-AF1(3,1))/E1
    DFST= FA3/E4-(AF2(4,1)-AF2(2,1))/E3 &
                +(AF2(4,1)-AF2(3,1))/E2-AF1(4,1)/E1
    DG  =-GA1/E1
    DGS = GA2/E3-(AG1(4,1)-AG1(2,1))/E1
    DGT =-GA2/E2-(AG1(4,1)-AG1(3,1))/E1
    DGST= GA3/E4-(AG2(4,1)-AG2(2,1))/E3 &
                +(AG2(4,1)-AG2(3,1))/E2-AG1(4,1)/E1
    QF  =-FA0/E1
    QFS = FA1/E3-(AF0(4,1)-AF0(2,1))/E1
    QFT =-FA1/E2-(AF0(4,1)-AF0(3,1))/E1
    QFST= FA2/E4-(AF1(4,1)-AF1(2,1))/E3 &
                +(AF1(4,1)-AF1(3,1))/E2-AF0(4,1)/E1
    QG  =-GA0/E1
    QGS = GA1/E3-(AG0(4,1)-AG0(2,1))/E1
    QGT =-GA1/E2-(AG0(4,1)-AG0(3,1))/E1
    QGST= GA2/E4-(AG1(4,1)-AG1(2,1))/E3 &
                +(AG1(4,1)-AG1(3,1))/E2-AG0(4,1)/E1
    QH  =-HA0/E1
    QHS = HA1/E3-(AF0(4,2)-AF0(2,2))/E1
    QHT =-HA1/E2-(AF0(4,2)-AF0(3,2))/E1
    QHST= HA2/E4-(AF1(4,2)-AF1(2,2))/E3 &
                +(AF1(4,2)-AF1(3,2))/E2-AF0(4,2)/E1

    IF(ABS(GX(1)).LT.1.D-70) THEN
       FF  =FF  -(AF1(1,1)-AF1(4,1))/E0
       FFS =FFS -(AF1(1,1)-AF1(4,1))/(2.D0*E0)
       FFT =FFT -(AF1(1,1)-AF1(4,1))/(2.D0*E0) &
                -(AF2(1,1)-AF2(4,1))/E1
       FFST=FFST-(AF1(1,1)-AF1(4,1))/(3.D0*E0) &
                -(AF2(1,1)-AF2(4,1))/(2.D0*E1)
       FG  =FG  -(AG1(1,1)-AG1(4,1))/E0
       FGS =FGS -(AG1(1,1)-AG1(4,1))/(2.D0*E0)
       FGT =FGT -(AG1(1,1)-AG1(4,1))/(2.D0*E0) &
                -(AG2(1,1)-AG2(4,1))/E1
       FGST=FGST-(AG1(1,1)-AG1(4,1))/(3.D0*E0) &
                -(AG2(1,1)-AG2(4,1))/(2.D0*E1)
       DFT =DFT -(AF1(1,1)-AF1(4,1))/E1
       DFST=DFST-(AF1(1,1)-AF1(4,1))/(2.D0*E1)
       DGT =DGT -(AG1(1,1)-AG1(4,1))/E1
       DGST=DGST-(AG1(1,1)-AG1(4,1))/(2.D0*E1)
    ENDIF

    AAF(1)=FF-FFS-FFT+FFST
    AAF(2)=FFS-FFST
    AAF(3)=FFT-FFST
    AAF(4)=FFST
    ABF(1)=FF-FFS
    ABF(2)=FFS
    ABF(3)=0.D0
    ABF(4)=0.D0
    ACF(1)=FF-FFT
    ACF(2)=0.D0
    ACF(3)=FFT
    ACF(4)=0.D0
    AAG(1)=FG-FGS-FGT+FGST
    AAG(2)=FGS-FGST
    AAG(3)=FGT-FGST
    AAG(4)=FGST
    ABG(1)=FG
    ABG(2)=0.D0
    ABG(3)=0.D0
    ABG(4)=0.D0
    ADF(1)=DF-DFS-DFT+DFST
    ADF(2)=DFS-DFST
    ADF(3)=DFT-DFST
    ADF(4)=DFST
    ADG(1)=DG-DGS
    ADG(2)=DGS
    ADG(3)=0.D0
    ADG(4)=0.D0
    BDG(1)=DG-DGT
    BDG(2)=0.D0
    BDG(3)=DGT
    BDG(4)=0.D0
    AQF(1)=QF-QFS-QFT+QFST
    AQF(2)=QFS-QFST
    AQF(3)=QFT-QFST
    AQF(4)=QFST
    AQG(1)=QG-QGS-QGT+QGST
    AQG(2)=QGS-QGST
    AQG(3)=QGT-QGST
    AQG(4)=QGST
    AQH(1)=QH-QHS-QHT+QHST
    AQH(2)=QHS-QHST
    AQH(3)=QHT-QHST
    AQH(4)=QHST

    DO IU=1,4
       CSB(1,1,IU)=NN*ABG(IU)*CT1
       CSB(1,2,IU)=(0.D0,-1.D0)*NN*ACF(IU)*CT1
       CSB(1,3,IU)= BDG(IU)*CT2
       CSB(2,1,IU)=(0.D0,1.D0)*NN*ABF(IU)*CT1
       CSB(2,2,IU)=(NN*AAG(IU)-2.D0*AQF(IU))*CT1
       CSB(2,3,IU)=-(0.D0,1.D0)*ADF(IU)*CT2
       CSB(3,1,IU)= ADG(IU)*CT2
       CSB(3,2,IU)= (0.D0,1.D0)*ADF(IU)*CT2
       CSB(3,3,IU)=-CT0*(AAF(IU)+NN*AAG(IU)-2.D0*AQF(IU)+2.D0*AQH(IU)) &
                   -CT3*AQG(IU)
    END DO
    RETURN
  END SUBROUTINE W1QCAL

!     ******** LINEAR INTERPOLATION 2 ******************************

  SUBROUTINE W1QLNI(V,VX,VK,CT0,CT1,CT2,CT3,NS)
    USE w1comm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: V
    REAL(rkind),INTENT(OUT):: VX,VK
    COMPLEX(rkind),INTENT(OUT):: CT0,CT1,CT2,CT3
    INTEGER:: NS,NHF,I,NPV
    REAL(rkind):: PVS

    IF(V.LE.XM(1)) THEN
       VX=YX(1,NS)
       VK=YK(1)
       CT0=CS0(1)
       CT1=CS1(1)
       CT2=CS2(1)
       CT3=CS3(1)
    ELSEIF(V.GE.XM(NXPMAX)) THEN
       VX=YX(NXPMAX,NS)
       VK=YK(NXPMAX)
       CT0=CS0(NXPMAX)
       CT1=CS1(NXPMAX)
       CT2=CS2(NXPMAX)
       CT3=CS3(NXPMAX)
    ELSE
       NHF=INT(NXPMAX*(V-XM(1))/(XM(NXPMAX)-XM(1)))
       IF(NHF.LT.1)      NHF=1
       IF(NHF.GT.NXPMAX) NHF=NXPMAX
       IF(V.LT.XM(NHF)) THEN
          DO I=NHF-1,1,-1
             IF(XM(I).LT.V) THEN
                NPV=I
                GOTO 100
             ENDIF
          END DO
       ELSE
          DO I=NHF+1,NXPMAX
             IF(XM(I).GT.V) THEN
                NPV=I-1
                GOTO 100
             ENDIF
          END DO
       ENDIF
100    CONTINUE
       PVS=(V-XM(NPV))/(XM(NPV+1)-XM(NPV))
       VX  =(1.D0-PVS)*YX(NPV,NS) +PVS*YX(NPV+1,NS)
       VK  =(1.D0-PVS)*YK(NPV) +PVS*YK(NPV+1)
       CT0 =(1.D0-PVS)*CS0(NPV)+PVS*CS0(NPV+1)
       CT1 =(1.D0-PVS)*CS1(NPV)+PVS*CS1(NPV+1)
       CT2 =(1.D0-PVS)*CS2(NPV)+PVS*CS2(NPV+1)
       CT3 =(1.D0-PVS)*CS3(NPV)+PVS*CS3(NPV+1)
    ENDIF
    RETURN
  END SUBROUTINE W1QLNI
END MODULE w1intg
