! w1exec7.f90

MODULE w1exec7

  USE w1comm,ONLY: rkind
  USE w1fflr,ONLY: w1fnmn

  PRIVATE
  PUBLIC w1_exec7,w1qtblx

  REAL(rkind):: DXD
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: CL
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NCLA
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XM,YK
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CS0,CS1,CS2,CS3
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: YX
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: SF,SG,AF,AG
  INTEGER:: NXDMIN,NXDMAX,NXWMAX,NXLMAX,NCLMAX

CONTAINS

  SUBROUTINE w1_exec7(NZ,IERR)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    CALL W1BCND
    CALL W1DSPQ
    CALL W1BNDQ(NZ,IERR)
       IF(IERR.NE.0) RETURN
    CALL W1EPWQ(NZ)
    CALL W1CLCD(NZ)
    CALL W1CLPW(NZ)
    RETURN
  END SUBROUTINE w1_exec7

!     ****** SET BOUNDARY CONDITIONS AT R=RA ******

  SUBROUTINE W1BCND
    USE w1comm
    IMPLICIT NONE
    INTEGER:: I,J
    REAL(rkind):: RW,RKV,RKPR,RNPR,RCE
    COMPLEX(rkind):: CKKV2,CKKV

    RW=2.D6*PI*RF
    RKV=RW/VC
    RKPR=RKZ
    RNPR=VC*RKPR/RW
    RCE=VC*EPS0
    CKKV2=RNPR*RNPR-1.D0
    CKKV =SQRT(CKKV2)

    IF(ABS(WALLR).GT.1.D-12) THEN
       CKKV2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
       CKKV =SQRT(CKKV2)
    ELSE 
       CKKV2=RNPR*RNPR-1.D0
       CKKV =SQRT(CKKV2)
    ENDIF

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

!     ******* LOCAL PLASMA PARAMETERS for INTG *******

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
    ALLOCATE(XM(NXMAX),YX(NXMAX,NSMAX),YK(NXMAX))
    ALLOCATE(CS0(NXMAX),CS1(NXMAX),CS2(NXMAX),CS3(NXMAX))
    ALLOCATE(CGZ(NXMAX,NCMAX),CZ(NXMAX,NCMAX),CDZ(NXMAX,NCMAX))

! Evaluation of NXWMAX (Maximum number of mesh to be integrated)

    NXDMIN=-1
    NXDMAX= 1
    DO NS=1,NSMAX
       FWC = AEE*PZ(NS)*BB/(AMP*PA(NS))
       FVT = AEE*1.D3/(AMP*PA(NS))
       DO NX=1,NXMAX
          WC = FWC*PROFB(NX)
          YX(NX,NS)=WC/SQRT(FVT*PROFTP(NX,NS))   ! inverse of Larmor radius
          DO NX1=1,NXMAX
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
    ALLOCATE(NCLA(NXWMAX,NXMAX,NSMAX))
    NCLA(1:NXWMAX,1:NXMAX,1:NSMAX)=0

! Evaluation of NCLMAX (Size of data array)

    NCL=0
    DO NS=1,NSMAX
       DO NX=1,NXMAX-1
          DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXMAX)
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

    IF(NZMAX.EQ.1) WRITE(6,'(A,I8)') 'NCLMAX=',NCLMAX

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
       DO NX = 1 , NXMAX
          RKPR=RKZ
          IF(ABS(RKPR).LE.1.D-6) RKPR=1.D-6
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
!          DO NX=1,NXMAX
!             CALL DSPFN(CGZ(NX,NC),CZ(NX,NC),CDZ(NX,NC))
!          END DO
          CALL DSPFNA(NXMAX,CGZ(1:NXMAX,NC),CZ(1:NXMAX,NC), &
                            CDZ(1:NXMAX,NC))
       END DO

!       NX=NXMAX/2
!      NC=3
!      WRITE(6,'(A,3I5,1P4E12.4)') 'DSP:',NS,NX,NC,CGZ(NX,NC),CZ(NX,NC)

       DO NC=1,NCMAXS
          NN=NC-ABS(IHARM(NS))-1
          DO NX=1,NXMAX
             RKPR = RKZ
             IF(ABS(RKPR).LE.1.D-6) RKPR=1.D-6
             WC=FWC*PROFB(NX)
             UD = SQRT(FVT*PROFPU(NX,NS))
             AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
             RT = PROFTP(NX,NS)/PROFTR(NX,NS)
             XM(NX)=XAM(NX)
!             XM(NX)=XA(NX)
             YK(NX)=FWP*PROFPN(NX,NS)*ABS(YX(NX,NS))*RW/AKPR
             CS0(NX)=CGZ(NX,NC)*CDZ(NX,NC)
             CS1(NX)=CZ(NX,NC)+0.5D0*(1.D0-RT)*AKPR*CDZ(NX,NC)/RW
             CS2(NX)=(RT+(1.D0-RT)*NN*WC/RW)*CDZ(NX,NC) &
                    /SQRT(2.D0*RT)
             CS3(NX)=(1.D0-1.D0/RT)*CS0(NX)*WC/RW
          END DO

          DO NX=1,NXMAX-1
             DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXMAX-1)
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

  SUBROUTINE W1BNDQ(NZ,IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: DS0(2,2),DS1(2,2),DS2(2,2),DS3(2,2)
    INTEGER:: I,J,L,N,MCEN,N1,N2,M,NX,NX1,NS,IA,IB,NCL,NXD,NQ
    REAL(rkind):: RW,RKV,DTT0,DSS0,DTS1,DTS2,DTTW,DX
    REAL(rkind):: RKPR,RNPR,FACT,RKZA
    COMPLEX(rkind):: CFJX1,CFJX2

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
    MLEN=3*NXMAX+4
    ALLOCATE(CF(MWID,MLEN))
    MCEN=3*NXDMAX+6

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

    DO I=1,2
       DO J=1,2
          L=MCEN+3*(J-I)-1
          DTT0=DS0(I,J)
          DSS0=DS1(I,J)
          DTS1=DS2(I,J)
          DTS2=DS2(J,I)
          DTTW=DS3(I,J)

          DO NX=1,NXMAX-1
             RKPR=RKZ
             IF(ABS(RKPR).LE.1.D-6) RKPR=1.D-6
             RNPR=VC*RKPR/RW
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

    
!    WRITE(6,*) 'MCEN=',MCEN
!    DO NQ=3*NXMAX/2+2-4,3*NXMAX/2+2+4
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN-3,NQ),CF(MCEN-2,NQ),CF(MCEN-1,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN  ,NQ),CF(MCEN+1,NQ),CF(MCEN+2,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN+3,NQ),CF(MCEN+4,NQ),CF(MCEN+5,NQ)
!    END DO

!    NS=1
!    NX=NXMAX/2
!    NXD=NXWMAX/2
!    NCL=NCLA(NXD,NX,NS)
!    WRITE(21,'(A,4I5)') 'NS,NX,NXD,NCL=',NS,NX,NXD,NCL
!    DO IA=1,3
!       DO IB=1,3
!          WRITE(21,'(I5,3I3,1P2E12.4)') NCL,IA,IA,1,CL(IB,IA,1,NCL)
!          WRITE(21,'(I5,3I3,1P2E12.4)') NCL,IA,IA,2,CL(IB,IA,2,NCL)
!          WRITE(21,'(I5,3I3,1P2E12.4)') NCL,IA,IA,3,CL(IB,IA,3,NCL)
!          WRITE(21,'(I5,3I3,1P2E12.4)')NCL,IA,IA,4,CL(IB,IA,4,NCL)
!       END DO
!    END DO

    DO NS=1,NSMAX
       DO NX=1,NXMAX
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

!    WRITE(6,*) 'MCEN=',MCEN
!    DO NQ=3*NXMAX/2+2-4,3*NXMAX/2+2+4
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN-3,NQ),CF(MCEN-2,NQ),CF(MCEN-1,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN  ,NQ),CF(MCEN+1,NQ),CF(MCEN+2,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN+3,NQ),CF(MCEN+4,NQ),CF(MCEN+5,NQ)
!    END DO

    IF(MDLWG.EQ.4) THEN
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

!   antenna current in HFS

    RKZA=RKZ
    IF(ABS(RKZA).LE.1.D-6) RKZA=1.D-6

    IF(NXANT1.GE.1) THEN
       CFJX1=CI*RKZA*CFJZ1
       DO NX=1,NXANT1-1
          DX=XA(NX+1)-XA(NX)
          CA(3*(NX-1)+2+1)=CI*CFJX1*DX/(RW*EPS0)
       END DO
       DX=-RD-XA(NXANT1)
       FACT=(-RD-XA(NXANT1))/(XA(NXANT1+1)-XA(NXANT1))
       CA(3*(NXANT1-1)+2+1)=            CI*CFJX1*DX/(RW*EPS0)
       CA(3*(NXANT1-1)+2+2)=(1.D0-FACT)*CI*CFJY1/(RW*EPS0)
       CA(3* NXANT1   +2+2)= FACT      *CI*CFJY1/(RW*EPS0)
       CA(3*(NXANT1-1)+2+3)=(1.D0-FACT)*CI*CFJZ1/(RW*EPS0)
       CA(3* NXANT1   +2+3)= FACT      *CI*CFJZ1/(RW*EPS0)
    END IF

!   antenna current in LFS

    IF(NXANT2.GE.0) THEN
       CFJX2=CI*RKZA*CFJZ2
       DO NX=NXANT2+1,NXMAX-1
          DX=XA(NX+1)-XA(NX)
          CA(3*(NX-1)+2+1)=CI*CFJX2*DX/(RW*EPS0)
       END DO
       DX=XA(NXANT2+1)-RD
       FACT=( RD-XA(NXANT2))/(XA(NXANT2+1)-XA(NXANT2))
       CA(3*(NXANT2-1)+2+1)=            CI*CFJX2*DX/(RW*EPS0)
       CA(3*(NXANT2-1)+2+2)=(1.D0-FACT)*CI*CFJY2/(RW*EPS0)
       CA(3* NXANT2   +2+2)= FACT      *CI*CFJY2/(RW*EPS0)
       CA(3*(NXANT2-1)+2+3)=(1.D0-FACT)*CI*CFJZ2/(RW*EPS0)
       CA(3* NXANT2   +2+3)= FACT      *CI*CFJZ2/(RW*EPS0)
    END IF

!    DO I=MLEN/2-3,MLEN/2+3
!       DO J=1,15,3
!          WRITE(6,'(A,2I4,1P6E11.3)') &
!               'CF:',I,J,CF(J,I),CF(J+1,I),CF(J+2,I)
!       END DO
!       J=16
!          WRITE(6,'(A,2I4,1P4E11.3)') &
!               'CF:',I,J,CF(J,I),CF(J+1,I)
!    END DO

!    DO NQ=1,8
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN-3,NQ),CF(MCEN-2,NQ),CF(MCEN-1,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN  ,NQ),CF(MCEN+1,NQ),CF(MCEN+2,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN+3,NQ),CF(MCEN+4,NQ),CF(MCEN+5,NQ)
!    END DO

!    DO NQ=MLEN-7,MLEN
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN-3,NQ),CF(MCEN-2,NQ),CF(MCEN-1,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN  ,NQ),CF(MCEN+1,NQ),CF(MCEN+2,NQ)
!       WRITE(6,'(I5,1P6E12.3)') &
!            NQ,CF(MCEN+3,NQ),CF(MCEN+4,NQ),CF(MCEN+5,NQ)
!    END DO

!    DO N=1,MLEN
!       IF(ABS(CA(N)).NE.0.D0) THEN
!          WRITE(6,'(A,I5,1P2E12.3)') 'I,CA=',N,CA(N)
!       END IF
!    END DO

    WRITE(6,'(A,3I5)') 'NZ,MLEN,MWID=',NZ,MLEN,MWID
    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR

!    DO NQ=1,8
!       WRITE(6,'(I5,1P2E12.3)') &
!            NQ,CA(NQ)
!    END DO
       
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
    REAL(rkind):: RW,RKV,RCE,DX,PABSL
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

    DO NS=1,NSMAX
       DO NX=1,NXMAX-1
          DX=RKV*(XA(NX+1)-XA(NX))
          DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXMAX-1)
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
!                         -CONJG(CA(NN+IA))*CL(IA,IB,IL,NCL)*CA(MM+IB)
                      END DO
                   END DO
                   PABSL=-CI*RCE*CABSL*RKV
!                   WRITE(6,'(A,2I5,1P5E12.4)') &
!                        'PABS:',NS,NX,RCE,CABSL,RKV,PABSL
                   PABS(NX,   NS)=PABS(NX,   NS)+0.5D0*PABSL
                   PABS(NX1,  NS)=PABS(NX1,  NS)+0.5D0*PABSL
                END DO
             END IF
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
    RETURN
  END SUBROUTINE W1EPWQ

!      ****** MAKE TABLE OF F,G,H ******

  SUBROUTINE W1QTBLX
    USE w1comm
    IMPLICIT NONE
    INTEGER,SAVE:: NDMAX_save=0
    INTEGER,SAVE:: NCMAX_save=0
    REAL(rkind),SAVE:: XDMAX_save=0.D0
    INTEGER:: NC,NN,ND
    
    WRITE(6,*) 'NDMAX,NCMAX=',NDMAX,NCMAX
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
  END SUBROUTINE W1QTBLX

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
    ELSEIF(V.GE.XM(NXMAX)) THEN
       VX=YX(NXMAX,NS)
       VK=YK(NXMAX)
       CT0=CS0(NXMAX)
       CT1=CS1(NXMAX)
       CT2=CS2(NXMAX)
       CT3=CS3(NXMAX)
    ELSE
       NHF=INT(NXMAX*(V-XM(1))/(XM(NXMAX)-XM(1)))
       IF(NHF.LT.1)     NHF=1
       IF(NHF.GT.NXMAX) NHF=NXMAX
       IF(V.LT.XM(NHF)) THEN
          DO I=NHF-1,1,-1
             IF(XM(I).LT.V) THEN
                NPV=I
                GOTO 100
             ENDIF
          END DO
       ELSE
          DO I=NHF+1,NXMAX
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

!     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******

  SUBROUTINE W1CLPW(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS,I
    REAL(rkind):: FACT,DX,RKZA
    COMPLEX(rkind):: CPANTKX1,CPANTKY1,CPANTKZ1,CPANTKX2,CPANTKY2,CPANTKZ2
    COMPLEX(rkind):: CFJX1,CFJX2

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

    CPANTKX1=0.D0
    CPANTKY1=0.D0
    CPANTKZ1=0.D0
    CPANTKX2=0.D0
    CPANTKY2=0.D0
    CPANTKZ2=0.D0

    RKZA=RKZ
    IF(ABS(RKZA).LE.1.D-6) RKZA=1.D-6

!   antenna power in HFS

    IF(NXANT1.GE.1) THEN
       CFJX1=CI*RKZA*CFJZ1
       DO NX=1,NXANT1-1
          DX=XA(NX+1)-XA(NX)
          CPANTKX1=CPANTKX1+RZ*DCONJG(CE2DA(NZ,NX,1))*CFJX1*DX
       END DO
       DX=-RD-XA(NXANT1)
       CPANTKX1=CPANTKX1+RZ*DCONJG(CE2DA(NZ,NXANT1,1))*CFJX1*DX
       FACT=(-RD-XA(NXANT1))/(XA(NXANT1+1)-XA(NXANT1))
       CPANTKY1=RZ*(1.D0-FACT)*DCONJG(CE2DA(NZ,NXANT1  ,2))*CFJY1 &
               +RZ*      FACT *DCONJG(CE2DA(NZ,NXANT1+1,2))*CFJY1
       CPANTKZ1=RZ*(1.D0-FACT)*DCONJG(CE2DA(NZ,NXANT1  ,3))*CFJZ1 &
               +RZ*      FACT *DCONJG(CE2DA(NZ,NXANT1+1,3))*CFJZ1
    END IF

!   antenna power in LFS

    IF(NXANT1.GE.2) THEN
       CFJX2=CI*RKZA*CFJZ2
       DO NX=NXANT2+1,NXMAX-1
          DX=XA(NX+1)-XA(NX)
          CPANTKX2=CPANTKX2+RZ*DCONJG(CE2DA(NZ,NX,1))*CFJX2*DX
       END DO
       DX=XA(NXANT2+1)-RD
       CPANTKX2=CPANTKX2+RZ*DCONJG(CE2DA(NZ,NXANT2,1))*CFJX2*DX
       FACT=( RD-XA(NXANT2))/(XA(NXANT2+1)-XA(NXANT2))
       CPANTKY2=RZ*(1.D0-FACT)*DCONJG(CE2DA(NZ,NXANT2  ,2))*CFJY2 &
               +RZ*      FACT *DCONJG(CE2DA(NZ,NXANT2+1,2))*CFJY2
       CPANTKZ2=RZ*(1.D0-FACT)*DCONJG(CE2DA(NZ,NXANT2  ,3))*CFJZ2 &
               +RZ*      FACT *DCONJG(CE2DA(NZ,NXANT2+1,3))*CFJZ2
    END IF
    
    WRITE(6,'(A,I5,1PE12.4)') 'NZ,RKZ,CPANTKX/Y/Z=',NZ,RKZ
    WRITE(6,'(A,1P6E12.4)') '  1:',CPANTKX1,CPANTKY1,CPANTKZ1
    WRITE(6,'(A,1P6E12.4)') '  2:',CPANTKX2,CPANTKY2,CPANTKZ2
    CPANTK1(NZ)=CPANTKX1+CPANTKY1+CPANTKZ1
    CPANTK2(NZ)=CPANTKX2+CPANTKY2+CPANTKZ2

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

END MODULE w1exec7
