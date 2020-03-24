! w1exec11.f90

MODULE w1exec11

  USE w1comm,ONLY: rkind
  USE w1fflr,ONLY: w1fnmn
  USE w1exec8,ONLY: w1qtblx,w1qcal,w1qlni
  USE w1qtbl

  PRIVATE
  PUBLIC w1_exec11

CONTAINS

  SUBROUTINE w1_exec11(NZ,IERR)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    CALL W1DSPQ
    CALL W1BNDQ(NZ,IERR)
       IF(IERR.NE.0) RETURN
    CALL W1EPWQ(NZ)
    CALL W1CLCD(NZ)
    CALL W1CLPW(NZ)
    RETURN
  END SUBROUTINE w1_exec11

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
    INTEGER:: I,J,L,N,N1,N2,M,NX,NX1,NS,IA,IB,NCL,NXD,NQ
    REAL(rkind):: RW,RKV,DTT0,DSS0,DTS1,DTS2,DTTW,DX,RCE
    REAL(rkind):: RKPR,RNPR,FACT,RKZA
    COMPLEX(rkind):: CKKV,CKKV2
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

!    CF(MCEN-5,MLEN-2)=CF(MCEN-5,MLEN-2)+CI*RNPR
    CF(MCEN-2,MLEN-2)=CF(MCEN-2,MLEN-2)+CI*RNPR
    CF(MCEN+2,MLEN-2)=CKKV
    CA(MLEN-2)=CKKV*CFWG1

    CF(MCEN-2,MLEN-1)=1.D0
    CF(MCEN,MLEN-1)=-1.D0
    CA(MLEN-1)=CFWG2

    CF(MCEN-2,MLEN)=1.D0
    CF(MCEN,MLEN)=-1.D0
    CA(MLEN)=CFWG1

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
    INTEGER:: NS,NX,NX1,NXD0,NCL0,IL,NN,MM,IA,IB,NCL1,NXD1
    REAL(rkind):: RW,RKV,RCE,DX,PABSL,RNZ
    REAL(rkind):: PIN1,PIN2,PIN3,PIN4,PIN,POUT1,POUT2,POUT3,POUT4,POUT,PCONV
    COMPLEX(rkind):: CPABSL,CDEY,CDEZ,CBY,CBZ,CLH0,CLH1,CAJ0L,CAJ1L
    COMPLEX(rkind):: CPABS0,CPABS1
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
         'HFS-X:',CFWG4,CA(1),ABS(CA(1))**2, &
         'HFS-O:',CFWG3,CA(2),ABS(CA(2))**2, &
         'LFS-X:',CFWG2,CA(MLEN-1),ABS(CA(MLEN-1))**2, &
         'LFS-O:',CFWG1,CA(MLEN),ABS(CA(MLEN))**2
    RNZ=AKZ(NZ)*VC/RW
    IF(ABS(RNZ).LT.1.D0) THEN
       PIN1=ABS(CFWG1)**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       PIN2=ABS(CFWG2)**2              *SQRT(1.D0-RNZ**2)
       PIN3=ABS(CFWG3)**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       PIN4=ABS(CFWG4)**2              *SQRT(1.D0-RNZ**2)
       PIN=PIN1+PIN2+PIN3+PIN4
       POUT1=ABS(CA(MLEN  ))**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       POUT2=ABS(CA(MLEN-1))**2              *SQRT(1.D0-RNZ**2)
       POUT3=ABS(CA(2     ))**2*(1.D0+RNZ**2)*SQRT(1.D0-RNZ**2)
       POUT4=ABS(CA(1     ))**2              *SQRT(1.D0-RNZ**2)
       POUT=POUT1+POUT2+POUT3+POUT4
       PCONV=PIN-POUT
    ELSE
       PIN=1.D0
       POUT=0.D0
       PCONV=PIN-POUT
    END IF
    WRITE(6,'(A,F7.2,F7.3,1P5E12.4)') &
         'REFL:',AKZ(NZ),RNZ, &
                 POUT1/PIN,POUT2/PIN,POUT3/PIN,POUT4/PIN,PCONV/PIN

    DO NS=1,NSMAX
       DO NX=1,NXMAX
          NN=3*NX-1
          DO IA=1,3
             CAJ0(NN+IA)=0.D0
             CAJ1(NN+IA)=0.D0
          END DO
       END DO
          
       DO NX=1,NXMAX-1
          DX=RKV*(XA(NX+1)-XA(NX))
          DO NX1=MAX(1,NX+NXDMIN),MIN(NX+NXDMAX,NXMAX-1)
             NXD0=NX1-NX+NXDMAX+1 ! positive
             NCL0=NCLA(NXD0,NX,NS)
             NXD1=NX-NX1+NXDMAX+1 ! positive
             NCL1=NCLA(NXD1,NX1,NS)
             CPABSL=0.D0
             IF(NCL0.NE.0.AND.NCL1.NE.0) THEN
                DO IL=1,4
                   SELECT CASE(IL)
                   CASE(1)
                      NN=3*NX-1
                      MM=3*NX1-1
                   CASE(2)
                      NN=3*NX+2
                      MM=3*NX1-1
                   CASE(3)
                      NN=3*NX-1
                      MM=3*NX1+2
                   CASE(4)
                      NN=3*NX+2
                      MM=3*NX1+2
                   END SELECT
                   CPABS0=0.D0
                   CPABS1=0.D0
                   DO IA=1,3
                      CAJ0L=0.D0
                      CAJ1L=0.D0
                      DO IB=1,3
!                         CLH0=0.5D0*(CL(IA,IB,IL,NCL0) &
!                                    -CONJG(CL(IB,IA,IL,NCL0)))
!                         CLH1=0.5D0*(CL(IA,IB,IL,NCL1) &
!                                    -CONJG(CL(IB,IA,IL,NCL1)))
                         CLH0=0.5D0*(CL(IA,IB,IL,NCL0) &
                                    -CONJG(CL(IB,IA,IL,NCL0)))
                         CLH1=0.5D0*(CL(IA,IB,IL,NCL1) &
                                    -CONJG(CL(IB,IA,IL,NCL1)))
                         CAJ0L=CAJ0L+CLH0*CA(MM+IB)
                         CAJ1L=CAJ1L+CLH1*CA(NN+IB)
                      END DO
                      CPABS0=CPABS0+0.5D0*CONJG(CA(NN+IA))*CAJ0L
                      CPABS1=CPABS1+0.5D0*CONJG(CA(MM+IA))*CAJ1L
!                      CABSL=CABSL &
!                           +0.5D0*CONJG(CA(NN+IA))*CAJ0L &
!                           +0.5D0*CONJG(CA(MM+IA))*CAJ1L
!                      CAJ0(NN+IA)=CAJ0(NN+IA)+CAJ0L
!                      CAJ1(MM+IA)=CAJ1(MM+IA)+CAJ1L
                   END DO
                   IF(NX.EQ.4000.AND.NX1.EQ.4001) THEN
                      WRITE(6,'(A,3I5,1P4E12.4)') &
                           'CPABS0:',NX,NX1,IL,CPABS0,CPABS1
                   END IF
                   IF(NX.EQ.4001.AND.NX1.EQ.4000) THEN
                      WRITE(6,'(A,3I5,1P4E12.4)') &
                           'CPABS1:',NX,NX1,IL,CPABS0,CPABS1
                   END IF
                   CPABSL=CPABSL+CPABS0+CPABS1
                END DO
             END IF
             PABSL=-CI*RCE*CPABSL*RKV
             PABS(NX,   NS)=PABS(NX,   NS)+0.5D0*PABSL
             PABS(NX1,  NS)=PABS(NX1,  NS)+0.5D0*PABSL
          END DO
          IF(ABS(PABSL).GT.1.D-8) WRITE(6,'(A,I6,1PE12.4)') 'NX,PABSL=',NX,PABSL
       END DO
       DO NX=1,NXMAX
          CJ2DA(NZ,NX,1,NS)=CAJ0(3*NX)
          CJ2DA(NZ,NX,2,NS)=CAJ0(3*NX+1)
          CJ2DA(NZ,NX,3,NS)=CAJ0(3*NX+2)
          CJ2DB(NZ,NX,1,NS)=CAJ1(3*NX)
          CJ2DB(NZ,NX,2,NS)=CAJ1(3*NX+1)
          CJ2DB(NZ,NX,3,NS)=CAJ1(3*NX+2)
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

!     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******

  SUBROUTINE W1CLPW(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS,I,NZ1,NZ2,NKZ1,NKZ2,NKZ,NZZ,IA,IB
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

  FUNCTION nkz_nz(nz)
    USE w1comm,ONLY: nzmax
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nz
    INTEGER:: nkz_nz

    IF(nz.GT.nzmax/2) THEN
       nkz_nz=nz-nzmax
    ELSE
       nkz_nz=nz-1
    END IF
    RETURN
  END FUNCTION nkz_nz

  FUNCTION nzz_nkz(nkz)
    USE w1comm,ONLY: nzmax
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nkz
    INTEGER:: nzz_nkz

    IF(nkz.GE.0) THEN
       nzz_nkz=nkz+1           ! nkz=0 nzz=1, nkz=nzmax-1 nzz=nzmax
    ELSE
       nzz_nkz=nkz+2*nzmax+1   ! nkz=-nzmax nzz=nzmax+1, nkz=-1 nzz=2*nzmax
    END IF
    RETURN
  END FUNCTION nzz_nkz

  FUNCTION nkz_nzz(nzz)
    USE w1comm,ONLY: nzmax
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nzz
    INTEGER:: nkz_nzz

    IF(nzz.LE.nzmax) THEN
       nkz_nzz=nzz-1           ! nkz=0 nzz=1, nkz=nzmax-1 nzz=nzmax
    ELSE
       nkz_nzz=nzz-2*nzmax-1   ! nkz=-nzmax nzz=nzmax+1, nkz=-1 nzz=2*nzmax
    END IF
    RETURN
  END FUNCTION nkz_nzz

END MODULE w1exec11
