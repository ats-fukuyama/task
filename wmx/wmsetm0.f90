! wmsetm0.f90

MODULE wmsetm0

  USE wmcomm,ONLY: rkind
  COMPLEX(rkind),ALLOCATABLE:: CEMP(:,:,:,:),CFVP(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGM12(:,:),CGM13(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGMH12(:,:),CGMH13(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGMH22(:,:),CGMH23(:,:),CGMH33(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGC11(:,:),CGC12(:,:),CGC13(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGPH22(:,:),CGPH23(:,:),CGPH33(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGP12(:,:),CGP13(:,:)

  PRIVATE
  PUBLIC wm_setm0,wm_setv0

CONTAINS

!     ****** CALCULATE LOCAL COEFFICIENT MATRIX ******

  SUBROUTINE wm_setm0(CA,CB,IGD,LA,NR_prev)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: IGD,LA
    INTEGER,INTENT(INOUT):: NR_prev
    COMPLEX(rkind),INTENT(OUT):: CA(LA),CB
    
    INTEGER:: ICOMP,IGD1,NR,MLXD,IGD2,NKXD,MB

    ! --- igd is the line number of matrix equation ---

    !     ICOMP=1 : R COMPONENT OF MAXWELL EQUATION
    !     ICOMP=2 : THETA COMPONENT OF MAXWELL EQUATION
    !     ICOMP=3 : PHI COMPONENT OF MAXWELL EQUATION

      ICOMP=MOD(IGD-1,3)+1

      IGD1=(IGD-ICOMP)/3+1
      MLXD=MOD(IGD1-1,MDSIZ)+1

      IGD2=(IGD1-MLXD)/MDSIZ+1
      NKXD=MOD(IGD2-1,NDSIZ)+1

      NR=(IGD2-NKXD)/NDSIZ+1

      ! --- if nr is changed, calculate corresponding part of
      !     matrix, right-had-side vector, and boudary condition

      IF(NR.NE.NR_prev) THEN
         NR_prev=NR
         IF(NR.EQ.nr_start) THEN
            IF(ALLOCATED(CEMP)) CALL wm_setm_deallocate
            CALL wm_setm_allocate
            CALL wm_setm_matrix(NR,1)   ! with initial setup
         ELSE
            CALL wm_setm_matrix(NR,0)
         END IF
         CALL wm_setm_vector(NR)
         CALL wm_setm_boundary(NR)
      END IF

      ! --- the data of the line in matrix and RHS vector is substituted ---

      DO MB=1,MBND
         CA(MB)=CEMP(MB,NKXD,MLXD,ICOMP)
      ENDDO
      CB=CFVP(NKXD,MLXD,ICOMP)

    RETURN
  END SUBROUTINE wm_setm0

!     ****** CALCULATE LOCAL right_hand_side vector ******

  SUBROUTINE wm_setv0(NR,CFVP_)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    COMPLEX(rkind),INTENT(OUT):: CFVP_(nhhmax,nthmax,3)
    INTEGER:: ndx,mdx,i
    
    CALL wm_setm_vector(NR)

    DO i=1,3
       DO mdx=1,mdsiz
          DO ndx=1,ndsiz
             CFVP_(ndx,mdx,i)=CFVP(ndx,mdx,i)
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE wm_setv0

  ! --- allocate wmsetm local data ---

  SUBROUTINE wm_setm_allocate

    USE wmcomm
    IMPLICIT NONE
    
    ALLOCATE(CEMP(mbnd,nhhmax,nthmax,3),CFVP(nhhmax,nthmax,3))
    ALLOCATE(CGM12(nthmax_f,nthmax_f),CGM13(nthmax_f,nhhmax_f))
    ALLOCATE(CGMH12(nthmax_f,nhhmax_f),CGMH13(nthmax_f,nhhmax_f))
    ALLOCATE(CGMH22(nthmax_f,nhhmax_f),CGMH23(nthmax_f,nhhmax_f))
    ALLOCATE(CGMH33(nthmax_f,nhhmax_f))
    ALLOCATE(CGC11(nthmax_f,nhhmax_f),CGC12(nthmax_f,nhhmax_f))
    ALLOCATE(CGC13(nthmax_f,nhhmax_f))
    ALLOCATE(CGPH22(nthmax_f,nhhmax_f),CGPH23(nthmax_f,nhhmax_f))
    ALLOCATE(CGPH33(nthmax_f,nhhmax_f))
    ALLOCATE(CGP12(nthmax_f,nhhmax_f),CGP13(nthmax_f,nhhmax_f))

    RETURN
  END SUBROUTINE wm_setm_allocate

  ! --- deallocate wmsetm local data ---

  SUBROUTINE wm_setm_deallocate

    USE wmcomm
    IMPLICIT NONE
    
    DEALLOCATE(CEMP,CFVP)
    DEALLOCATE(CGM12,CGM13,CGMH12,CGMH13)
    DEALLOCATE(CGMH22,CGMH23,CGMH33)
    DEALLOCATE(CGC11,CGC12)
    DEALLOCATE(CGC13,CGPH22,CGPH23)
    DEALLOCATE(CGPH33,CGP12,CGP13)

    RETURN
  END SUBROUTINE wm_setm_deallocate

!     ****** ASSEMBLE TOTAL ELEMENT COEFFICIENT MATRIX ******

  SUBROUTINE wm_setm_matrix(NR,IND)

    USE wmcomm
    USE wmsetf0
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,IND
      COMPLEX(rkind):: CDVM(3,3),CDVC(3,3),CDVP(3,3)
      INTEGER:: KDX,LDX,J,I,KDXF,LDXF,LBAND,NKX,MLX
      INTEGER:: ND,NDX,NKD,NN,NK,KKD,KKDX,KD,MD,MDX,MLD,MM,ML,LD
      INTEGER:: ID,LBND
      REAL(rkind):: XRHOM,XRHOC,XRHOP,XRHOMH,XRHOPH,DRHOM,DRHOP,DRHOPM
      REAL(rkind):: QPMH,QPC,QPPH,DPSIPDRHOMH,DPSIPDRHOPH,DPSIPDRHOC
      REAL(rkind):: FMHM,FMHC,FPHC,FPHP,DFMHM,DFMHC
      REAL(rkind):: DFCM,DFCC,DFCP,DDFMHM,DDFMHC,DDFPHC,DDFPHP
      REAL(rkind):: FCMH,FCPH,DFCMH,DFCPH
      REAL(rkind):: FACT1M,FACT1C,FACT1P,FACT2M,FACT2C,FACT2P
      REAL(rkind):: FACT3M,FACT3C,FACT3P
      
      ! initial setting to set up 2) data for NR,  3) data for NR+1

      IF(IND.EQ.1) THEN

         CALL wm_setf0(NR,0) !  wm_setf0 setup 3) data for NR
                             !     CGD: metric tensor
                             !     CMA: conversion tensor
                             !     CGD: dielectric tensor 

         ! shift 3) data to 2) data for NR
         
         DO KDX=1,KDSIZ_F
            DO LDX=1,LDSIZ_F
               CGF11(LDX,KDX,2)=CGF11(LDX,KDX,3)
               CGF12(LDX,KDX,2)=CGF12(LDX,KDX,3)
               CGF13(LDX,KDX,2)=CGF13(LDX,KDX,3)
               CGF22(LDX,KDX,2)=CGF22(LDX,KDX,3)
               CGF23(LDX,KDX,2)=CGF23(LDX,KDX,3)
               CGF33(LDX,KDX,2)=CGF33(LDX,KDX,3)
            ENDDO
         ENDDO

         DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ_F
               DO MDX=1,MDSIZ
                  DO LDX=1,LDSIZ_F
                     DO J=1,3
                        DO I=1,3
                           CGD(I,J,LDX,MDX,KDX,NDX,2) &
                                =CGD(I,J,LDX,MDX,KDX,NDX,3)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         CALL wm_setf0(NR+1,0) !  wm_setf setup 3) data for NR+1

      ENDIF

      ! --- now shift 2) data for NR to 1) data,
      ! --- and shift 3) data for NR+1 to 2) daata

      DO KDX=1,KDSIZ_F
         DO LDX=1,LDSIZ_F
            CGF11(LDX,KDX,1)=CGF11(LDX,KDX,2)
            CGF12(LDX,KDX,1)=CGF12(LDX,KDX,2)
            CGF13(LDX,KDX,1)=CGF13(LDX,KDX,2)
            CGF22(LDX,KDX,1)=CGF22(LDX,KDX,2)
            CGF23(LDX,KDX,1)=CGF23(LDX,KDX,2)
            CGF33(LDX,KDX,1)=CGF33(LDX,KDX,2)
            CGF11(LDX,KDX,2)=CGF11(LDX,KDX,3)
            CGF12(LDX,KDX,2)=CGF12(LDX,KDX,3)
            CGF13(LDX,KDX,2)=CGF13(LDX,KDX,3)
            CGF22(LDX,KDX,2)=CGF22(LDX,KDX,3)
            CGF23(LDX,KDX,2)=CGF23(LDX,KDX,3)
            CGF33(LDX,KDX,2)=CGF33(LDX,KDX,3)
         ENDDO
      ENDDO

      DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ_F
            DO MDX=1,MDSIZ
               DO LDX=1,LDSIZ_F
                  DO J=1,3
                     DO I=1,3
                        CGD(I,J,LDX,MDX,KDX,NDX,1)=CGD(I,J,LDX,MDX,KDX,NDX,2)
                        CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF(NR.LT.NRMAX) CALL wm_setf0(NR+2,0) !  wm_setf setup 3) data for NR+2

      IF(NR.EQ.1) THEN
         XRHOM = XRHO(2)/1.D6
      ELSE
         XRHOM = XRHO(NR)
      ENDIF
      XRHOC = XRHO(NR+1)
      IF(NR.EQ.NRMAX) THEN
         XRHOP=XRHO(NR+1)
      ELSE
         XRHOP=XRHO(NR+2)
      ENDIF
      XRHOMH=0.5D0*(XRHOM+XRHOC)
      XRHOPH=0.5D0*(XRHOC+XRHOP)

      DRHOM =XRHOC-XRHOM
      IF(NR.EQ.NRMAX) THEN
         DRHOP =XRHOC-XRHOM
      ELSE
         DRHOP =XRHOP-XRHOC
      ENDIF
      DRHOPM=DRHOM+DRHOP

      IF(MODELG.EQ.3) THEN
         QPMH=0.5D0*(QPS(NR)+QPS(NR+1))
         QPC=QPS(NR+1)
         IF(NR.EQ.NRMAX) THEN
            QPPH=QPS(NR+1)
         ELSE
            QPPH=0.5D0*(QPS(NR+1)+QPS(NR+2))
         ENDIF
         DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
         DPSIPDRHOPH=2.D0*PSITA*XRHOPH/QPPH
         DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
      ELSE
         DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
         DPSIPDRHOPH=2.D0*PSIPA*XRHOPH
         DPSIPDRHOC =2.D0*PSIPA*XRHOC
      ENDIF

      FMHM=0.5D0
      FMHC=0.5D0
      FPHC=0.5D0
      FPHP=0.5D0
      DFMHM=-1.0D0/DRHOM   /DPSIPDRHOMH
      DFMHC= 1.0D0/DRHOM   /DPSIPDRHOMH

      DFCM= -(DRHOP*DPSIPDRHOPH)/(DRHOM*DRHOPM) &
             /(DPSIPDRHOMH*DPSIPDRHOC)
      DFCC=  (DRHOP*DPSIPDRHOPH)-(DRHOM*DPSIPDRHOMH)/(DRHOP*DRHOM) &
             /(DPSIPDRHOPH*DPSIPDRHOMH)
      DFCP=  (DRHOM*DPSIPDRHOMH)/(DRHOP*DRHOPM) &
             /(DPSIPDRHOPH*DPSIPDRHOC)   !
      DDFMHM= 2.D0/(DRHOM*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOMH)
      DDFMHC=-2.D0/(DRHOM*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOMH)
      DDFPHC=-2.D0/(DRHOP*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOPH)
      DDFPHP= 2.D0/(DRHOP*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOPH)

      FCMH  = (DRHOP*DPSIPDRHOPH)/(DRHOPM*DPSIPDRHOC)
      FCPH  = (DRHOM*DPSIPDRHOMH)/(DRHOPM*DPSIPDRHOC)
      DFCMH =-2.D0 /DRHOPM /DPSIPDRHOC
      DFCPH = 2.D0 /DRHOPM /DPSIPDRHOC

      FACT1M=XRHOMH/XRHOM
      FACT1C=XRHOMH/XRHOC
      FACT1P=XRHOPH/XRHOC

      FACT2M=XRHOMH/XRHOM
      FACT2C=XRHOMH/XRHOC
      FACT2P=XRHOPH/XRHOC

      FACT3M=XRHOMH/XRHOM
      FACT3C=XRHOMH/XRHOC
      FACT3P=XRHOPH/XRHOC

      DO KDX=1,KDSIZ_F
         DO LDX=1,LDSIZ_F

            LDXF=LDX
            KDXF=KDX

            CGM12(LDX,KDX)=       CGF12(LDXF,KDXF,1)
            CGM13(LDX,KDX)=       CGF13(LDXF,KDXF,1)/XRHOM

            CGMH12(LDX,KDX)= FMHM*CGF12(LDXF,KDXF,1) &
                            +FMHC*CGF12(LDXF,KDXF,2)
            CGMH13(LDX,KDX)=(FMHM*CGF13(LDXF,KDXF,1) &
                            +FMHC*CGF13(LDXF,KDXF,2))/XRHOMH
            CGMH22(LDX,KDX)=(FMHM*CGF22(LDXF,KDXF,1) &
                            +FMHC*CGF22(LDXF,KDXF,2))*XRHOMH**2
            CGMH23(LDX,KDX)=(FMHM*CGF23(LDXF,KDXF,1) &
                            +FMHC*CGF23(LDXF,KDXF,2))*XRHOMH
            CGMH33(LDX,KDX)= FMHM*CGF33(LDXF,KDXF,1) &
                            +FMHC*CGF33(LDXF,KDXF,2)

            CGC11(LDX,KDX)=       CGF11(LDXF,KDXF,2)/XRHOC**2
            CGC12(LDX,KDX)=       CGF12(LDXF,KDXF,2)
            CGC13(LDX,KDX)=       CGF13(LDXF,KDXF,2)/XRHOC

            CGPH22(LDX,KDX)=(FPHC*CGF22(LDXF,KDXF,2) &
                            +FPHP*CGF22(LDXF,KDXF,3))*XRHOPH**2
            CGPH23(LDX,KDX)=(FPHC*CGF23(LDXF,KDXF,2) &
                            +FPHP*CGF23(LDXF,KDXF,3))*XRHOPH
            CGPH33(LDX,KDX)= FPHC*CGF33(LDXF,KDXF,2) &
                            +FPHP*CGF33(LDXF,KDXF,3)

            CGP12(LDX,KDX)=       CGF12(LDXF,KDXF,3)
            CGP13(LDX,KDX)=       CGF13(LDXF,KDXF,3)/XRHOP

         ENDDO
      ENDDO

!        ND : (n - n0) / Np
!        KD : n'              = (k - n) / Np
!        NN : n               = n0 +  ND       * Np
!        NK : k = n + n' * Np = n0 + (ND + KD) * Np

!        MD : m - m0
!        LD : m'              = l - m
!        MM : m               = m0 +  MD
!        ML : l = m + m'      = m0 + (MD + LD)

      DO LBAND=1,MBND
         DO NKX=1,NDSIZ
            DO MLX=1,MDSIZ
               CEMP(LBAND,NKX,MLX,1)=0.D0
               CEMP(LBAND,NKX,MLX,2)=0.D0
               CEMP(LBAND,NKX,MLX,3)=0.D0
            ENDDO
         ENDDO
      ENDDO

      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
         DO NKD=NDMIN,NDMAX
            NKX=NKD-NDMIN+1
            KD=NKD-ND
            IF(MODELK.EQ.0.OR. &
                 (KD.GE.KDMIN_F.AND. &
                 (KD.LE.KDMAX_F))) THEN
               KDX=KD-KDMIN_F + 1
               KDXF=KDX

               NN=NPH0+NHC*ND
               NK=NPH0+NHC*NKD

               DO MD=MDMIN,MDMAX
                  MDX=MD-MDMIN+1
                  DO MLD=MDMIN,MDMAX
                     MLX=MLD-MDMIN+1
                     LD=MLD-MD
                     IF(LD.GE.LDMIN_F.AND. &
                       (LD.LE.LDMAX_F.OR.LDMAX_F==0))THEN
                        LDX =LD-LDMIN_F+1
                        MM=NTH0+MD
                        ML=NTH0+MLD
                        DO J=1,3
                           DO I=1,3
                              CDVM(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,1)
                              CDVC(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,2)
                              CDVP(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,3)
                           ENDDO
                        ENDDO

                        ID=3*MDSIZ*NDSIZ

!     --- R COMPONENT OF MAXWELL EQUATION ---

                        LBND=MCENT-3*KD*MDSIZ-3*LD-1

                        CEMP(LBND+1    ,NKX,MLX,1) &
                             =CEMP(LBND+1   ,NKX,MLX,1) &
                             +( &
                             -ML*NN*CGMH23(LDX,KDX) &
                             +ML*MM*CGMH33(LDX,KDX) &
                             +NK*NN*CGMH22(LDX,KDX) &
                             -NK*MM*CGMH23(LDX,KDX) &
                             +0.5D0*CDVM(1,1)*FACT1M &
                             +0.5D0*CDVC(1,1)*FACT1C &
                             )/XRHOMH

                        CEMP(LBND+2-ID,NKX,MLX,1) &
                             =CEMP(LBND+2-ID,NKX,MLX,1) &
                             +( &
                             +ML*NN*CGMH13(LDX,KDX)*FMHM &
                             -NK*NN*CGMH12(LDX,KDX)*FMHM &
                             +CI*ML*CGMH33(LDX,KDX)*DFMHM &
                             -CI*NK*CGMH23(LDX,KDX)*DFMHM &
                             +CDVM(1,2)*FMHM &
                             )*XRHOM

                        CEMP(LBND+2   ,NKX,MLX,1) &
                             =CEMP(LBND+2   ,NKX,MLX,1) &
                             +( &
                             +ML*NN*CGMH13(LDX,KDX)*FMHC &
                             -NK*NN*CGMH12(LDX,KDX)*FMHC &
                             +CI*ML*CGMH33(LDX,KDX)*DFMHC &
                             -CI*NK*CGMH23(LDX,KDX)*DFMHC &
                             +CDVC(1,2)*FMHC &
                             )*XRHOC

                        CEMP(LBND+3-ID,NKX,MLX,1) &
                             =CEMP(LBND+3-ID,NKX,MLX,1) &
                             +( &
                             -ML*MM*CGMH13(LDX,KDX)*FMHM &
                             +NK*MM*CGMH12(LDX,KDX)*FMHM &
                             -CI*ML*CGMH23(LDX,KDX)*DFMHM &
                             +CI*NK*CGMH22(LDX,KDX)*DFMHM &
                             +CDVM(1,3)*FMHM &
                             )

                        CEMP(LBND+3   ,NKX,MLX,1) &
                             =CEMP(LBND+3   ,NKX,MLX,1) &
                             +( &
                             -ML*MM*CGMH13(LDX,KDX)*FMHC &
                             +NK*MM*CGMH12(LDX,KDX)*FMHC &
                             -CI*ML*CGMH23(LDX,KDX)*DFMHC &
                             +CI*NK*CGMH22(LDX,KDX)*DFMHC &
                             +CDVC(1,3)*FMHC &
                             )
!     --- THETA COMPONENT OF MAXWELL EQUATION ---

                        LBND=MCENT-3*KD*MDSIZ-3*LD-2

                        CEMP(LBND+1   ,NKX,MLX,2) &
                             =CEMP(LBND+1   ,NKX,MLX,2) &
                             +( &
                             -NK*NN*CGC12(LDX,KDX)*FCMH &
                             +NK*MM*CGC13(LDX,KDX)*FCMH &
                             -CI*NN*CGMH23(LDX,KDX)*DFCMH &
                             +CI*MM*CGMH33(LDX,KDX)*DFCMH &
                             +CDVC(2,1)*FCMH*FACT2C &
                             )/XRHOMH

                        CEMP(LBND+1+ID,NKX,MLX,2) &
                             =CEMP(LBND+1+ID,NKX,MLX,2) &
                             +( &
                             -NK*NN*CGC12(LDX,KDX)*FCPH &
                             +NK*MM*CGC13(LDX,KDX)*FCPH &
                             -CI*NN*CGPH23(LDX,KDX)*DFCPH &
                             +CI*MM*CGPH33(LDX,KDX)*DFCPH &
                             +CDVC(2,1)*FCPH*FACT2P &
                             )/XRHOPH

                        CEMP(LBND+2-ID,NKX,MLX,2) &
                             =CEMP(LBND+2-ID,NKX,MLX,2) &
                             +( &
                             +CI*NK*CGC13(LDX,KDX)*DFCM &
                             +CI*NN*CGM13(LDX,KDX)*DFCM &
                             -      CGMH33(LDX,KDX)*DDFMHM &
                             )*XRHOM

                        CEMP(LBND+2   ,NKX,MLX,2) &
                             =CEMP(LBND+2   ,NKX,MLX,2) &
                             +( &
                             +NK*NN*CGC11(LDX,KDX) &
                             +CI*NK*CGC13(LDX,KDX)*DFCC &
                             +CI*NN*CGC13(LDX,KDX)*DFCC &
                             -      CGMH33(LDX,KDX)*DDFMHC &
                             -      CGPH33(LDX,KDX)*DDFPHC &
                             +CDVC(2,2) &
                             )*XRHOC

                        CEMP(LBND+2+ID,NKX,MLX,2) &
                             =CEMP(LBND+2+ID,NKX,MLX,2) &
                             +( &
                             +CI*NK*CGC13(LDX,KDX)*DFCP &
                             +CI*NN*CGP13(LDX,KDX)*DFCP &
                             -      CGPH33(LDX,KDX)*DDFPHP &
                             )*XRHOP

                        CEMP(LBND+3-ID,NKX,MLX,2) &
                             =CEMP(LBND+3-ID,NKX,MLX,2) &
                             +( &
                             -CI*NK*CGC12(LDX,KDX)*DFCM &
                             -CI*MM*CGM13(LDX,KDX)*DFCM &
                             +      CGMH23(LDX,KDX)*DDFMHM &
                             )

                        CEMP(LBND+3   ,NKX,MLX,2) &
                             =CEMP(LBND+3   ,NKX,MLX,2) &
                             +( &
                             -NK*MM*CGC11(LDX,KDX) &
                             -CI*NK*CGC12(LDX,KDX)*DFCC &
                             -CI*MM*CGC13(LDX,KDX)*DFCC &
                             +      CGMH23(LDX,KDX)*DDFMHC &
                             +      CGPH23(LDX,KDX)*DDFPHC &
                             +CDVC(2,3) &
                             )

                        CEMP(LBND+3+ID,NKX,MLX,2) &
                             =CEMP(LBND+3+ID,NKX,MLX,2) &
                             +( &
                             -CI*NK*CGC12(LDX,KDX)*DFCP &
                             -CI*MM*CGP13(LDX,KDX)*DFCP &
                             +      CGPH23(LDX,KDX)*DDFPHP &
                             )

!     --- PHI COMPONENT OF MAXWELL EQUATION ---

                        LBND=MCENT-3*KD*MDSIZ-3*LD-3

                        CEMP(LBND+1   ,NKX,MLX,3) &
                             =CEMP(LBND+1   ,NKX,MLX,3) &
                             +( &
                             +ML*NN*CGC12(LDX,KDX)*FCMH &
                             -ML*MM*CGC13(LDX,KDX)*FCMH &
                             +CI*NN*CGMH22(LDX,KDX)*DFCMH &
                             -CI*MM*CGMH23(LDX,KDX)*DFCMH &
                             +CDVC(3,1)*FCMH*FACT3C &
                             )/XRHOMH

                        CEMP(LBND+1+ID,NKX,MLX,3) &
                             =CEMP(LBND+1+ID,NKX,MLX,3) &
                             +( &
                             +ML*NN*CGC12(LDX,KDX)*FCPH &
                             -ML*MM*CGC13(LDX,KDX)*FCPH &
                             +CI*NN*CGPH22(LDX,KDX)*DFCPH &
                             -CI*MM*CGPH23(LDX,KDX)*DFCPH &
                             +CDVC(3,1)*FCPH*FACT3P &
                             )/XRHOPH

                        CEMP(LBND+2-ID,NKX,MLX,3) &
                             =CEMP(LBND+2-ID,NKX,MLX,3) &
                             +( &
                             -CI*ML*CGC13(LDX,KDX)*DFCM &
                             -CI*NN*CGM12(LDX,KDX)*DFCM &
                             +      CGMH23(LDX,KDX)*DDFMHM &
                             )*XRHOM

                        CEMP(LBND+2   ,NKX,MLX,3) &
                             =CEMP(LBND+2   ,NKX,MLX,3) &
                             +( &
                             -ML*NN*CGC11(LDX,KDX) &
                             -CI*ML*CGC13(LDX,KDX)*DFCC &
                             -CI*NN*CGC12(LDX,KDX)*DFCC &
                             +      CGMH23(LDX,KDX)*DDFMHC &
                             +      CGPH23(LDX,KDX)*DDFPHC &
                             +CDVC(3,2) &
                             )*XRHOC

                        CEMP(LBND+2+ID,NKX,MLX,3) &
                             =CEMP(LBND+2+ID,NKX,MLX,3) &
                             +( &
                             -CI*ML*CGC13(LDX,KDX)*DFCP &
                             -CI*NN*CGP12(LDX,KDX)*DFCP &
                             +      CGPH23(LDX,KDX)*DDFPHP &
                             )*XRHOP

                        CEMP(LBND+3-ID,NKX,MLX,3) &
                             =CEMP(LBND+3-ID,NKX,MLX,3) &
                             +( &
                             +CI*ML*CGC12(LDX,KDX)*DFCM &
                             +CI*MM*CGM12(LDX,KDX)*DFCM &
                             -      CGMH22(LDX,KDX)*DDFMHM &
                             )

                        CEMP(LBND+3   ,NKX,MLX,3) &
                             =CEMP(LBND+3   ,NKX,MLX,3) &
                             +( &
                             +ML*MM*CGC11(LDX,KDX) &
                             +CI*ML*CGC12(LDX,KDX)*DFCC &
                             +CI*MM*CGC12(LDX,KDX)*DFCC &
                             -      CGMH22(LDX,KDX)*DDFMHC &
                             -      CGPH22(LDX,KDX)*DDFPHC &
                             +CDVC(3,3) &
                             )

                        CEMP(LBND+3+ID,NKX,MLX,3) &
                             =CEMP(LBND+3+ID,NKX,MLX,3) &
                             +( &
                             +CI*ML*CGC12(LDX,KDX)*DFCP &
                             +CI*MM*CGP12(LDX,KDX)*DFCP &
                             -      CGPH22(LDX,KDX)*DDFPHP &
                             )
                     END IF
                  ENDDO ! MLD
               ENDDO ! MD
            END IF
         ENDDO ! NKD
      ENDDO ! ND
    END SUBROUTINE wm_setm_matrix

!     ****** ASSEMBLE TOTAL ELEMENT FREE VECTOR ******

    SUBROUTINE wm_setm_vector(NR)

      USE wmcomm
      USE plprof,ONLY: pl_prof2
      IMPLICIT NONE

      INTEGER,INTENT(IN):: NR
      REAL(rkind):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
      INTEGER:: MDX,NDX,ND,NN,MD,MM,NRI,NRANT
      REAL(rkind):: DPH,DTH,XRHO1,XRHO2,DRHO,XRHOC,FACTM,FACTP
      REAL(rkind):: RT,RJFACT,QPC,DPSIPDRHOC
      COMPLEX(rkind):: CW,CC,CJTHM,CJPHM,CJR,CJTHP,CJPHP

      DO MDX=1,MDSIZ
         DO NDX=1,NDSIZ
            CFVP(NDX,MDX,1)=0.D0
            CFVP(NDX,MDX,2)=0.D0
            CFVP(NDX,MDX,3)=0.D0
         ENDDO
      ENDDO

      IF(MODEEG.EQ.0) THEN
         DO NRI=1,NRMAX 
            IF(XR(NRI)/RD.LT.1.D0) NRANT=NRI
         ENDDO

         IF(NR+1.GE.NRANT) THEN
            CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
            CC=CI*CW*RMU0
            DPH=2.D0*PI/NHHMAX
            DTH=2.D0*PI/NTHMAX

            IF(NR+1.EQ.NRANT.OR.NR+1.EQ.NRANT+1) THEN
               XRHO1=XRHO(NRANT)
               XRHO2=XRHO(NRANT+1)
               DRHO=XRHO2-XRHO1
               XRHOC=0.5D0*(XRHO2+XRHO1)

               IF(MODELG.EQ.3) THEN
                  QPC=0.5D0*(QPS(NRANT)+QPS(NRANT+1))
                  DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
               ELSE
                  DPSIPDRHOC =2.D0*PSIPA*XRHOC
               ENDIF

               FACTM=(XRHO2-RD/RA)/DRHO
               FACTP=(RD/RA-XRHO1)/DRHO

               DO ND=NDMIN,NDMAX
                  NDX=ND-NDMIN+1
                  NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*FACTM*XRHOC &
                                              /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*FACTM*XRHOC &
                                              /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM) &
                                              *DPSIPDRHOC*DRHO/XRHOC
                     CJTHP=CC*CJANT(2,MDX,NDX)*FACTP*XRHOC  &
                                              /(DPSIPDRHOC*DRHO*DPH)
                     CJPHP=CC*CJANT(3,MDX,NDX)*FACTP*XRHOC &
                                              /(DPSIPDRHOC*DRHO*DTH)
                     CFVP(NDX,MDX,1)=0.D0
                     CFVP(NDX,MDX,2)=0.D0
                     CFVP(NDX,MDX,3)=0.D0
                     IF(NR+1.EQ.NRANT) THEN
                        CFVP(NDX,MDX,2)=CJTHM
                        CFVP(NDX,MDX,3)=CJPHM
                     ELSE IF(NR+1.EQ.NRANT+1) THEN
                        CFVP(NDX,MDX,1)=CJR
                        CFVP(NDX,MDX,2)=CJTHP
                        CFVP(NDX,MDX,3)=CJPHP
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               XRHO1=XRHO(NR+1)
               IF(NR.LT.NRMAX) THEN
                  XRHO2=XRHO(NR+2)
               ELSE
                  XRHO2=2*XRHO(NR+1)-XRHO(NR)
               ENDIF
               DRHO=XRHO2-XRHO1
               XRHOC=0.5D0*(XRHO2+XRHO1)

               IF(MODELG.EQ.3) THEN
                  IF(NR.LT.NRMAX) THEN
                     QPC=0.5D0*(QPS(NR+1)+QPS(NR+2))
                  ELSE
                     QPC=0.5D0*(3*QPS(NR+1)-QPS(NR))
                  ENDIF
                  DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
               ELSE
                  DPSIPDRHOC =2.D0*PSIPA*XRHOC
               ENDIF

               DO ND=NDMIN,NDMAX
                  NDX=ND-NDMIN+1
                  NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*XRHOC &
                          /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*XRHOC &
                          /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DPSIPDRHOC*DRHO &
                          /XRHOC
                     CFVP(NDX,MDX,1)=CJR
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ELSE
         DPH=2.D0*PI/NHHMAX
         DTH=2.D0*PI/NTHMAX
         XRHO1=XRHO(NR+1)
         IF(NR.LT.NRMAX) THEN
            XRHO2=XRHO(NR+2)
         ELSE
            XRHO2=2*XRHO(NR+1)-XRHO(NR)
         ENDIF
         DRHO=XRHO2-XRHO1
         XRHOC=0.5D0*(XRHO2+XRHO1)

         IF(MODELG.EQ.3) THEN
            IF(NR.LT.NRMAX) THEN
               QPC=0.5D0*(QPS(NR+1)+QPS(NR+2))
            ELSE
               QPC=0.5D0*(3*QPS(NR+1)-QPS(NR))
            ENDIF
            DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
         ELSE
            DPSIPDRHOC =2.D0*PSIPA*XRHOC
         ENDIF

         CALL PL_PROF2(XRHO1,RN,RTPR,RTPP,RU)
         RT=(RTPR(1)+2*RTPP(1))
         RJFACT=RN(1)*RT
         RJFACT=RJFACT*XRHO(NR+1)
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+NHC*ND
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               MM=NTH0+MD
               CJTHM=RJFACT*XRHOC/(DPSIPDRHOC*DRHO*DPH)
               CJPHM=RJFACT*XRHOC/(DPSIPDRHOC*DRHO*DTH)
               CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DPSIPDRHOC*DRHO/XRHOC
               CFVP(NDX,MDX,1)=0.D0
               CFVP(NDX,MDX,2)=CJTHM
               CFVP(NDX,MDX,3)=CJPHM
            ENDDO
         ENDDO
      ENDIF

      RETURN
    END SUBROUTINE wm_setm_vector

!     ****** SET BOUNDARY CONDITION ******

    SUBROUTINE wm_setm_boundary(NR)

      USE wmcomm
      IMPLICIT NONE
      
      INTEGER,INTENT(IN):: NR
      INTEGER:: ID0,ID,ND,NDX,NKD,NKX,KD,KDX,MD,MDX,MLD,MLX,LDX,MM,LBND,MB,LD
      REAL(rkind):: DRHO1,DRHO2,A1,A2

      DRHO1=(XRHO(2)-XRHO(1))**2
      DRHO2=(XRHO(3)-XRHO(1))**2
      A1= DRHO2/(DRHO2-DRHO1)
      A2=-DRHO1/(DRHO2-DRHO1)

!        ****** R=0 ******

      ID0=MDSIZ*NDSIZ

      IF(NR.EQ.1) THEN

         ID=3*MDSIZ*NDSIZ         
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            DO NKD=NDMIN,NDMAX
               NKX=NKD-NDMIN+1
               KD=NKD-ND
               IF(MODELK.EQ.0.OR. &
                    (KD.GE.KDMIN.AND.KD.LE.KDMAX)) THEN
                  KDX=MOD(KD-KDMIN+2*KDSIZ,KDSIZ)+1+KDMIN-KDMIN_F
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     DO MLD=MDMIN,MDMAX
                        MLX=MLD-MDMIN+1
                        LD=MLD-MD
                        IF(MODELK.EQ.0.OR. &
                             (LD.GE.LDMIN.AND.LD.LE.LDMAX)) THEN
                           LDX=MOD(LD-LDMIN+2*LDSIZ,LDSIZ)+1+LDMIN-LDMIN_F
                           MM=NTH0+MD

!        ****** EPH'(0) = 0 FOR MM.EQ.0 ******

                           IF(MM.EQ.0 .and. &
                                CMAF(3,3,-MDMIN+1,-NDMIN+1,1) .NE.0 ) THEN
                      
                              LBND=MCENT-3*KD*MDSIZ-3*LD-1
                              CEMP(LBND   +3,NKX,MLX,1) &
                                   =CEMP(LBND   +3,NKX,MLX,1) &
                                   +CEMP(LBND-ID+3,NKX,MLX,1)*A1

                              CEMP(LBND+ID+3,NKX,MLX,1) &
                                   =CEMP(LBND+ID+3,NKX,MLX,1) &
                                   +CEMP(LBND-ID+3,NKX,MLX,1)*A2

                              CEMP(LBND-ID+3,NKX,MLX,1)=0.D0

                              LBND=MCENT-3*KD*MDSIZ-3*LD-2

                              CEMP(LBND   +3,NKX,MLX,2) &
                                   =CEMP(LBND   +3,NKX,MLX,2) &
                                   +CEMP(LBND-ID+3,NKX,MLX,2)*A1

                              CEMP(LBND+ID+3,NKX,MLX,2) &
                                   =CEMP(LBND+ID+3,NKX,MLX,2) &
                                   +CEMP(LBND-ID+3,NKX,MLX,2)*A2

                              CEMP(LBND-ID+3,NKX,MLX,2)=0.D0

                              LBND=MCENT-3*KD*MDSIZ-3*LD-3

                              CEMP(LBND   +3,NKX,MLX,3) &
                                   =CEMP(LBND   +3,NKX,MLX,3) &
                                   +CEMP(LBND-ID+3,NKX,MLX,3)*A1

                              CEMP(LBND+ID+3,NKX,MLX,3) &
                                   =CEMP(LBND+ID+3,NKX,MLX,3) &
                                   +CEMP(LBND-ID+3,NKX,MLX,3)*A2

                              CEMP(LBND-ID+3,NKX,MLX,3)=0.D0


!        ****** ETH'(0) = 0 FOR ABS(MM).EQ.1 ******

                           ELSEIF(ABS(MM).EQ.1) THEN

                              LBND=MCENT-3*KD*MDSIZ-3*LD-1
                              CEMP(LBND   +2,NKX,MLX,1) &
                                   =CEMP(LBND   +2,NKX,MLX,1) &
                                   +CEMP(LBND-ID+2,NKX,MLX,1)*A1

                              CEMP(LBND+ID+2,NKX,MLX,1) &
                                   =CEMP(LBND+ID+2,NKX,MLX,1) &
                                   +CEMP(LBND-ID+2,NKX,MLX,1)*A2

                              LBND=MCENT-3*KD*MDSIZ-3*LD-2

                              CEMP(LBND   +2,NKX,MLX,2) &
                                   =CEMP(LBND   +2,NKX,MLX,2) &
                                   +CEMP(LBND-ID+2,NKX,MLX,2)*A1
                              CEMP(LBND+ID+2,NKX,MLX,2) &
                                   =CEMP(LBND+ID+2,NKX,MLX,2) &
                                   +CEMP(LBND-ID+2,NKX,MLX,2)*A2
                              LBND=MCENT-3*KD*MDSIZ-3*LD-3

                              CEMP(LBND   +2,NKX,MLX,3) &
                                   =CEMP(LBND   +2,NKX,MLX,3) &
                                   +CEMP(LBND-ID+2,NKX,MLX,3)*A1
                              CEMP(LBND+ID+2,NKX,MLX,3) &
                                   =CEMP(LBND+ID+2,NKX,MLX,3) &
                                   +CEMP(LBND-ID+2,NKX,MLX,3)*A2
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!        ****** ETH = 0, EPH =0 AT R=RA ******

      IF(NR.EQ.NRMAX) THEN

         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               DO MB=1,MBND
                  CEMP(MB,NDX,MDX,1)= 0.D0
                  CEMP(MB,NDX,MDX,2)= 0.D0
                  CEMP(MB,NDX,MDX,3)= 0.D0
               ENDDO
               CEMP(MBND,NDX,MDX,1)= 1.D0
               CEMP(MBND-3*MDSIZ*NDSIZ,NDX,MDX,1)= -1.D0
               CEMP(MCENT,NDX,MDX,2)= 1.D0
               CEMP(MCENT,NDX,MDX,3)= 1.D0
               CFVP(NDX,MDX,1)= 0.D0
               CFVP(NDX,MDX,2)= CEWALL(2,MDX,NDX)
               CFVP(NDX,MDX,3)= CEWALL(3,MDX,NDX)
               WRITE(22,'(A,2I8,4ES12.4)') &
                    'CFVP:',MDX,NDX,CEWALL(2,MDX,NDX),CEWALL(3,MDX,NDX)
            ENDDO
         ENDDO

      ENDIF

!     ****** ELIMINATE EFLD AT MDMAX ******

      IF(MDSIZ.GT.1) THEN
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            MD=MDMAX
            MDX=MD-MDMIN+1
            DO MB=1,MBND
               CEMP(MB,NDX,MDX,1)= 0.D0
               CEMP(MB,NDX,MDX,2)= 0.D0
               CEMP(MB,NDX,MDX,3)= 0.D0
            ENDDO
            CEMP(MCENT,NDX,MDX,1)= 1.D0
            CEMP(MCENT,NDX,MDX,2)= 1.D0
            CEMP(MCENT,NDX,MDX,3)= 1.D0
            CFVP(NDX,MDX,1)= 0.D0
            CFVP(NDX,MDX,2)= 0.D0
            CFVP(NDX,MDX,3)= 0.D0
         ENDDO
      ENDIF

!     ****** ELIMINATE EFLD AT NDMAX ******

      IF(NDSIZ.GT.1) THEN
         ND=NDMAX
         NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            DO MB=1,MBND
               CEMP(MB,NDX,MDX,1)= 0.D0
               CEMP(MB,NDX,MDX,2)= 0.D0
               CEMP(MB,NDX,MDX,3)= 0.D0
            ENDDO
            CEMP(MCENT,NDX,MDX,1)= 1.D0
            CEMP(MCENT,NDX,MDX,2)= 1.D0
            CEMP(MCENT,NDX,MDX,3)= 1.D0
            CFVP(NDX,MDX,1)= 0.D0
            CFVP(NDX,MDX,2)= 0.D0
            CFVP(NDX,MDX,3)= 0.D0
         ENDDO
      ENDIF

      RETURN
    END SUBROUTINE wm_setm_boundary
  END MODULE wmsetm0
