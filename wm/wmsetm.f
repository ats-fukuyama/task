C     $Id$
C
C     ****** CALCULATE LOCAL COEFFICIENT MATRIX ******
C
      SUBROUTINE WMSETM(CA,CB,IGD,LA,NRP)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CA(LA)
      COMMON /WMCMV1/ CEMP(MBNDM,NDM,MDM,3),CFVP(NDM,MDM,3)
C
C     ICOMP=1 : R COMPONENT OF MAXWELL EQUATION
C     ICOMP=2 : THETA COMPONENT OF MAXWELL EQUATION
C     ICOMP=3 : PHI COMPONENT OF MAXWELL EQUATION
C
      ICOMP=MOD(IGD-1,3)+1
C
      IGD1=(IGD-ICOMP)/3+1
      MLXD=MOD(IGD1-1,MDSIZ)+1
C
      IGD2=(IGD1-MLXD)/MDSIZ+1
      NKXD=MOD(IGD2-1,NDSIZ)+1
C
      NR=(IGD2-NKXD)/NDSIZ+1
C
      IF(NR.NE.NRP) THEN
         NRP=NR
         CALL WMSETM_M(NR,CEMP)
         CALL WMSETM_V(NR,CFVP)
         CALL WMSETM_B(NR,CEMP,CFVP)
      ENDIF
C
      DO MB=1,2*MBND-1
         CA(MB)=CEMP(MB,NKXD,MLXD,ICOMP)
      ENDDO
      CB=CFVP(NKXD,MLXD,ICOMP)
C      WRITE(6,*) NR,NKXD,MLXD,ICOMP,CB
C
      RETURN
      END
C
C     ****** ASSEMBLE TOTAL ELEMENT COEFFICIENT MATRIX ******
C
      SUBROUTINE WMSETM_M(NR,CEMP)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CEMP(MBNDM,NDM,MDM,3)
      DIMENSION CGM12(MDM,NDM),CGM13(MDM,NDM)
      DIMENSION CGMH12(MDM,NDM),CGMH13(MDM,NDM)
      DIMENSION CGMH22(MDM,NDM),CGMH23(MDM,NDM),CGMH33(MDM,NDM)
      DIMENSION CGC11(MDM,NDM),CGC12(MDM,NDM),CGC13(MDM,NDM)
      DIMENSION CGPH22(MDM,NDM),CGPH23(MDM,NDM),CGPH33(MDM,NDM)
      DIMENSION CGP12(MDM,NDM),CGP13(MDM,NDM)
      DIMENSION CDVM(3,3),CDVC(3,3),CDVP(3,3)
C
         IF(NR.EQ.1) THEN
C
            CALL WMSETF(1,0)
C
            DO KDX=1,KDSIZ
               DO LDX=1,LDSIZ
                  CGF11(LDX,KDX,2)=CGF11(LDX,KDX,3)
                  CGF12(LDX,KDX,2)=CGF12(LDX,KDX,3)
                  CGF13(LDX,KDX,2)=CGF13(LDX,KDX,3)
                  CGF22(LDX,KDX,2)=CGF22(LDX,KDX,3)
                  CGF23(LDX,KDX,2)=CGF23(LDX,KDX,3)
                  CGF33(LDX,KDX,2)=CGF33(LDX,KDX,3)
               ENDDO
            ENDDO
C
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ
            DO J=1,3
            DO I=1,3
               CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
C
            CALL WMSETF(2,0)
C
         ENDIF
C
         DO KDX=1,KDSIZ
         DO LDX=1,LDSIZ
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
C
         DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ
         DO MDX=1,MDSIZ
         DO LDX=1,LDSIZ
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
C
         IF(NR.LT.NRMAX) CALL WMSETF(NR+2,0)
C
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
C
         DRHOM =XRHOC-XRHOM
         IF(NR.EQ.NRMAX) THEN
            DRHOP =XRHOC-XRHOM
         ELSE
            DRHOP =XRHOP-XRHOC
         ENDIF
         DRHOPM=DRHOM+DRHOP
C
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
C
         FMHM=0.5D0
         FMHC=0.5D0
         FPHC=0.5D0
         FPHP=0.5D0
         DFMHM=-1.0D0/DRHOM   /DPSIPDRHOMH
         DFMHC= 1.0D0/DRHOM   /DPSIPDRHOMH
C
CCC         DFCM=  -DRHOP       /(DRHOM*DRHOPM)/DPSIPDRHOC
         DFCM= -(DRHOP*DPSIPDRHOPH)/(DRHOM*DRHOPM)   !
     &          /(DPSIPDRHOMH*DPSIPDRHOC)   !
CCC         DFCC=  (DRHOP-DRHOM)/(DRHOP*DRHOM )/DPSIPDRHOC
         DFCC=  (DRHOP*DPSIPDRHOPH-DRHOM*DPSIPDRHOMH)/(DRHOP*DRHOM)   !
     &          /(DPSIPDRHOPH*DPSIPDRHOMH)   !
CCC         DFCP=   DRHOM       /(DRHOP*DRHOPM)/DPSIPDRHOC
         DFCP=  (DRHOM*DPSIPDRHOMH)/(DRHOP*DRHOPM) !
     &          /(DPSIPDRHOPH*DPSIPDRHOC)   !
         DDFMHM= 2.D0/(DRHOM*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOMH)
         DDFMHC=-2.D0/(DRHOM*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOMH)
         DDFPHC=-2.D0/(DRHOP*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOPH)
         DDFPHP= 2.D0/(DRHOP*DRHOPM)/(DPSIPDRHOC*DPSIPDRHOPH)
C
CCC         FCMH  = DRHOP/DRHOPM
         FCMH  = (DRHOP*DPSIPDRHOPH)/(DRHOPM*DPSIPDRHOC)   !
CCC         FCPH  = DRHOM/DRHOPM
         FCPH  = (DRHOM*DPSIPDRHOMH)/(DRHOPM*DPSIPDRHOC)   !
         DFCMH =-2.D0 /DRHOPM /DPSIPDRHOC
         DFCPH = 2.D0 /DRHOPM /DPSIPDRHOC
C
         FACT1M=XRHOMH/XRHOM
         FACT1C=XRHOMH/XRHOC
         FACT1P=XRHOPH/XRHOC
C
         FACT2M=XRHOMH/XRHOM
         FACT2C=XRHOMH/XRHOC
         FACT2P=XRHOPH/XRHOC
C
         FACT3M=XRHOMH/XRHOM
         FACT3C=XRHOMH/XRHOC
         FACT3P=XRHOPH/XRHOC
C
         DO KDX=1,KDSIZ
         DO LDX=1,LDSIZ
C
            CGM12(LDX,KDX)=       CGF12(LDX,KDX,1)
            CGM13(LDX,KDX)=       CGF13(LDX,KDX,1)/XRHOM
C
            CGMH12(LDX,KDX)= FMHM*CGF12(LDX,KDX,1)
     &                      +FMHC*CGF12(LDX,KDX,2)
            CGMH13(LDX,KDX)=(FMHM*CGF13(LDX,KDX,1)
     &                      +FMHC*CGF13(LDX,KDX,2))/XRHOMH
            CGMH22(LDX,KDX)=(FMHM*CGF22(LDX,KDX,1)
     &                      +FMHC*CGF22(LDX,KDX,2))*XRHOMH**2
            CGMH23(LDX,KDX)=(FMHM*CGF23(LDX,KDX,1)
     &                      +FMHC*CGF23(LDX,KDX,2))*XRHOMH
            CGMH33(LDX,KDX)= FMHM*CGF33(LDX,KDX,1)
     &                      +FMHC*CGF33(LDX,KDX,2)
C
            CGC11(LDX,KDX)=       CGF11(LDX,KDX,2)/XRHOC**2
            CGC12(LDX,KDX)=       CGF12(LDX,KDX,2)
            CGC13(LDX,KDX)=       CGF13(LDX,KDX,2)/XRHOC
C
            CGPH22(LDX,KDX)=(FPHC*CGF22(LDX,KDX,2)
     &                      +FPHP*CGF22(LDX,KDX,3))*XRHOPH**2
            CGPH23(LDX,KDX)=(FPHC*CGF23(LDX,KDX,2)
     &                      +FPHP*CGF23(LDX,KDX,3))*XRHOPH
            CGPH33(LDX,KDX)= FPHC*CGF33(LDX,KDX,2)
     &                      +FPHP*CGF33(LDX,KDX,3)
C
            CGP12(LDX,KDX)=       CGF12(LDX,KDX,3)
            CGP13(LDX,KDX)=       CGF13(LDX,KDX,3)/XRHOP
         ENDDO
         ENDDO
C
C            IF(NR.EQ.1) THEN
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF11(1,1,1),CGF11(1,1,2),CGF11(1,1,3)
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF12(1,1,1),CGF12(1,1,2),CGF12(1,1,3)
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF13(1,1,1),CGF13(1,1,2),CGF13(1,1,3)
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF22(1,1,1),CGF22(1,1,2),CGF22(1,1,3)
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF23(1,1,1),CGF23(1,1,2),CGF23(1,1,3)
C               WRITE(6,'(1P6E12.4)') 
C     &                 CGF33(1,1,1),CGF33(1,1,2),CGF33(1,1,3)
C            ENDIF
C
C        ND : (n - n0) / Np
C        KD : n'              = (k - n) / Np
C        NN : n               = n0 +  ND       * Np
C        NK : k = n + n' * Np = n0 + (ND + KD) * Np
C
C        MD : m - m0
C        LD : m'              = l - m
C        MM : m               = m0 +  MD
C        ML : l = m + m'      = m0 + (MD + LD)
C
         DO LBAND=1,2*MBND-1
         DO NKX=1,NDSIZ
         DO MLX=1,MDSIZ
            CEMP(LBAND,NKX,MLX,1)=0.D0
            CEMP(LBAND,NKX,MLX,2)=0.D0
            CEMP(LBAND,NKX,MLX,3)=0.D0
         ENDDO
         ENDDO
         ENDDO
C
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO NKD=NDMIN,NDMAX
            NKX=NKD-NDMIN+1
            KD=NKD-ND
            IF(MODELK.EQ.0.OR.
     &         (KD.GE.KDMIN.AND.KD.LE.KDMAX)) THEN
            KDX=MOD(KD-KDMIN+2*KDSIZ,KDSIZ)+1
            NN=NPH0+NHC*ND
            NK=NPH0+NHC*NKD
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
         DO MLD=MDMIN,MDMAX
            MLX=MLD-MDMIN+1
            LD=MLD-MD
            IF(MODELK.EQ.0.OR.
     &         (LD.GE.LDMIN.AND.LD.LE.LDMAX)) THEN
            LDX=MOD(LD-LDMIN+2*LDSIZ,LDSIZ)+1
            MM=NTH0+MD
            ML=NTH0+MLD
C
            DO J=1,3
            DO I=1,3
               CDVM(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,1)
               CDVC(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,2)
               CDVP(I,J)=CGD(I,J,LDX,MLX,KDX,NKX,3)
            ENDDO
            ENDDO
C
C            IF(NR.EQ.1) THEN
C               DO J=1,3
C               DO I=1,3
C                  WRITE(6,'(1P6E12.4)') CDVM(I,J),CDVC(I,J),CDVP(I,J)
C               ENDDO
C               ENDDO
C            ENDIF
C     
            ID=3*MDSIZ*NDSIZ
C
C     --- R COMPONENT OF MAXWELL EQUATION ---
C
            LBND=MBND-3*KD*MDSIZ-3*LD-1
C
            CEMP(LBND+1    ,NKX,MLX,1)
     &                         =CEMP(LBND+1   ,NKX,MLX,1)
     &                         +(
     &                         -ML*NN*CGMH23(LDX,KDX)
     &                         +ML*MM*CGMH33(LDX,KDX)
     &                         +NK*NN*CGMH22(LDX,KDX)
     &                         -NK*MM*CGMH23(LDX,KDX)
     &                         +0.5D0*CDVM(1,1)*FACT1M
     &                         +0.5D0*CDVC(1,1)*FACT1C
     &                          )/XRHOMH
C
            CEMP(LBND+2-ID,NKX,MLX,1)
     &                         =CEMP(LBND+2-ID,NKX,MLX,1)
     &                         +(
     &                         +ML*NN*CGMH13(LDX,KDX)*FMHM
     &                         -NK*NN*CGMH12(LDX,KDX)*FMHM
     &                         +CI*ML*CGMH33(LDX,KDX)*DFMHM
     &                         -CI*NK*CGMH23(LDX,KDX)*DFMHM
     &                         +CDVM(1,2)*FMHM
     &                          )*XRHOM
C
            CEMP(LBND+2   ,NKX,MLX,1)
     &                         =CEMP(LBND+2   ,NKX,MLX,1)
     &                         +(
     &                         +ML*NN*CGMH13(LDX,KDX)*FMHC
     &                         -NK*NN*CGMH12(LDX,KDX)*FMHC
     &                         +CI*ML*CGMH33(LDX,KDX)*DFMHC
     &                         -CI*NK*CGMH23(LDX,KDX)*DFMHC
     &                         +CDVC(1,2)*FMHC
     &                          )*XRHOC
C
            CEMP(LBND+3-ID,NKX,MLX,1)
     &                         =CEMP(LBND+3-ID,NKX,MLX,1)
     &                         +(
     &                         -ML*MM*CGMH13(LDX,KDX)*FMHM
     &                         +NK*MM*CGMH12(LDX,KDX)*FMHM
     &                         -CI*ML*CGMH23(LDX,KDX)*DFMHM
     &                         +CI*NK*CGMH22(LDX,KDX)*DFMHM
     &                         +CDVM(1,3)*FMHM
     &                          )
C
            CEMP(LBND+3   ,NKX,MLX,1)
     &                         =CEMP(LBND+3   ,NKX,MLX,1)
     &                         +(
     &                         -ML*MM*CGMH13(LDX,KDX)*FMHC
     &                         +NK*MM*CGMH12(LDX,KDX)*FMHC
     &                         -CI*ML*CGMH23(LDX,KDX)*DFMHC
     &                         +CI*NK*CGMH22(LDX,KDX)*DFMHC
     &                         +CDVC(1,3)*FMHC
     &                          )
C
C     --- THETA COMPONENT OF MAXWELL EQUATION ---
C
            LBND=MBND-3*KD*MDSIZ-3*LD-2
C
            CEMP(LBND+1   ,NKX,MLX,2)
     &                         =CEMP(LBND+1   ,NKX,MLX,2)
     &                         +(
     &                         -NK*NN*CGC12(LDX,KDX)*FCMH
     &                         +NK*MM*CGC13(LDX,KDX)*FCMH
     &                         -CI*NN*CGMH23(LDX,KDX)*DFCMH
     &                         +CI*MM*CGMH33(LDX,KDX)*DFCMH
     &                         +CDVC(2,1)*FCMH*FACT2C
     &                          )/XRHOMH
C
            CEMP(LBND+1+ID,NKX,MLX,2)
     &                         =CEMP(LBND+1+ID,NKX,MLX,2)
     &                         +(
     &                         -NK*NN*CGC12(LDX,KDX)*FCPH
     &                         +NK*MM*CGC13(LDX,KDX)*FCPH
     &                         -CI*NN*CGPH23(LDX,KDX)*DFCPH
     &                         +CI*MM*CGPH33(LDX,KDX)*DFCPH
     &                         +CDVC(2,1)*FCPH*FACT2P
     &                          )/XRHOPH
C
            CEMP(LBND+2-ID,NKX,MLX,2)
     &                         =CEMP(LBND+2-ID,NKX,MLX,2)
     &                         +(
     &                         +CI*NK*CGC13(LDX,KDX)*DFCM
     &                         +CI*NN*CGM13(LDX,KDX)*DFCM
     &                         -      CGMH33(LDX,KDX)*DDFMHM
     &                          )*XRHOM
C
            CEMP(LBND+2   ,NKX,MLX,2)
     &                         =CEMP(LBND+2   ,NKX,MLX,2)
     &                         +(
     &                         +NK*NN*CGC11(LDX,KDX)
     &                         +CI*NK*CGC13(LDX,KDX)*DFCC
     &                         +CI*NN*CGC13(LDX,KDX)*DFCC
     &                         -      CGMH33(LDX,KDX)*DDFMHC
     &                         -      CGPH33(LDX,KDX)*DDFPHC
     &                         +CDVC(2,2)
     &                          )*XRHOC
C
            CEMP(LBND+2+ID,NKX,MLX,2)
     &                         =CEMP(LBND+2+ID,NKX,MLX,2)
     &                         +(
     &                         +CI*NK*CGC13(LDX,KDX)*DFCP
     &                         +CI*NN*CGP13(LDX,KDX)*DFCP
     &                         -      CGPH33(LDX,KDX)*DDFPHP
     &                          )*XRHOP
C
            CEMP(LBND+3-ID,NKX,MLX,2)
     &                         =CEMP(LBND+3-ID,NKX,MLX,2)
     &                         +(
     &                         -CI*NK*CGC12(LDX,KDX)*DFCM
     &                         -CI*MM*CGM13(LDX,KDX)*DFCM
     &                         +      CGMH23(LDX,KDX)*DDFMHM
     &                          )
C
            CEMP(LBND+3   ,NKX,MLX,2)
     &                         =CEMP(LBND+3   ,NKX,MLX,2)
     &                         +(
     &                         -NK*MM*CGC11(LDX,KDX)
     &                         -CI*NK*CGC12(LDX,KDX)*DFCC
     &                         -CI*MM*CGC13(LDX,KDX)*DFCC
     &                         +      CGMH23(LDX,KDX)*DDFMHC
     &                         +      CGPH23(LDX,KDX)*DDFPHC
     &                         +CDVC(2,3)
     &                          )
C
            CEMP(LBND+3+ID,NKX,MLX,2)
     &                         =CEMP(LBND+3+ID,NKX,MLX,2)
     &                         +(
     &                         -CI*NK*CGC12(LDX,KDX)*DFCP
     &                         -CI*MM*CGP13(LDX,KDX)*DFCP
     &                         +      CGPH23(LDX,KDX)*DDFPHP
     &                          )
C
C     --- PHI COMPONENT OF MAXWELL EQUATION ---
C
            LBND=MBND-3*KD*MDSIZ-3*LD-3
C
            CEMP(LBND+1   ,NKX,MLX,3)
     &                         =CEMP(LBND+1   ,NKX,MLX,3)
     &                         +(
     &                         +ML*NN*CGC12(LDX,KDX)*FCMH
     &                         -ML*MM*CGC13(LDX,KDX)*FCMH
     &                         +CI*NN*CGMH22(LDX,KDX)*DFCMH
     &                         -CI*MM*CGMH23(LDX,KDX)*DFCMH
     &                         +CDVC(3,1)*FCMH*FACT3C
     &                          )/XRHOMH
C
            CEMP(LBND+1+ID,NKX,MLX,3)
     &                         =CEMP(LBND+1+ID,NKX,MLX,3)
     &                         +(
     &                         +ML*NN*CGC12(LDX,KDX)*FCPH
     &                         -ML*MM*CGC13(LDX,KDX)*FCPH
     &                         +CI*NN*CGPH22(LDX,KDX)*DFCPH
     &                         -CI*MM*CGPH23(LDX,KDX)*DFCPH
     &                         +CDVC(3,1)*FCPH*FACT3P
     &                          )/XRHOPH
C
            CEMP(LBND+2-ID,NKX,MLX,3)
     &                         =CEMP(LBND+2-ID,NKX,MLX,3)
     &                         +(
     &                         -CI*ML*CGC13(LDX,KDX)*DFCM
     &                         -CI*NN*CGM12(LDX,KDX)*DFCM
     &                         +      CGMH23(LDX,KDX)*DDFMHM
     &                          )*XRHOM
C
            CEMP(LBND+2   ,NKX,MLX,3)
     &                         =CEMP(LBND+2   ,NKX,MLX,3)
     &                         +(
     &                         -ML*NN*CGC11(LDX,KDX)
     &                         -CI*ML*CGC13(LDX,KDX)*DFCC
     &                         -CI*NN*CGC12(LDX,KDX)*DFCC
     &                         +      CGMH23(LDX,KDX)*DDFMHC
     &                         +      CGPH23(LDX,KDX)*DDFPHC
     &                         +CDVC(3,2)
     &                          )*XRHOC
C
            CEMP(LBND+2+ID,NKX,MLX,3)
     &                         =CEMP(LBND+2+ID,NKX,MLX,3)
     &                         +(
     &                         -CI*ML*CGC13(LDX,KDX)*DFCP
     &                         -CI*NN*CGP12(LDX,KDX)*DFCP
     &                         +      CGPH23(LDX,KDX)*DDFPHP
     &                          )*XRHOP
C
            CEMP(LBND+3-ID,NKX,MLX,3)
     &                         =CEMP(LBND+3-ID,NKX,MLX,3)
     &                         +(
     &                         +CI*ML*CGC12(LDX,KDX)*DFCM
     &                         +CI*MM*CGM12(LDX,KDX)*DFCM
     &                         -      CGMH22(LDX,KDX)*DDFMHM
     &                          )
C
            CEMP(LBND+3   ,NKX,MLX,3)
     &                         =CEMP(LBND+3   ,NKX,MLX,3)
     &                         +(
     &                         +ML*MM*CGC11(LDX,KDX)
     &                         +CI*ML*CGC12(LDX,KDX)*DFCC
     &                         +CI*MM*CGC12(LDX,KDX)*DFCC
     &                         -      CGMH22(LDX,KDX)*DDFMHC
     &                         -      CGPH22(LDX,KDX)*DDFPHC
     &                         +CDVC(3,3)
     &                          )
C
            CEMP(LBND+3+ID,NKX,MLX,3)
     &                         =CEMP(LBND+3+ID,NKX,MLX,3)
     &                         +(
     &                         +CI*ML*CGC12(LDX,KDX)*DFCP
     &                         +CI*MM*CGP12(LDX,KDX)*DFCP
     &                         -      CGPH22(LDX,KDX)*DDFPHP
     &                          )
            ENDIF
         ENDDO
         ENDDO
            ENDIF
         ENDDO
         ENDDO
      RETURN
      END
C
C     ****** ASSEMBLE TOTAL ELEMENT FREE VECTOR ******
C
      SUBROUTINE WMSETM_V(NR,CFVP)
C
      USE plprof,ONLY: pl_prof2
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CFVP(NDM,MDM,3)
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
C
         DO MDX=1,MDSIZ
            DO NDX=1,NDSIZ
               CFVP(NDX,MDX,1)=0.D0
               CFVP(NDX,MDX,2)=0.D0
               CFVP(NDX,MDX,3)=0.D0
            ENDDO
         ENDDO
C
         IF(MODEEG.EQ.0) THEN
            DO NRI=1,NRMAX 
               IF(XR(NRI)/RD.LT.1.D0) NRANT=NRI
            ENDDO
C
            IF(NR+1.GE.NRANT) THEN
               CW=2.D0*PI*CRF*1.D6
               CC=CI*CW*RMU0
               DPH=2.D0*PI/NHHMAX
               DTH=2.D0*PI/NTHMAX
C
               IF(NR+1.EQ.NRANT.OR.NR+1.EQ.NRANT+1) THEN
                  XRHO1=XRHO(NRANT)
                  XRHO2=XRHO(NRANT+1)
                  DRHO=XRHO2-XRHO1
                  XRHOC=0.5D0*(XRHO2+XRHO1)
C
                  IF(MODELG.EQ.3) THEN
                     QPC=0.5D0*(QPS(NRANT)+QPS(NRANT+1))
                     DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
                  ELSE
                     DPSIPDRHOC =2.D0*PSIPA*XRHOC
                  ENDIF
C
                  FACTM=(XRHO2-RD/RA)/DRHO
                  FACTP=(RD/RA-XRHO1)/DRHO
C
                  DO ND=NDMIN,NDMAX
                     NDX=ND-NDMIN+1
                     NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*FACTM*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*FACTM*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)
     &                                         *DPSIPDRHOC*DRHO/XRHOC
                     CJTHP=CC*CJANT(2,MDX,NDX)*FACTP*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DPH)
                     CJPHP=CC*CJANT(3,MDX,NDX)*FACTP*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DTH)
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
C                     WRITE(6,*) 'CFVP(',NDX,MDX,1,')',CFVP(NDX,MDX,1)
C                     WRITE(6,*) 'CFVP(',NDX,MDX,2,')',CFVP(NDX,MDX,2)
C                     WRITE(6,*) 'CFVP(',NDX,MDX,3,')',CFVP(NDX,MDX,3)
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
C
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
C
                  DO ND=NDMIN,NDMAX
                     NDX=ND-NDMIN+1
                     NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*XRHOC
     &                             /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*XRHOC
     &                             /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DPSIPDRHOC*DRHO
     &                             /XRHOC
                     CFVP(NDX,MDX,1)=CJR
C                     WRITE(6,*) 'CFVP(',NDX,MDX,1,')',CFVP(NDX,MDX,1)
C                     WRITE(6,*) 'CFVP(',NDX,MDX,2,')',CFVP(NDX,MDX,2)
C                     WRITE(6,*) 'CFVP(',NDX,MDX,3,')',CFVP(NDX,MDX,3)
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
C
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
C
C            CALL WMCDEN(NR+1,RN,RTPR,RTPP,RU)
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
C               CFVP(NDX,MDX,1)=CJR
               CFVP(NDX,MDX,1)=0.D0
               CFVP(NDX,MDX,2)=CJTHM
               CFVP(NDX,MDX,3)=CJPHM
            ENDDO
            ENDDO
         ENDIF
C
      RETURN
      END
C
C     ****** SET BOUNDARY CONDITION ******
C
      SUBROUTINE WMSETM_B(NR,CEMP,CFVP)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CEMP(MBNDM,NDM,MDM,3),CFVP(NDM,MDM,3)
C
         DRHO1=(XRHO(2)-XRHO(1))**2
         DRHO2=(XRHO(3)-XRHO(1))**2
         A1= DRHO2/(DRHO2-DRHO1)
         A2=-DRHO1/(DRHO2-DRHO1)
C
C        ****** R=0 ******
C
         IF(NR.EQ.1) THEN
            ID=3*MDSIZ*NDSIZ         
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO NKD=NDMIN,NDMAX
               NKX=NKD-NDMIN+1
               KD=NKD-ND
               IF(MODELK.EQ.0.OR.
     &              (KD.GE.KDMIN.AND.KD.LE.KDMAX)) THEN
C                  KDX=MOD(KD-KDMIN+2*KDSIZ,KDSIZ)+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
            DO MLD=MDMIN,MDMAX
               MLX=MLD-MDMIN+1
               LD=MLD-MD
               IF(MODELK.EQ.0.OR.
     &              (LD.GE.LDMIN.AND.LD.LE.LDMAX)) THEN
C                  LDX=MOD(LD-LDMIN+2*LDSIZ,LDSIZ)+1
                  MM=NTH0+MD
C
C        ****** EPH'(0) = 0 FOR MM.EQ.0 ******
C
                  IF(MM.EQ.0) THEN
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-1
C
                     CEMP(LBND   +3,NKX,MLX,1)
     &                               =CEMP(LBND   +3,NKX,MLX,1)
     &                               +CEMP(LBND-ID+3,NKX,MLX,1)*A1
                     CEMP(LBND+ID+3,NKX,MLX,1)
     &                               =CEMP(LBND+ID+3,NKX,MLX,1)
     &                               +CEMP(LBND-ID+3,NKX,MLX,1)*A2
                     CEMP(LBND-ID+3,NKX,MLX,1)=0.D0
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-2
C
                     CEMP(LBND   +3,NKX,MLX,2)
     &                               =CEMP(LBND   +3,NKX,MLX,2)
     &                               +CEMP(LBND-ID+3,NKX,MLX,2)*A1
                     CEMP(LBND+ID+3,NKX,MLX,2)
     &                               =CEMP(LBND+ID+3,NKX,MLX,2)
     &                               +CEMP(LBND-ID+3,NKX,MLX,2)*A2
                     CEMP(LBND-ID+3,NKX,MLX,2)=0.D0
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-3
C
                     CEMP(LBND   +3,NKX,MLX,3)
     &                               =CEMP(LBND   +3,NKX,MLX,3)
     &                               +CEMP(LBND-ID+3,NKX,MLX,3)*A1
                     CEMP(LBND+ID+3,NKX,MLX,3)
     &                               =CEMP(LBND+ID+3,NKX,MLX,3)
     &                               +CEMP(LBND-ID+3,NKX,MLX,3)*A2
                     CEMP(LBND-ID+3,NKX,MLX,3)=0.D0
C
C        ****** ETH'(0) = 0 FOR ABS(MM).EQ.1 ******
C
                  ELSEIF(ABS(MM).EQ.1) THEN
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-1
C
                     CEMP(LBND   +2,NKX,MLX,1)
     &                               =CEMP(LBND   +2,NKX,MLX,1)
     &                               +CEMP(LBND-ID+2,NKX,MLX,1)*A1
                     CEMP(LBND+ID+2,NKX,MLX,1)
     &                               =CEMP(LBND+ID+2,NKX,MLX,1)
     &                               +CEMP(LBND-ID+2,NKX,MLX,1)*A2
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-2
C
                     CEMP(LBND   +2,NKX,MLX,2)
     &                               =CEMP(LBND   +2,NKX,MLX,2)
     &                               +CEMP(LBND-ID+2,NKX,MLX,2)*A1
                     CEMP(LBND+ID+2,NKX,MLX,2)
     &                               =CEMP(LBND+ID+2,NKX,MLX,2)
     &                               +CEMP(LBND-ID+2,NKX,MLX,2)*A2
C
                     LBND=MBND-3*KD*MDSIZ-3*LD-3
C
                     CEMP(LBND   +2,NKX,MLX,3)
     &                               =CEMP(LBND   +2,NKX,MLX,3)
     &                               +CEMP(LBND-ID+2,NKX,MLX,3)*A1
                     CEMP(LBND+ID+2,NKX,MLX,3)
     &                               =CEMP(LBND+ID+2,NKX,MLX,3)
     &                               +CEMP(LBND-ID+2,NKX,MLX,3)*A2
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
               ENDIF
            ENDDO
            ENDDO
         ENDIF
C
C        ****** ETH = 0, EPH =0 AT R=RA ******
C
         IF(NR.EQ.NRMAX) THEN
C
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               DO MB=1,2*MBND-1
                  CEMP(MB,NDX,MDX,1)= 0.D0
                  CEMP(MB,NDX,MDX,2)= 0.D0
                  CEMP(MB,NDX,MDX,3)= 0.D0
               ENDDO
               CEMP(MBND,NDX,MDX,1)= 1.D0
               CEMP(MBND-3*MDSIZ*NDSIZ,NDX,MDX,1)= -1.D0
               CEMP(MBND,NDX,MDX,2)= 1.D0
               CEMP(MBND,NDX,MDX,3)= 1.D0
               CFVP(NDX,MDX,1)= 0.D0
               CFVP(NDX,MDX,2)= CEWALL(MDX,NDX,2)
               CFVP(NDX,MDX,3)= CEWALL(MDX,NDX,3)
               WRITE(22,'(A,2I8,4ES12.4)') 
     &               'CFVP:',MDX,NDX,CEWALL(MDX,NDX,2),CEWALL(MDX,NDX,3)
            ENDDO
            ENDDO
C
         ENDIF
C
C     ****** ELIMINATE EFLD AT MDMAX ******
C
         IF(MDSIZ.GT.1) THEN
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
               MD=MDMAX
               MDX=MD-MDMIN+1
               DO MB=1,2*MBND-1
                  CEMP(MB,NDX,MDX,1)= 0.D0
                  CEMP(MB,NDX,MDX,2)= 0.D0
                  CEMP(MB,NDX,MDX,3)= 0.D0
               ENDDO
               CEMP(MBND,NDX,MDX,1)= 1.D0
               CEMP(MBND,NDX,MDX,2)= 1.D0
               CEMP(MBND,NDX,MDX,3)= 1.D0
               CFVP(NDX,MDX,1)= 0.D0
               CFVP(NDX,MDX,2)= 0.D0
               CFVP(NDX,MDX,3)= 0.D0
            ENDDO
         ENDIF
C
C     ****** ELIMINATE EFLD AT NDMAX ******
C
         IF(NDSIZ.GT.1) THEN
            ND=NDMAX
            NDX=ND-NDMIN+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               DO MB=1,2*MBND-1
                  CEMP(MB,NDX,MDX,1)= 0.D0
                  CEMP(MB,NDX,MDX,2)= 0.D0
                  CEMP(MB,NDX,MDX,3)= 0.D0
               ENDDO
               CEMP(MBND,NDX,MDX,1)= 1.D0
               CEMP(MBND,NDX,MDX,2)= 1.D0
               CEMP(MBND,NDX,MDX,3)= 1.D0
               CFVP(NDX,MDX,1)= 0.D0
               CFVP(NDX,MDX,2)= 0.D0
               CFVP(NDX,MDX,3)= 0.D0
            ENDDO
         ENDIF
C
      RETURN
      END

