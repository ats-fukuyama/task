C     $Id$
C
C     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
C
      COMPLEX*16 FUNCTION CFDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
C
      USE plcomm
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CDET(3,3)
C
      CALL DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      CDET11=CDET(1,1)
      CDET12=CDET(1,2)
      CDET13=CDET(1,3)
      CDET21=CDET(2,1)
      CDET22=CDET(2,2)
      CDET23=CDET(2,3)
      CDET31=CDET(3,1)
      CDET32=CDET(3,2)
      CDET33=CDET(3,3)
C
      CFDISP =CDET33*(CDET11*CDET22-CDET12*CDET21)
     &       +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11
     &       +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
C
      RETURN
      END
C
C     ****** CALCULATE REAL PART OF DETERMINANT OF DISPERSION TENSOR ******
C
      COMPLEX*16 FUNCTION CFDISPR(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
C
      USE plcomm
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CDET(3,3)
C
      CALL DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      CDET11=0.5D0*(CDET(1,1)+DCONJG(CDET(1,1)))
      CDET12=0.5D0*(CDET(1,2)+DCONJG(CDET(2,1)))
      CDET13=0.5D0*(CDET(1,3)+DCONJG(CDET(3,1)))
      CDET21=0.5D0*(CDET(2,1)+DCONJG(CDET(1,2)))
      CDET22=0.5D0*(CDET(2,2)+DCONJG(CDET(2,2)))
      CDET23=0.5D0*(CDET(2,3)+DCONJG(CDET(3,2)))
      CDET31=0.5D0*(CDET(3,1)+DCONJG(CDET(1,3)))
      CDET32=0.5D0*(CDET(3,2)+DCONJG(CDET(2,3)))
      CDET33=0.5D0*(CDET(3,3)+DCONJG(CDET(3,3)))
C
      CFDISPR=CDET33*(CDET11*CDET22-CDET12*CDET21)
     &       +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11
     &       +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
C
      RETURN
      END
C
C     ****** CALCULATE DISPERSION TENSOR ******
C
      SUBROUTINE DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      USE plcomm
      USE pllocal
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CDTNS(3,3),CDET(3,3)
C
      CW=2.D0*PI*1.D6*CRF
      CL2=VC*VC/(CW*CW)
C
      CALL DPEXEC(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,0,CDTNS)
C
      CDET(1,1)=-CL2*(CKY*CKY+CKZ*CKZ)+CDTNS(1,1)
      CDET(1,2)= CL2* CKX*CKY         +CDTNS(1,2)
      CDET(1,3)= CL2* CKX*CKZ         +CDTNS(1,3)
      CDET(2,1)= CL2* CKY*CKX         +CDTNS(2,1)
      CDET(2,2)=-CL2*(CKZ*CKZ+CKX*CKX)+CDTNS(2,2)
      CDET(2,3)= CL2* CKY*CKZ         +CDTNS(2,3)
      CDET(3,1)= CL2* CKZ*CKX         +CDTNS(3,1)
      CDET(3,2)= CL2* CKZ*CKY         +CDTNS(3,2)
      CDET(3,3)=-CL2*(CKX*CKX+CKY*CKY)+CDTNS(3,3)
C
      RETURN
      END
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE DPEXEC(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,NS,CDTNS)
C
      USE plcomm
      USE pllocal
      USE plprof
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CDTNS(3,3)
C
      CALL PL_MAG_OLD(XPOS,YPOS,ZPOS,RHON)
      IF(MODELG.LE.1.OR.MODELG.GT.10) THEN
         CALL PL_PROF3D_OLD(XPOS,YPOS,ZPOS)
      ELSE
         CALL PL_PROF_OLD(RHON)
      END IF
C
      CW=2.D0*PI*1.D6*CRF
      CKPR=BNX*CKX+BNY*CKY+BNZ*CKZ
      IF(ABS(CKPR).LE.1.D-8) CKPR=1.D-8
      CKPP=SQRT(CKX**2+CKY**2+CKZ**2-CKPR**2)
C
      CALL DPCALC(CW,CKPR,CKPP,NS,CDTNS)
C
      IF(ABS(CKPP).EQ.0.D0) THEN
         BNXY=SQRT(BNX*BNX+BNY*BNY)
         IF(BNXY.EQ.0.D0) THEN
            CU11= 1.D0
            CU12= 0.D0
            CU13= 0.D0
            CU21= 0.D0
            CU22= 1.D0
            CU23= 0.D0
            CU31= 1.D0
            CU32= 0.D0
            CU33= 0.D0
         ELSE
            CU11=-BNZ*BNX/BNXY
            CU12=-BNZ*BNY/BNXY
            CU13= BNXY
            CU21= BNY/BNXY
            CU22=-BNX/BNXY
            CU23= 0.D0
            CU31= BNX
            CU32= BNY
            CU33= BNZ
         ENDIF
      ELSE
         CU11=( (1.D0-BNX**2)*CKX-BNX*BNY*CKY-BNX*BNZ*CKZ)/CKPP
         CU12=(-BNY*BNX*CKX+(1.D0-BNY**2)*CKY-BNY*BNZ*CKZ)/CKPP
         CU13=(-BNZ*BNX*CKX-BNZ*BNY*CKY+(1.D0-BNZ**2)*CKZ)/CKPP
         CU21=(BNY*CKZ-BNZ*CKY)/CKPP
         CU22=(BNZ*CKX-BNX*CKZ)/CKPP
         CU23=(BNX*CKY-BNY*CKX)/CKPP
         CU31=BNX
         CU32=BNY
         CU33=BNZ
      ENDIF
C
      CDET11=CDTNS(1,1)*CU11+CDTNS(1,2)*CU21+CDTNS(1,3)*CU31
      CDET12=CDTNS(1,1)*CU12+CDTNS(1,2)*CU22+CDTNS(1,3)*CU32
      CDET13=CDTNS(1,1)*CU13+CDTNS(1,2)*CU23+CDTNS(1,3)*CU33
      CDET21=CDTNS(2,1)*CU11+CDTNS(2,2)*CU21+CDTNS(2,3)*CU31
      CDET22=CDTNS(2,1)*CU12+CDTNS(2,2)*CU22+CDTNS(2,3)*CU32
      CDET23=CDTNS(2,1)*CU13+CDTNS(2,2)*CU23+CDTNS(2,3)*CU33
      CDET31=CDTNS(3,1)*CU11+CDTNS(3,2)*CU21+CDTNS(3,3)*CU31
      CDET32=CDTNS(3,1)*CU12+CDTNS(3,2)*CU22+CDTNS(3,3)*CU32
      CDET33=CDTNS(3,1)*CU13+CDTNS(3,2)*CU23+CDTNS(3,3)*CU33
c
      CDTNS(1,1)=CU11*CDET11+CU21*CDET21+CU31*CDET31
      CDTNS(1,2)=CU11*CDET12+CU21*CDET22+CU31*CDET32
      CDTNS(1,3)=CU11*CDET13+CU21*CDET23+CU31*CDET33
      CDTNS(2,1)=CU12*CDET11+CU22*CDET21+CU32*CDET31
      CDTNS(2,2)=CU12*CDET12+CU22*CDET22+CU32*CDET32
      CDTNS(2,3)=CU12*CDET13+CU22*CDET23+CU32*CDET33
      CDTNS(3,1)=CU13*CDET11+CU23*CDET21+CU33*CDET31
      CDTNS(3,2)=CU13*CDET12+CU23*CDET22+CU33*CDET32
      CDTNS(3,3)=CU13*CDET13+CU23*CDET23+CU33*CDET33
C
      RETURN
      END
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE DPCALC(CW,CKPR,CKPP,NS,CDTNS)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: pl_prof_old
      INCLUDE '../dp/dpcomm.inc'
      COMPLEX(8),INTENT(IN):: CW,CKPR,CKPP
      INTEGER(4),INTENT(IN):: NS
      COMPLEX(8),DIMENSION (3,3),INTENT(OUT):: CDTNS
      COMPLEX(8),DIMENSION(6):: CDISP,CLDISP
C
      IF(NS.EQ.0) THEN
         CDISP(1)=1.D0
         DO I=2,6
            CDISP(I)=0.D0
         ENDDO
         DO NS1=1,NSMAX
            IF(modelp(ns1).EQ.5.OR.
     &         modelp(ns1).EQ.6.OR.
     &         modelp(ns1).eq.15) THEN
               CALL DPCOLD_RKPERP(cw,ckpr,ckppf,ckpps)
!               write(6,'(1P6E12.4)') ckpr,ckppf,ckpps 
               IF(real(ckppf**2).GT.0.d0) THEN
                  ckpp1=ckppf
               ELSE
                  ckpp1=ckpp
               ENDIF
            ELSE
               ckpp1=ckpp
            ENDIF
            CALL DPTENS(CW,CKPR,CKPP1,NS1,CLDISP)
            DO I=1,6
               CDISP(I)=CDISP(I)+CLDISP(I)
            ENDDO
         ENDDO
      ELSE
         IF(modelp(ns).EQ.5.OR.
     &      modelp(ns).eq.6.OR.
     &      modelp(ns).eq.15) THEN
            CALL DPCOLD_RKPERP(cw,ckpr,ckppf,ckpps)
!            write(6,'(1P6E12.4)') ckpr,ckppf,ckpps 
            IF(real(ckppf**2).GT.0.d0) THEN
               ckpp1=ckppf
            ELSE
               ckpp1=ckpp
            ENDIF
         ELSE
            ckpp1=ckpp
         ENDIF
         CALL DPTENS(CW,CKPR,CKPP1,NS,CDISP)
      ENDIF
C

      CDTNS(1,1)= CDISP(1)
      CDTNS(1,2)= CDISP(5)
      CDTNS(1,3)= CDISP(4)
      CDTNS(2,1)=-CDISP(5)
      CDTNS(2,2)= CDISP(1)+CDISP(3)
      CDTNS(2,3)= CDISP(6)
      CDTNS(3,1)= CDISP(4)
      CDTNS(3,2)=-CDISP(6)
      CDTNS(3,3)= CDISP(1)+CDISP(2)
C
      RETURN
      END
C
C     ****** CALCULATE COLD KPERP ******
C
      SUBROUTINE DPCOLD_RKPERP(CW,CKPR,CKPPF,CKPPS)
C
      USE plcomm
      USE pllocal
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CDISP(6),CLDISP(6)
C
      CDISP(1)=1.D0
      DO I=2,6
         CDISP(I)=0.D0
      ENDDO
      CKPP=(0.D0,0.D0)
      CNPR=CKPR*VC/CW

      DO NS1=1,NSMAX
         MODELP_SAVE=MODELP(NS1)
         MODELP(NS1)=0
         CALL DPTENS(CW,CKPR,CKPP,NS1,CLDISP)
         DO I=1,6
            CDISP(I)=CDISP(I)+CLDISP(I)
         ENDDO
         MODELP(NS1)=MODELP_SAVE
      ENDDO

      CCS= CDISP(1)
      CCD= -CI*CDISP(5)
      CCP= CDISP(1)+CDISP(2)

      CCA=CCS
      CCB=(CCP+CCS)*CNPR**2-(CCS**2-CCD**2+CCS*CCP)
      CCC=((CCS-CNPR**2)**2-CCD**2)*CCP
      CCD=SQRT(CCB**2-4.D0*CCA*CCC)
      CKPPA=(-CCB+CCD)/(2.D0*CCA)
      CKPPB=(-CCB-CCD)/(2.D0*CCA)
      RKPPA2=REAL(CKPPA**2)
      RKPPB2=REAL(CKPPB**2)
      IF(RKPPA2.GE.0.D0) THEN
         IF(RKPPB2.GE.0.D0) THEN
            IF(RKPPA2.GT.RKPPB2) THEN
               CTEMP=CKPPA
               CKPPA=CKPPB
               CKPPB=CTEMP
            ENDIF
         ELSE
            CKPPB=0.D0
         ENDIF
      ELSE
         IF(RKPPB2.GE.0.D0) THEN
            CKPPA=CKPPB
            CKPPB=0.D0
         ELSE
            CKPPA=0.D0
            CKPPB=0.D0
         ENDIF
      ENDIF
      CKPPF=SQRT(CKPPA)*CW/VC
      CKPPS=SQRT(CKPPB)*CW/VC
      RETURN
      END
