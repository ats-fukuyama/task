C     $Id$
C
C     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
C
      COMPLEX*16 FUNCTION CFDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
C
      INCLUDE '../dp/dpcomm.h'
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
      INCLUDE '../dp/dpcomm.h'
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
      INCLUDE '../dp/dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CDISP(6),CLDISP(6),CDET(3,3)
C
      CALL PLMAG(XPOS,YPOS,ZPOS,PSIN)
      CALL PLPROF(PSIN)
C
      CW=2.D0*PI*1.D6*CRF
      CKPR=BNX*CKX+BNY*CKY+BNZ*CKZ
      IF(ABS(CKPR).LE.1.D-8) CKPR=1.D-8
      CKPP=SQRT(CKX**2+CKY**2+CKZ**2-CKPR**2)
      IF(ABS(CKPP).LE.1.D-8) CKPP=1.D-8
C
      CL2=VC*VC/(CW*CW)
      CDV11 =1.D0-CL2*(CKY*CKY+CKZ*CKZ)
      CDV12 =     CL2* CKX*CKY
      CDV13 =     CL2* CKX*CKZ
      CDV21 =     CL2* CKY*CKX
      CDV22 =1.D0-CL2*(CKZ*CKZ+CKX*CKX)
      CDV23 =     CL2* CKY*CKZ
      CDV31 =     CL2* CKZ*CKX
      CDV32 =     CL2* CKZ*CKY
      CDV33 =1.D0-CL2*(CKX*CKX+CKY*CKY)
C
      DO I=1,6
         CDISP(I)=0.D0
      ENDDO
C
      DO NS=1,NSMAX
         CALL DPTENS(NS,CLDISP)
         DO I=1,6
            CDISP(I)=CDISP(I)+CLDISP(I)
         ENDDO
      ENDDO
C
      CU11=( (1.D0-BNX**2)*CKX-BNX*BNY*CKY-BNX*BNZ*CKZ)/CKPP
      CU12=(-BNY*BNX*CKX+(1.D0-BNY**2)*CKY-BNY*BNZ*CKZ)/CKPP
      CU13=(-BNZ*BNX*CKX-BNZ*BNY*CKY+(1.D0-BNZ**2)*CKZ)/CKPP
      CU21=(BNY*CKZ-BNZ*CKY)/CKPP
      CU22=(BNZ*CKX-BNX*CKZ)/CKPP
      CU23=(BNX*CKY-BNY*CKX)/CKPP
      CU31=BNX
      CU32=BNY
      CU33=BNZ
C
      CDET11=CDV11+CDISP(1)
     &            +     CU31*CU31*CDISP(2)
     &            +     CU21*CU21*CDISP(3)
     &            +2.D0*CU31*CU11*CDISP(4)
      CDET22=CDV22+CDISP(1)
     &            +     CU32*CU32*CDISP(2)
     &            +     CU22*CU22*CDISP(3)
     &            +2.D0*CU32*CU12*CDISP(4)
      CDET33=CDV33+CDISP(1)
     &            +     CU33*CU33*CDISP(2)
     &            +     CU23*CU23*CDISP(3)
     &            +2.D0*CU33*CU13*CDISP(4)
      CDET1P=           CU31*CU32*CDISP(2)
     &            +     CU21*CU22*CDISP(3)
     &            +     CU31*CU12*CDISP(4)
     &            +     CU32*CU11*CDISP(4)
      CDET2P=           CU32*CU33*CDISP(2)
     &            +     CU22*CU23*CDISP(3)
     &            +     CU32*CU13*CDISP(4)
     &            +     CU33*CU12*CDISP(4)
      CDET3P=           CU33*CU31*CDISP(2)
     &            +     CU23*CU21*CDISP(3)
     &            +     CU33*CU11*CDISP(4)
     &            +     CU31*CU13*CDISP(4)
      CDET1M=           CU33*CDISP(5)
     &            +     CU13*CDISP(6)
      CDET2M=           CU31*CDISP(5)
     &            +     CU11*CDISP(6)
      CDET3M=           CU32*CDISP(5)
     &            +     CU12*CDISP(6)
      CDET(1,1)=CDET11
      CDET(2,2)=CDET22
      CDET(3,3)=CDET33
      CDET(1,2)=CDV12+CDET1P+CDET1M
      CDET(2,1)=CDV21+CDET1P-CDET1M
      CDET(2,3)=CDV23+CDET2P+CDET2M
      CDET(3,2)=CDV32+CDET2P-CDET2M
      CDET(3,1)=CDV31+CDET3P+CDET3M
      CDET(1,3)=CDV13+CDET3P-CDET3M
C
      RETURN
      END
