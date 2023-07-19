MODULE dpdisp

  PRIVATE
  PUBLIC CFDISP,CFDISPR,DP_DISP,DP_DTNS_XYZ,DP_DTNS_PZP,DPCOLD_RKPERP

CONTAINS

!     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******

  FUNCTION CFDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)

    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind):: CFDISP
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    COMPLEX(rkind):: CDET11,CDET12,CDET13,CDET21,CDET22,CDET23, &
                     CDET31,CDET32,CDET33
    COMPLEX(rkind):: CDET(3,3),CDTNS(3,3)

    SELECT CASE(MODEL_ES)
    CASE(0)
       CALL DP_DISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)

       CDET11=CDET(1,1)
       CDET12=CDET(1,2)
       CDET13=CDET(1,3)
       CDET21=CDET(2,1)
       CDET22=CDET(2,2)
       CDET23=CDET(2,3)
       CDET31=CDET(3,1)
       CDET32=CDET(3,2)
       CDET33=CDET(3,3)

       CFDISP =CDET33*(CDET11*CDET22-CDET12*CDET21) &
              +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11 &
              +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
    CASE(1)
       CALL DP_DTNS_XYZ(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,0,CDTNS)
!       RK=SQRT(CKX*CONJG(CKX)+CKY*CONJG(CKY)+CKX*CONJG(CKY))
!       CKXN=CKX/RK
!       CKYN=CKY/RK
!       CKZN=CKZ/RK
       CFDISP=CKX*(CDTNS(1,1)*CKX+CDTNS(1,2)*CKY+CDTNS(1,3)*CKZ) &
             +CKY*(CDTNS(2,1)*CKX+CDTNS(2,2)*CKY+CDTNS(2,3)*CKZ) &
             +CKZ*(CDTNS(3,1)*CKX+CDTNS(3,2)*CKY+CDTNS(3,3)*CKZ)
    END SELECT
    RETURN
  END FUNCTION CFDISP

!     ****** CALCULATE REAL PART OF DETERMINANT OF DISPERSION TENSOR ******


  FUNCTION CFDISPR(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind):: CFDISPR
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    COMPLEX(rkind):: CDET11,CDET12,CDET13,CDET21,CDET22,CDET23, &
                     CDET31,CDET32,CDET33
    COMPLEX(rkind):: CDTNS11,CDTNS12,CDTNS13,CDTNS21,CDTNS22,CDTNS23, &
                     CDTNS31,CDTNS32,CDTNS33
    COMPLEX(rkind):: CDET(3,3),CDTNS(3,3)
    COMPLEX(rkind):: CKXN,CKYN,CKZN
    REAL(rkind):: RK

    SELECT CASE(MODEL_ES)
    CASE(0)
    
       CALL DP_DISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)

       CDET11=0.5D0*(CDET(1,1)+DCONJG(CDET(1,1)))
       CDET12=0.5D0*(CDET(1,2)+DCONJG(CDET(2,1)))
       CDET13=0.5D0*(CDET(1,3)+DCONJG(CDET(3,1)))
       CDET21=0.5D0*(CDET(2,1)+DCONJG(CDET(1,2)))
       CDET22=0.5D0*(CDET(2,2)+DCONJG(CDET(2,2)))
       CDET23=0.5D0*(CDET(2,3)+DCONJG(CDET(3,2)))
       CDET31=0.5D0*(CDET(3,1)+DCONJG(CDET(1,3)))
       CDET32=0.5D0*(CDET(3,2)+DCONJG(CDET(2,3)))
       CDET33=0.5D0*(CDET(3,3)+DCONJG(CDET(3,3)))

       CFDISPR=CDET33*(CDET11*CDET22-CDET12*CDET21) &
              +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11 &
              +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
    CASE(1)
       CALL DP_DTNS_XYZ(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,0,CDTNS)
       CDTNS11=0.5D0*(CDTNS(1,1)+DCONJG(CDTNS(1,1)))
       CDTNS12=0.5D0*(CDTNS(1,2)+DCONJG(CDTNS(2,1)))
       CDTNS13=0.5D0*(CDTNS(1,3)+DCONJG(CDTNS(3,1)))
       CDTNS21=0.5D0*(CDTNS(2,1)+DCONJG(CDTNS(1,2)))
       CDTNS22=0.5D0*(CDTNS(2,2)+DCONJG(CDTNS(2,2)))
       CDTNS23=0.5D0*(CDTNS(2,3)+DCONJG(CDTNS(3,2)))
       CDTNS31=0.5D0*(CDTNS(3,1)+DCONJG(CDTNS(1,3)))
       CDTNS32=0.5D0*(CDTNS(3,2)+DCONJG(CDTNS(2,3)))
       CDTNS33=0.5D0*(CDTNS(3,3)+DCONJG(CDTNS(3,3)))
       RK=SQRT(CKX*CONJG(CKX)+CKY*CONJG(CKY)+CKX*CONJG(CKY))
       CKXN=CKX/RK
       CKYN=CKY/RK
       CKZN=CKZ/RK

       CFDISPR=CKXN*(CDTNS11*CKXN+CDTNS12*CKYN+CDTNS13*CKZN) &
              +CKYN*(CDTNS21*CKXN+CDTNS22*CKYN+CDTNS23*CKZN) &
              +CKZN*(CDTNS31*CKXN+CDTNS32*CKYN+CDTNS33*CKZN)
    END SELECT
       

    RETURN
  END FUNCTION CFDISPR

!     ****** CALCULATE DISPERSION TENSOR ******

  SUBROUTINE DP_DISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)

    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    COMPLEX(rkind),INTENT(OUT):: CDET(3,3)
    COMPLEX(rkind):: CDTNS(3,3)
    COMPLEX(rkind):: CW,CL2

    CW=2.D0*PI*1.D6*CRF
    IF(ABS(CW).LE.1.D-8) CW=(1.D-8,0.D0)
    CL2=VC*VC/(CW*CW)

    CALL DP_DTNS_XYZ(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,0,CDTNS)

    CDET(1,1)=-CL2*(CKY*CKY+CKZ*CKZ)+CDTNS(1,1)
    CDET(1,2)= CL2* CKX*CKY         +CDTNS(1,2)
    CDET(1,3)= CL2* CKX*CKZ         +CDTNS(1,3)
    CDET(2,1)= CL2* CKY*CKX         +CDTNS(2,1)
    CDET(2,2)=-CL2*(CKZ*CKZ+CKX*CKX)+CDTNS(2,2)
    CDET(2,3)= CL2* CKY*CKZ         +CDTNS(2,3)
    CDET(3,1)= CL2* CKZ*CKX         +CDTNS(3,1)
    CDET(3,2)= CL2* CKZ*CKY         +CDTNS(3,2)
    CDET(3,3)=-CL2*(CKX*CKX+CKY*CKY)+CDTNS(3,3)

    RETURN
  END SUBROUTINE DP_DISP

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE DP_DTNS_XYZ(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,NS,CDTNS)

    USE dpcomm
    USE plprof
    USE plprofw
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    INTEGER,INTENT(IN):: NS
    COMPLEX(rkind),INTENT(OUT):: CDTNS(3,3)
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd
    REAL(rkind):: RHON,BNXY
    COMPLEX(rkind):: CW,CKPR,CKPP
    COMPLEX(rkind):: CU11,CU12,CU13,CU21,CU22,CU23,CU31,CU32,CU33
    COMPLEX(rkind):: CDET11,CDET12,CDET13, &
                     CDET21,CDET22,CDET23, &
                     CDET31,CDET32,CDET33

    CALL PL_MAG(XPOS,YPOS,ZPOS,mag)
    RHON=mag%rhon
    IF(MODELG.LE.1.OR.MODELG.GT.10) THEN
       CALL PL_PROFW3D(XPOS,YPOS,ZPOS,plfw)
    ELSE
       CALL PL_PROFW(RHON,plfw)
    END IF
    CALL PL_GRAD(RHON,grd)

    CW=2.D0*PI*1.D6*CRF
    IF(ABS(CW).LE.1.D-8) CW=(1.D-8,0.D0)
    CKPR=mag%BNX*CKX+mag%BNY*CKY+mag%BNZ*CKZ
    IF(ABS(CKPR).LE.1.D-8) CKPR=1.D-8
    CKPP=SQRT(CKX**2+CKY**2+CKZ**2-CKPR**2)

    CALL DP_DTNS_PZP(CW,CKPR,CKPP,NS,mag,plfw,CDTNS,GRD=grd)

    IF(ABS(CKPP).EQ.0.D0) THEN
       BNXY=SQRT(mag%BNX*mag%BNX+mag%BNY*mag%BNY)
       IF(BNXY.EQ.0.D0) THEN
          CU11= 1.D0
          CU12= 0.D0
          CU13= 0.D0
          CU21= 0.D0
          CU22= 1.D0
          CU23= 0.D0
          CU31= 0.D0
          CU32= 0.D0
          CU33= 1.D0
       ELSE
          CU11=-mag%BNZ*mag%BNX/BNXY
          CU12=-mag%BNZ*mag%BNY/BNXY
          CU13= BNXY
          CU21= mag%BNY/BNXY
          CU22=-mag%BNX/BNXY
          CU23= 0.D0
          CU31= mag%BNX
          CU32= mag%BNY
          CU33= mag%BNZ
       ENDIF
    ELSE
       CU11=( (1.D0-mag%BNX**2)*CKX &
             -mag%BNX*mag%BNY*CKY &
             -mag%BNX*mag%BNZ*CKZ )/CKPP
       CU12=(-mag%BNY*mag%BNX*CKX &
             +(1.D0-mag%BNY**2)*CKY &
             -mag%BNY*mag%BNZ*CKZ )/CKPP
       CU13=(-mag%BNZ*mag%BNX*CKX &
             -mag%BNZ*mag%BNY*CKY &
             +(1.D0-mag%BNZ**2)*CKZ )/CKPP
       CU21=(mag%BNY*CKZ-mag%BNZ*CKY)/CKPP
       CU22=(mag%BNZ*CKX-mag%BNX*CKZ)/CKPP
       CU23=(mag%BNX*CKY-mag%BNY*CKX)/CKPP
       CU31=mag%BNX
       CU32=mag%BNY
       CU33=mag%BNZ
    ENDIF

    CDET11=CDTNS(1,1)*CU11+CDTNS(1,2)*CU21+CDTNS(1,3)*CU31
    CDET12=CDTNS(1,1)*CU12+CDTNS(1,2)*CU22+CDTNS(1,3)*CU32
    CDET13=CDTNS(1,1)*CU13+CDTNS(1,2)*CU23+CDTNS(1,3)*CU33
    CDET21=CDTNS(2,1)*CU11+CDTNS(2,2)*CU21+CDTNS(2,3)*CU31
    CDET22=CDTNS(2,1)*CU12+CDTNS(2,2)*CU22+CDTNS(2,3)*CU32
    CDET23=CDTNS(2,1)*CU13+CDTNS(2,2)*CU23+CDTNS(2,3)*CU33
    CDET31=CDTNS(3,1)*CU11+CDTNS(3,2)*CU21+CDTNS(3,3)*CU31
    CDET32=CDTNS(3,1)*CU12+CDTNS(3,2)*CU22+CDTNS(3,3)*CU32
    CDET33=CDTNS(3,1)*CU13+CDTNS(3,2)*CU23+CDTNS(3,3)*CU33

    CDTNS(1,1)=CU11*CDET11+CU21*CDET21+CU31*CDET31
    CDTNS(1,2)=CU11*CDET12+CU21*CDET22+CU31*CDET32
    CDTNS(1,3)=CU11*CDET13+CU21*CDET23+CU31*CDET33
    CDTNS(2,1)=CU12*CDET11+CU22*CDET21+CU32*CDET31
    CDTNS(2,2)=CU12*CDET12+CU22*CDET22+CU32*CDET32
    CDTNS(2,3)=CU12*CDET13+CU22*CDET23+CU32*CDET33
    CDTNS(3,1)=CU13*CDET11+CU23*CDET21+CU33*CDET31
    CDTNS(3,2)=CU13*CDET12+CU23*CDET22+CU33*CDET32
    CDTNS(3,3)=CU13*CDET13+CU23*CDET23+CU33*CDET33

    RETURN
  END SUBROUTINE DP_DTNS_XYZ

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE DP_DTNS_PZP(CW,CKPR,CKPP,NS,mag,plfw,CDTNS,grd)

    USE dpcomm
    USE plprof
    USE plprofw
    USE dptnsr0
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax),OPTIONAL:: grd
    COMPLEX(rkind),INTENT(OUT):: CDTNS(3,3)
    COMPLEX(rkind):: CDISP(6),CLDISP(6)
    INTEGER:: I,NS1

    IF(NS.EQ.0) THEN
       CDISP(1)=1.D0
       DO I=2,6
          CDISP(I)=0.D0
       ENDDO

       IF(PRESENT(grd)) THEN
          DO ns1=1,nsmax
             grd(ns1)%grdn=0.D0
             grd(ns1)%grdtpr=0.D0
             grd(ns1)%grdtpp=0.D0
             grd(ns1)%grdu=0.D0
          END DO
       END IF
          
       DO NS1=1,NSMAX

          CALL DP_TNSR0(CW,CKPR,CKPP,NS1,mag,plfw,grd,CLDISP)
          DO I=1,6
             CDISP(I)=CDISP(I)+CLDISP(I)
          ENDDO
       ENDDO
    ELSE
       IF(PRESENT(grd)) THEN
          grd(ns)%grdn=0.D0
          grd(ns)%grdtpr=0.D0
          grd(ns)%grdtpp=0.D0
          grd(ns)%grdu=0.D0
       END IF

       CALL DP_TNSR0(CW,CKPR,CKPP,NS,mag,plfw,grd,CDISP)
    ENDIF

    CDTNS(1,1)= CDISP(1)
    CDTNS(1,2)= CDISP(5)
    CDTNS(1,3)= CDISP(4)
    CDTNS(2,1)=-CDISP(5)
    CDTNS(2,2)= CDISP(1)+CDISP(3)
    CDTNS(2,3)= CDISP(6)
    CDTNS(3,1)= CDISP(4)
    CDTNS(3,2)=-CDISP(6)
    CDTNS(3,3)= CDISP(1)+CDISP(2)

    RETURN
  END SUBROUTINE DP_DTNS_PZP

!     ****** CALCULATE COLD KPERP ******

  SUBROUTINE DPCOLD_RKPERP(CW,CKPR,mag,plfw,CKPPF,CKPPS)

    USE dpcomm
    USE plprof
    USE plprofw
    USE dptnsr0
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR
    TYPE(pl_mag_type),INTENT(IN):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd
    COMPLEX(rkind),INTENT(OUT):: CKPPF,CKPPS
    COMPLEX(rkind):: CDISP(6),CLDISP(6)
    REAL(rkind):: RKPPA2,RKPPB2
    COMPLEX(rkind):: CKPP,CNPR,CCS,CCD,CCP,CCA,CCB,CCC,CKPPA,CKPPB,CTEMP
    INTEGER:: I,NS1,MODELP_SAVE

    DO ns1=1,nsmax
       grd(ns1)%grdn=0.D0
       grd(ns1)%grdtpr=0.D0
       grd(ns1)%grdtpp=0.D0
       grd(ns1)%grdu=0.D0
    END DO

    CDISP(1)=1.D0
    DO I=2,6
       CDISP(I)=0.D0
    ENDDO
    CKPP=(0.D0,0.D0)
    CNPR=CKPR*VC/CW

    DO NS1=1,NSMAX
       MODELP_SAVE=MODELP(NS1)
       MODELP(NS1)=0
       CALL DP_TNSR0(CW,CKPR,CKPP,NS1,mag,plfw,grd,CLDISP)
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

    RKPPA2=REAL(CKPPA)
    RKPPB2=REAL(CKPPB)
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
  END SUBROUTINE DPCOLD_RKPERP
END MODULE dpdisp
