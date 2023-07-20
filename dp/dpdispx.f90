MODULE dpdispx

  PRIVATE
  PUBLIC cf_dispx
  PUBLIC dp_dtns

CONTAINS

!     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******

  FUNCTION CF_DISPX(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,NS_,MODE_ES_)

    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind):: CF_DISPX
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    INTEGER,OPTIONAL:: NS_  ! Particles species
                            ! NS=0  : H(all)+AH(all)
                            ! NS!=0 : H(all)+AH(NS)
    INTEGER,OPTIONAL:: MODE_ES_
    INTEGER:: NS,MODE_ES
    COMPLEX(rkind):: CW,CL2,CDET(3,3),CDTNS(3,3)
    INTEGER:: I,J

    IF(PRESENT(NS_)) THEN
       NS=NS_
    ELSE
       NS=0
    END IF
    IF(PRESENT(MODE_ES_)) THEN
       MODE_ES=MODE_ES_
    ELSE
       MODE_ES=0
    END IF
       
    CALL DP_DTNS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,NS)
    
    SELECT CASE(MODEL_ES)
    CASE(0)
       CW=2.D0*PI*1.D6*CRF
       IF(ABS(CW).LE.1.D-8) CW=(1.D-8,0.D0)
       CL2=VC*VC/(CW*CW)
       CDET(1,1)=-CL2*(CKY*CKY+CKZ*CKZ)+CDTNS(1,1)
       CDET(1,2)= CL2* CKX*CKY         +CDTNS(1,2)
       CDET(1,3)= CL2* CKX*CKZ         +CDTNS(1,3)
       CDET(2,1)= CL2* CKY*CKX         +CDTNS(2,1)
       CDET(2,2)=-CL2*(CKZ*CKZ+CKX*CKX)+CDTNS(2,2)
       CDET(2,3)= CL2* CKY*CKZ         +CDTNS(2,3)
       CDET(3,1)= CL2* CKZ*CKX         +CDTNS(3,1)
       CDET(3,2)= CL2* CKZ*CKY         +CDTNS(3,2)
       CDET(3,3)=-CL2*(CKX*CKX+CKY*CKY)+CDTNS(3,3)

       CF_DISPX=CDET(1,1)*CDET(2,2)*CDET(3,3)-CDET(1,1)*CDET(2,3)*CDET(3,2) &
               +CDET(1,2)*CDET(2,3)*CDET(3,1)-CDET(1,2)*CDET(2,1)*CDET(3,3) &
               +CDET(1,3)*CDET(2,1)*CDET(3,2)-CDET(1,3)*CDET(2,1)*CDET(3,1)
    CASE(1)
       CF_DISPX=CKX*(CDTNS(1,1)*CKX+CDTNS(1,2)*CKY+CDTNS(1,3)*CKZ) &
               +CKY*(CDTNS(2,1)*CKX+CDTNS(2,2)*CKY+CDTNS(2,3)*CKZ) &
               +CKZ*(CDTNS(3,1)*CKX+CDTNS(3,2)*CKY+CDTNS(3,3)*CKZ)
    END SELECT
    RETURN
  END FUNCTION CF_DISPX

  !     ****** CALCULATE DIELECTRIC TENSOR :

  SUBROUTINE DP_DTNS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,NS_)

    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    COMPLEX(rkind),INTENT(OUT):: CDTNS(3,3)
    INTEGER,OPTIONAL:: NS_  ! Particles species
                            ! NS=0  : H(all)+AH(all)
                            ! NS!=0 : H(all)+AH(NS)
    COMPLEX(rkind):: CDTNS_H(3,3),CDTNS_AH(3,3)
    REAL(rkind):: RK
    INTEGER:: I,J,NS

    IF(PRESENT(NS_)) THEN
       NS=NS_
    ELSE
       NS=0
    END IF

    IF(NS.EQ.0) THEN ! All plasma component
       CALL DP_DTNS_XYZ_NS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,0)
    ELSE ! Single Anti-Hermite component
       CALL DP_DTNS_XYZ_NS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,0)
       DO J=1,3
          DO I=1,3
             CDTNS_H(I,J)=0.5D0*(CDTNS(I,J)+DCONJG(CDTNS(J,I))) ! All Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS0_H:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(1,1),CDTNS_H(1,2),CDTNS_H(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(2,1),CDTNS_H(2,2),CDTNS_H(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(3,1),CDTNS_H(3,2),CDTNS_H(3,3)
       DO J=1,3
          DO I=1,3
             CDTNS_AH(I,J)=0.5D0*(CDTNS(I,J)-DCONJG(CDTNS(J,I))) ! Anti-Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS0_AH:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(1,1),CDTNS_AH(1,2),CDTNS_AH(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(2,1),CDTNS_AH(2,2),CDTNS_AH(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(3,1),CDTNS_AH(3,2),CDTNS_AH(3,3)

       CALL DP_DTNS_XYZ_NS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,1)
       DO J=1,3
          DO I=1,3
             CDTNS_H(I,J)=0.5D0*(CDTNS(I,J)+DCONJG(CDTNS(J,I))) ! All Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS1_H:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(1,1),CDTNS_H(1,2),CDTNS_H(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(2,1),CDTNS_H(2,2),CDTNS_H(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(3,1),CDTNS_H(3,2),CDTNS_H(3,3)
       DO J=1,3
          DO I=1,3
             CDTNS_AH(I,J)=0.5D0*(CDTNS(I,J)-DCONJG(CDTNS(J,I))) ! Anti-Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS1_AH:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(1,1),CDTNS_AH(1,2),CDTNS_AH(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(2,1),CDTNS_AH(2,2),CDTNS_AH(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(3,1),CDTNS_AH(3,2),CDTNS_AH(3,3)

       CALL DP_DTNS_XYZ_NS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,2)
       DO J=1,3
          DO I=1,3
             CDTNS_H(I,J)=0.5D0*(CDTNS(I,J)+DCONJG(CDTNS(J,I))) ! All Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS1_H:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(1,1),CDTNS_H(1,2),CDTNS_H(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(2,1),CDTNS_H(2,2),CDTNS_H(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_H(3,1),CDTNS_H(3,2),CDTNS_H(3,3)
       DO J=1,3
          DO I=1,3
             CDTNS_AH(I,J)=0.5D0*(CDTNS(I,J)-DCONJG(CDTNS(J,I))) ! Anti-Hermite
          END DO
       END DO
       WRITE(6,'(A)'        ) '## DP_DTNS1_AH:'
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(1,1),CDTNS_AH(1,2),CDTNS_AH(1,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(2,1),CDTNS_AH(2,2),CDTNS_AH(2,3)
       WRITE(6,'(A,6ES12.4)') '   ',CDTNS_AH(3,1),CDTNS_AH(3,2),CDTNS_AH(3,3)
    END IF
    RETURN
  END SUBROUTINE DP_DTNS

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE DP_DTNS_XYZ_NS(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDTNS,NS)

    USE dpcomm
    USE plprof
    USE plprofw
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CRF,CKX,CKY,CKZ
    REAL(rkind),INTENT(IN):: XPOS,YPOS,ZPOS
    COMPLEX(rkind),INTENT(OUT):: CDTNS(3,3)
    INTEGER,OPTIONAL:: NS   ! Particles species
                            ! NS=0  : all
                            ! NS!=0 : only NS
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd
    REAL(rkind):: RHON,BNXY
    COMPLEX(rkind):: CW,CKPR,CKPP
    COMPLEX(rkind):: CU(3,3),CDTNS1(3,3)
    INTEGER:: I,J,NS1

    IF(.NOT.PRESENT(NS)) NS=0

    CALL PL_MAG(XPOS,YPOS,ZPOS,mag)
    RHON=mag%rhon
    IF(MODELG.LE.1.OR.MODELG.GT.10) THEN
       CALL PL_PROFW3D(XPOS,YPOS,ZPOS,plfw)
    ELSE
       CALL PL_PROFW(RHON,plfw)
    END IF
    IF(MODELG.LE.1.OR.MODELG.GT.10) THEN ! not supported
       DO ns1=1,nsmax
          grd(ns1)%grdn=0.D0
          grd(ns1)%grdtpr=0.D0
          grd(ns1)%grdtpp=0.D0
          grd(ns1)%grdu=0.D0
       END DO
    ELSE
       CALL PL_GRAD(RHON,grd)
    END IF

    CW=2.D0*PI*1.D6*CRF
    IF(ABS(CW).LE.1.D-8) CW=(1.D-8,0.D0)
    CKPR=mag%BNX*CKX+mag%BNY*CKY+mag%BNZ*CKZ
    IF(ABS(CKPR).LE.1.D-8) CKPR=1.D-8
    CKPP=SQRT(CKX**2+CKY**2+CKZ**2-CKPR**2)

    CALL DP_DTNS_123_NS(CW,CKPR,CKPP,NS,mag,plfw,CDTNS,GRD=grd)

    IF(ABS(CKPP).EQ.0.D0) THEN
       BNXY=SQRT(mag%BNX*mag%BNX+mag%BNY*mag%BNY)
       IF(BNXY.EQ.0.D0) THEN
          CU(1,1)= 1.D0
          CU(1,2)= 0.D0
          CU(1,3)= 0.D0
          CU(2,1)= 0.D0
          CU(2,2)= 1.D0
          CU(2,3)= 0.D0
          CU(3,1)= 0.D0
          CU(3,2)= 0.D0
          CU(3,3)= 1.D0
       ELSE
          CU(1,1)=-mag%BNZ*mag%BNX/BNXY
          CU(1,2)=-mag%BNZ*mag%BNY/BNXY
          CU(1,3)= BNXY
          CU(2,1)= mag%BNY/BNXY
          CU(2,2)=-mag%BNX/BNXY
          CU(2,3)= 0.D0
          CU(3,1)= mag%BNX
          CU(3,2)= mag%BNY
          CU(3,3)= mag%BNZ
       ENDIF
    ELSE
       CU(1,1)=( (1.D0-mag%BNX**2)    *CKX &
                      -mag%BNX*mag%BNY*CKY &
                      -mag%BNX*mag%BNZ*CKZ )/CKPP
       CU(1,2)=(      -mag%BNY*mag%BNX*CKX &
                +(1.D0-mag%BNY**2)    *CKY &
                      -mag%BNY*mag%BNZ*CKZ )/CKPP
       CU(1,3)=(      -mag%BNZ*mag%BNX*CKX &
                      -mag%BNZ*mag%BNY*CKY &
                +(1.D0-mag%BNZ**2)    *CKZ )/CKPP
       CU(2,1)=       (mag%BNY*CKZ-mag%BNZ*CKY)/CKPP
       CU(2,2)=       (mag%BNZ*CKX-mag%BNX*CKZ)/CKPP
       CU(2,3)=       (mag%BNX*CKY-mag%BNY*CKX)/CKPP
       CU(3,1)=        mag%BNX
       CU(3,2)=        mag%BNY
       CU(3,3)=        mag%BNZ
    ENDIF

    DO J=1,3
       DO I=1,3
          CDTNS1(I,J)=CDTNS(I,1)*CU(1,J) &
                     +CDTNS(I,2)*CU(2,J) &
                     +CDTNS(I,3)*CU(3,J)
       END DO
    END DO
    DO J=1,3
       DO I=1,3
          CDTNS(I,J)=CU(I,1)*CDTNS1(1,J) &
                    +CU(I,2)*CDTNS1(2,J) &
                    +CU(I,3)*CDTNS1(3,J)
       END DO
    END DO

    RETURN
  END SUBROUTINE DP_DTNS_XYZ_NS

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE DP_DTNS_123_NS(CW,CKPR,CKPP,NS,mag,plfw,CDTNS,grd)

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

       IF(.NOT.PRESENT(grd)) THEN
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
  END SUBROUTINE DP_DTNS_123_NS

END MODULE dpdispx
