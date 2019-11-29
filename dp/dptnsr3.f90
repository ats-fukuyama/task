MODULE DPTNSR3

!     ***************************************
!         COMPONENTS OF DIELECTRIC TENSOR
!             MAGNETIC FIELD     (0,   0,   B)
!             WAVE NUMBER VECTOR (k_x, 0, k_z)

!             CLDISP(1)=EPS_XX
!             CLDISP(2)=EPS_ZZ - EPS_XX
!             CLDISP(3)=EPS_YY - EPS_XX
!             CLDISP(4)=EPS_ZX
!             CLDISP(5)=EPS_XY
!             CLDISP(6)=EPS_YZ

CONTAINS

!     ****** KINETIC MODEL WITH FLR and Drift (by Swanson) ******

    SUBROUTINE DPTNKD(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE libbes,ONLY: BESEINX
      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)

      COMPLEX(rkind):: CWP,CWC,CUP,CBETA,CGZ,CZ,CDZ,CKPRL,CSIGN
      COMPLEX(rkind):: CGX,CWX,CTX,CKX,CK0,CK1,CK2,CK3,CK4,CK5
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM
      REAL(rkind):: WTPR,VTPR,WTPP,VTPP
      INTEGER:: I,NHMAX,NC,NH

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))     ! CWC is not devided by CW
      CUP=plfw(NS)%RUPR                        ! u_para

      IF(ABS(CWC).GT.1.D-16) THEN
         CSIGN=CWC/ABS(CWC)
         CWC=ABS(CWC)
      ELSE
         CSIGN=1.D0
         CWC=0.D0
      END IF


      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))  ! T_para/m
      VTPR=SQRT(2.D0*WTPR)                     ! vt_para
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))  ! T_perp/m
      VTPP=SQRT(2.D0*WTPP)                     ! vt_perp

      CBETA=CKPP*CKPP*WTPP/(CWC*CWC)           ! l.c. lambda
      DO NH=0,NHMAX
         CALAM(NH)=BESEINX(NH,DBLE(CBETA))     ! I_n(cbeta) EXP(-cbeta)
      END DO

      DO NC=NCMIN(NS),NCMAX(NS)

         IF(ABS(CKPR).LE.0.D0) THEN
            CKPRL=1.D-4
         ELSE
            CKPRL=CKPR
         END IF
         CGZ= (CW+NC*CWC-CKPR*CUP)/(ABS(CKPRL)*VTPR)
         CALL DSPFN(CGZ,CZ,CDZ)

         CGX= 1.D0-CKPR*CUP/CW
         CTX= 1.D0-WTPP/WTPR
         CWX= CW/(CKPRL*VTPR)
         CKX= CKPRL*VTPR/CW
         CSIGN=CWC/ABS(CWC)

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM      ! d Lambda/dx
         CBLAM=0.5D0*(CLAMM-CLAMP)           ! n Lambda/x
!         IF(NC.LT.0) CBLAM=-CBLAM

         CK0=-2.D0*CWX   *CBETA*CDLAM*(CGX*CZ+0.5D0*CKX*CTX*CDZ)
         CK1=      CWX*NC      *CBLAM*(CGX*CZ+0.5D0*CKX*CTX*CDZ)
         CK2=  -CI*CWX*NC*CSIGN*CDLAM*(CGX*CZ+0.5D0*CKX*CTX*CDZ)
         CK3=-     CWX*CWX*(1.D0+NC*CWC/CW)*CLAM &
                          *((1.D0+NC*CWC/CW*(1.D0-WTPR/WTPP))*CDZ &
                            +2.D0*NC*CWC*WTPR*CUP*CZ/(CW*WTPP*VTPR))
         CK4=      CKPP/CKPRL*CW/CWC*CBLAM &
                          *(NC*CWC*CUP*CZ/(CW*VTPR) &
                            +0.5D0*(WTPP/WTPR-NC*CWC/CW*CTX)*CDZ)
         CK5=  -CI*CSIGN*CKPP*CW*CDLAM/(CKPR*CWC) &
                          *(NC*CWC*CUP*CZ/(CW*VTPR) &
                            +0.5D0*(WTPP/WTPR-NC*CWC*CTX*CDZ/CW))

         CLDISP(1)=CLDISP(1)+CWP*CK1
         CLDISP(2)=CLDISP(2)+CWP*(CK3-CK1)
         CLDISP(3)=CLDISP(3)+CWP*CK0
         CLDISP(4)=CLDISP(4)+CWP*CK4
         CLDISP(5)=CLDISP(5)+CWP*CK2
         CLDISP(6)=CLDISP(6)-CWP*CK5

      ENDDO
      RETURN
    END SUBROUTINE DPTNKD
END MODULE DPTNSR3
